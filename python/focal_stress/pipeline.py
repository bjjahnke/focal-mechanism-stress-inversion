"""Run all three steps from one config file and write every output to disk."""

from __future__ import annotations

import json
import logging
import os
import platform
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .catalog import load_catalog
from .iterative_inversion import run_stress_inversions
from .nodal_planes import compute_nodal_planes
from .summarize import summarize_inversions

log = logging.getLogger(__name__)

OUTPUT_FILES = {
    "focal_mechanisms": "focal_mechanisms.csv",
    "inversion_runs": "inversion_runs.csv",
    "stress_solutions": "stress_solutions.csv",
    "solution_planes": "solution_planes.csv",
    "nodal_plane_frequency": "nodal_plane_frequency.csv",
    "metadata": "run_metadata.json",
}


@dataclass
class Config:
    input_catalog: Path
    output_dir: Path
    friction: float = 0.6
    n_runs: int = 1000
    max_iterations: int = 20
    random_seed: int | None = None
    significance_level: float = 0.05
    make_plots: bool = True
    extra: dict = field(default_factory=dict)

    @classmethod
    def from_file(cls, path) -> "Config":
        """Load a JSON config. Relative paths are resolved from the config file's folder."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")
        try:
            raw = json.loads(path.read_text())
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path}: not valid JSON ({exc})") from exc
        for key in ("input_catalog", "output_dir"):
            if key not in raw:
                raise ValueError(f"{path}: config is missing required key '{key}'")
        known = set(cls.__dataclass_fields__) - {"extra"}
        unknown = {k: v for k, v in raw.items() if k not in known}
        if unknown:
            log.warning("Ignoring unknown config keys: %s", sorted(unknown))
        base = path.parent
        kwargs = {k: v for k, v in raw.items() if k in known}
        kwargs["input_catalog"] = (base / raw["input_catalog"]).resolve()
        kwargs["output_dir"] = (base / raw["output_dir"]).resolve()
        return cls(**kwargs, extra=unknown)


def run_pipeline(config: Config) -> dict[str, Path]:
    """Steps 1-3 plus plots. Returns the paths of every file written."""
    started = time.perf_counter()
    out = Path(config.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = {}

    log.info("Step 1/3: computing nodal planes")
    catalog = load_catalog(config.input_catalog)
    mechanisms = compute_nodal_planes(catalog)
    written["focal_mechanisms"] = _write(mechanisms, out, "focal_mechanisms")

    log.info("Step 2/3: %d iterative inversions (friction %.2f)", config.n_runs, config.friction)
    runs = run_stress_inversions(mechanisms, config.friction, config.n_runs,
                                 config.max_iterations, config.random_seed)
    written["inversion_runs"] = _write(runs, out, "inversion_runs")

    log.info("Step 3/3: summarizing and validating")
    summary = summarize_inversions(mechanisms, runs, config.friction, config.significance_level)
    for name in ("stress_solutions", "solution_planes", "nodal_plane_frequency"):
        written[name] = _write(getattr(summary, name), out, name)

    if config.make_plots:
        from .plots import make_plots
        for p in make_plots(summary, out):
            written[f"plot_{p.stem}"] = p

    metadata = {
        "package_version": __version__,
        "run_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_seconds": round(time.perf_counter() - started, 2),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        # Paths are stored relative to the output folder so metadata is portable.
        "config": {k: os.path.relpath(v, out) if isinstance(v, Path) else v
                   for k, v in config.__dict__.items() if k != "extra"},
        "n_events": len(mechanisms),
        "run_status_counts": {k: int(v) for k, v in runs["status"].value_counts().items()},
        "n_stress_solutions": len(summary.stress_solutions),
        "outputs": {k: str(v.relative_to(out)) for k, v in written.items()},
    }
    meta_path = out / OUTPUT_FILES["metadata"]
    meta_path.write_text(json.dumps(metadata, indent=2))
    written["metadata"] = meta_path

    top = summary.stress_solutions.iloc[0]
    log.info("Done in %.1fs. Top solution reached by %.0f%% of runs: sigma1 %.0f/%.0f, R = %.2f",
             metadata["runtime_seconds"], 100 * top["fraction_of_runs"],
             top["sigma1_trend"], top["sigma1_plunge"], top["shape_ratio"])
    return written


def _write(df: pd.DataFrame, out: Path, name: str) -> Path:
    path = out / OUTPUT_FILES[name]
    df.to_csv(path, index=False)
    log.info("  wrote %s (%d rows)", path.name, len(df))
    return path

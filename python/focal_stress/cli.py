"""Command-line interface.

    focal-stress run CONFIG.json
    focal-stress nodal-planes CATALOG.csv OUTPUT.csv
    focal-stress invert FOCAL_MECHANISMS.csv OUTPUT.csv [--friction --n-runs --max-iterations --seed]
    focal-stress summarize FOCAL_MECHANISMS.csv INVERSION_RUNS.csv OUTPUT_DIR [--alpha --no-plots]
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from .catalog import CatalogError, load_catalog
from .iterative_inversion import run_stress_inversions
from .nodal_planes import compute_nodal_planes
from .pipeline import Config, run_pipeline
from .summarize import summarize_inversions


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="focal-stress", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-q", "--quiet", action="store_true", help="only print warnings and errors")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("run", help="run all three steps from a config file")
    p.add_argument("config", type=Path)

    p = sub.add_parser("nodal-planes", help="step 1: add the second nodal plane and P/T/N axes")
    p.add_argument("catalog", type=Path)
    p.add_argument("output", type=Path)

    p = sub.add_parser("invert", help="step 2: repeated iterative stress inversions")
    p.add_argument("mechanisms", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--friction", type=float, default=0.6)
    p.add_argument("--n-runs", type=int, default=1000)
    p.add_argument("--max-iterations", type=int, default=20)
    p.add_argument("--seed", type=int, default=None)

    p = sub.add_parser("summarize", help="step 3: group runs into stress solutions and validate")
    p.add_argument("mechanisms", type=Path)
    p.add_argument("runs", type=Path)
    p.add_argument("output_dir", type=Path)
    p.add_argument("--alpha", type=float, default=0.05, help="significance level for the variance test")
    p.add_argument("--no-plots", action="store_true")

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                        format="%(levelname)s %(message)s")
    try:
        return _dispatch(args)
    except (CatalogError, ValueError, FileNotFoundError) as exc:
        logging.error("%s", exc)
        return 1


def _dispatch(args) -> int:
    if args.command == "run":
        run_pipeline(Config.from_file(args.config))
    elif args.command == "nodal-planes":
        _save(compute_nodal_planes(load_catalog(args.catalog)), args.output)
    elif args.command == "invert":
        mech = pd.read_csv(args.mechanisms, dtype={"event_id": str})
        runs = run_stress_inversions(mech, args.friction, args.n_runs, args.max_iterations, args.seed)
        _save(runs, args.output)
    elif args.command == "summarize":
        mech = pd.read_csv(args.mechanisms, dtype={"event_id": str})
        runs = pd.read_csv(args.runs, dtype={"plane_choice": str})
        summary = summarize_inversions(mech, runs, alpha=args.alpha)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name in ("stress_solutions", "solution_planes", "nodal_plane_frequency"):
            _save(getattr(summary, name), args.output_dir / f"{name}.csv")
        if not args.no_plots:
            from .plots import make_plots
            make_plots(summary, args.output_dir)
    return 0


def _save(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    logging.info("wrote %s (%d rows)", path, len(df))


if __name__ == "__main__":
    sys.exit(main())

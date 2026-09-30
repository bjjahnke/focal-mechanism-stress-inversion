"""Load and validate focal mechanism catalogs."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

REQUIRED_COLUMNS = ("strike", "dip", "rake")
OPTIONAL_COLUMNS = ("event_id", "uncertainty_deg")
VALID_RANGES = {
    "strike": (0.0, 360.0),
    "dip": (0.0, 90.0),
    "rake": (-180.0, 180.0),
}
MIN_EVENTS_FOR_INVERSION = 3
RECOMMENDED_EVENTS = 10


class CatalogError(ValueError):
    """Raised when an input catalog cannot be used."""


def validate_catalog(catalog: pd.DataFrame, source: str = "catalog", min_events: int = 1) -> pd.DataFrame:
    """Check a catalog and return a cleaned copy.

    * ``strike``, ``dip``, ``rake`` must exist, be numeric, have no gaps and
      sit inside their valid ranges.
    * ``event_id`` is created (EV001, EV002, ...) if missing and must be unique.
    * ``uncertainty_deg``, if present, must be positive with no gaps.
    * Any other columns are kept unchanged.
    """
    df = catalog.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise CatalogError(
            f"{source}: missing required column(s) {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    problems = []
    for col in REQUIRED_COLUMNS:
        values = pd.to_numeric(df[col], errors="coerce")
        bad = values.isna()
        if bad.any():
            problems.append(f"'{col}' has missing or non-numeric values in rows {_rows(bad)}")
        lo, hi = VALID_RANGES[col]
        out = (values < lo) | (values > hi)
        if out.any():
            problems.append(f"'{col}' must be between {lo:g} and {hi:g}; rows {_rows(out)} are not")
        df[col] = values

    if "uncertainty_deg" in df.columns:
        unc = pd.to_numeric(df["uncertainty_deg"], errors="coerce")
        bad = unc.isna() | (unc <= 0)
        if bad.any():
            problems.append(f"'uncertainty_deg' must be a positive number; rows {_rows(bad)} are not")
        df["uncertainty_deg"] = unc

    if problems:
        raise CatalogError(f"{source}: " + "; ".join(problems))

    if len(df) < min_events:
        raise CatalogError(f"{source}: need at least {min_events} events, found {len(df)}")

    if "event_id" not in df.columns:
        df.insert(0, "event_id", [f"EV{i + 1:03d}" for i in range(len(df))])
    df["event_id"] = df["event_id"].astype(str)
    dupes = df["event_id"].duplicated()
    if dupes.any():
        raise CatalogError(f"{source}: duplicate event_id values {sorted(set(df.loc[dupes, 'event_id']))}")

    return df.reset_index(drop=True)


def load_catalog(path: str | Path) -> pd.DataFrame:
    """Read a CSV catalog from disk and validate it."""
    path = Path(path)
    if not path.exists():
        raise CatalogError(f"Catalog file not found: {path}")
    df = pd.read_csv(path, dtype={"event_id": str})
    log.info("Loaded %d events from %s", len(df), path)
    return validate_catalog(df, source=str(path))


def _rows(mask) -> str:
    idx = np.flatnonzero(np.asarray(mask)) + 2  # +2: header line and 1-based rows
    shown = ", ".join(str(i) for i in idx[:10])
    return shown + (" ..." if len(idx) > 10 else "")

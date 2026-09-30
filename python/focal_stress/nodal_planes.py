"""Step 1: compute the second nodal plane and P, T, N axes.

A focal mechanism has two possible fault planes (nodal planes). Catalogs
usually list only one. Given that plane's strike, dip and rake, the second
plane is fixed by geometry: its normal is the first plane's slip direction and
its slip direction is the first plane's normal.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .catalog import validate_catalog
from .geometry import plane_vectors, to_trend_plunge, vectors_to_plane

NODAL_PLANE_COLUMNS = [
    "strike_1", "dip_1", "rake_1",
    "strike_2", "dip_2", "rake_2",
    "p_trend", "p_plunge", "t_trend", "t_plunge", "n_trend", "n_plunge",
]


def compute_nodal_planes(catalog: pd.DataFrame, decimals: int = 2) -> pd.DataFrame:
    """Add both nodal planes and the P, T, N axes to a catalog.

    Parameters
    ----------
    catalog : DataFrame with at least ``strike``, ``dip``, ``rake``.
    decimals : rounding for the new angle columns.

    Returns
    -------
    DataFrame with ``event_id``, the new columns in ``NODAL_PLANE_COLUMNS``,
    ``uncertainty_deg`` if it was supplied, then every other input column.
    """
    cat = validate_catalog(catalog)
    normal_1, slip_1 = plane_vectors(cat["strike"], cat["dip"], cat["rake"])

    # Plane 2 swaps the roles of normal and slip.
    strike_2, dip_2, rake_2 = vectors_to_plane(slip_1, normal_1)

    # Pressure (P), tension (T) and null (N) axes.
    p_axis = (normal_1 - slip_1) / np.sqrt(2)
    t_axis = (normal_1 + slip_1) / np.sqrt(2)
    n_axis = np.cross(normal_1, slip_1)

    new = pd.DataFrame({
        "strike_1": cat["strike"].to_numpy(float),
        "dip_1": cat["dip"].to_numpy(float),
        "rake_1": cat["rake"].to_numpy(float),
        "strike_2": strike_2, "dip_2": dip_2, "rake_2": rake_2,
    })
    for name, axis in (("p", p_axis), ("t", t_axis), ("n", n_axis)):
        trend, plunge = to_trend_plunge(axis)
        new[f"{name}_trend"], new[f"{name}_plunge"] = trend, plunge
    new = new.round(decimals)

    passthrough = [c for c in cat.columns if c not in ("event_id", "strike", "dip", "rake", "uncertainty_deg")]
    parts = [cat[["event_id"]], new]
    if "uncertainty_deg" in cat.columns:
        parts.append(cat[["uncertainty_deg"]])
    parts.append(cat[passthrough])
    return pd.concat(parts, axis=1)


def validate_nodal_planes(mechanisms: pd.DataFrame, source: str = "focal mechanisms") -> pd.DataFrame:
    """Check a table produced by :func:`compute_nodal_planes` (or by hand)."""
    required = ["event_id", "strike_1", "dip_1", "rake_1", "strike_2", "dip_2", "rake_2"]
    missing = [c for c in required if c not in mechanisms.columns]
    if missing:
        raise ValueError(
            f"{source}: missing column(s) {missing}. Run step 1 (compute_nodal_planes) first."
        )
    df = mechanisms.copy()
    df["event_id"] = df["event_id"].astype(str)
    for c in required[1:]:
        df[c] = pd.to_numeric(df[c], errors="raise")
    return df

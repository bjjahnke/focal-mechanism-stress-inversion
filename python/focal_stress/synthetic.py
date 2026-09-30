"""Generate synthetic focal mechanisms from a known stress state.

Used by the tests and the synthetic example to check that the pipeline gets
back the stress it was given.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .geometry import from_trend_plunge, vectors_to_plane
from .linear_inversion import StressResult, fault_instability


def stress_tensor(sigma1, sigma2, shape_ratio):
    """Build a stress tensor (NED, compression negative) from sigma_1 and sigma_2 [trend, plunge].

    sigma_3 is set perpendicular to both. Principal values are -1, -1 + 2R, 1.
    """
    v1 = from_trend_plunge(*sigma1)[0]
    v2 = from_trend_plunge(*sigma2)[0]
    v2 = v2 - (v2 @ v1) * v1
    v2 /= np.linalg.norm(v2)
    v3 = np.cross(v1, v2)
    values = np.array([-1.0, -1.0 + 2 * shape_ratio, 1.0])
    axes = np.vstack([v1, v2, v3])
    return axes.T @ np.diag(values) @ axes, axes, values


def make_synthetic_catalog(
    n_events: int = 40,
    sigma1=(0.0, 90.0),
    sigma2=(0.0, 0.0),
    shape_ratio: float = 0.5,
    friction: float = 0.6,
    rake_noise_deg: float = 0.0,
    min_instability: float = 0.7,
    seed: int = 0,
) -> pd.DataFrame:
    """Catalog of faults that slip in the direction of maximum shear under a known stress.

    Only faults with instability >= ``min_instability`` are kept (so they are
    plausible slip planes), and the catalog lists the true fault plane or the
    auxiliary plane at random, like a real catalog.
    """
    rng = np.random.default_rng(seed)
    tensor, axes, values = stress_tensor(sigma1, sigma2, shape_ratio)
    truth = StressResult(tensor, values, axes, shape_ratio, np.zeros(0), 0.0)

    rows = []
    while len(rows) < n_events:
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        if n[2] > 0:
            n = -n
        traction = tensor @ n
        shear = traction - (n @ traction) * n
        if np.linalg.norm(shear) < 1e-6:
            continue
        slip = shear / np.linalg.norm(shear)
        inst_fault = fault_instability(n[None], truth, friction)[0]
        inst_aux = fault_instability(slip[None], truth, friction)[0]
        if inst_fault < min_instability or inst_fault <= inst_aux:
            continue
        listed_is_fault = rng.random() < 0.5
        normal, s = (n, slip) if listed_is_fault else (slip, n)
        strike, dip, rake = (v[0] for v in vectors_to_plane(normal[None], s[None]))
        rake = (rake + rng.normal(0, rake_noise_deg) + 180) % 360 - 180 if rake_noise_deg else rake
        rows.append({
            "event_id": f"SYN{len(rows) + 1:03d}",
            "strike": round(float(strike), 1), "dip": round(float(dip), 1), "rake": round(float(rake), 1),
            "uncertainty_deg": max(rake_noise_deg, 5.0),
            "true_plane": 1 if listed_is_fault else 2,
        })
    return pd.DataFrame(rows)

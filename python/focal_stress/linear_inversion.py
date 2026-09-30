"""Linear stress inversion of fault-slip data (Michael, 1984).

Finds the deviatoric stress tensor whose shear traction on each fault plane
best lines up with that fault's slip direction, in the least-squares sense.
Assumes one uniform stress field and equal shear-traction magnitude on every
fault. Compression is negative, so sigma_1 (most compressive) is the smallest
eigenvalue.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import to_trend_plunge


@dataclass
class StressResult:
    tensor: np.ndarray            # (3, 3) deviatoric stress, NED, compression negative
    principal_values: np.ndarray  # (3,) sigma_1 <= sigma_2 <= sigma_3
    principal_axes: np.ndarray    # (3, 3) rows are unit vectors of sigma_1, sigma_2, sigma_3
    shape_ratio: float            # R = (s1 - s2) / (s1 - s3), 0..1
    misfit_deg: np.ndarray        # (n,) angle between observed slip and predicted shear traction
    variance: float               # sum of squared residuals / (3n - 5)

    def axes_trend_plunge(self) -> np.ndarray:
        """(3, 2) array of [trend, plunge] for sigma_1, sigma_2, sigma_3."""
        trend, plunge = to_trend_plunge(self.principal_axes)
        return np.column_stack([trend, plunge])


def _design_matrix(normals: np.ndarray) -> np.ndarray:
    """Stack the 3x5 blocks that map [s11, s12, s13, s22, s23] to shear traction."""
    n1, n2, n3 = normals[:, 0], normals[:, 1], normals[:, 2]
    rows = np.stack([
        np.stack([n1 - n1**3 + n1 * n3**2, n2 - 2 * n2 * n1**2, n3 - 2 * n3 * n1**2,
                  -n1 * (n2**2 - n3**2), -2 * n1 * n2 * n3], axis=1),
        np.stack([-n2 * (n1**2 - n3**2), n1 - 2 * n1 * n2**2, -2 * n1 * n2 * n3,
                  n2 - n2**3 + n2 * n3**2, n3 - 2 * n3 * n2**2], axis=1),
        np.stack([-n3 * n1**2 - n3 + n3**3, -2 * n1 * n2 * n3, n1 - 2 * n1 * n3**2,
                  -n2**2 * n3 - n3 + n3**3, n2 - 2 * n2 * n3**2], axis=1),
    ], axis=1)  # (n, 3, 5)
    return rows.reshape(-1, 5)


def invert_stress(normals, slips) -> StressResult:
    """Invert fault normals and slip vectors (NED unit vectors) for stress."""
    normals = np.asarray(normals, dtype=float)
    slips = np.asarray(slips, dtype=float)
    n = len(normals)
    if n < 2:
        raise ValueError("Stress inversion needs at least 2 faults (5 unknowns, 3 equations per fault)")

    G = _design_matrix(normals)
    d = slips.reshape(-1)
    m, *_ = np.linalg.lstsq(G, d, rcond=None)

    tensor = np.array([
        [m[0], m[1], m[2]],
        [m[1], m[3], m[4]],
        [m[2], m[4], -m[0] - m[3]],
    ])
    values, vectors = np.linalg.eigh(tensor)  # ascending: sigma_1 first
    s1, s2, s3 = values
    shape_ratio = float((s1 - s2) / (s1 - s3)) if s1 != s3 else np.nan

    predicted = (G @ m).reshape(n, 3)
    cos_fit = np.sum(predicted * slips, axis=1) / np.linalg.norm(predicted, axis=1)
    misfit = np.degrees(np.arccos(np.clip(cos_fit, -1.0, 1.0)))
    dof = max(3 * n - 5, 1)
    variance = float(np.sum((predicted.reshape(-1) - d) ** 2) / dof)

    return StressResult(tensor, values, vectors.T, shape_ratio, misfit, variance)


def fault_instability(normals, result: StressResult, friction: float) -> np.ndarray:
    """Fault instability I (Vavrycuk, 2014): 1 = optimally oriented for slip, 0 = most stable.

    Uses stress normalised so sigma_1 = 1, sigma_3 = -1 (compression positive)
    and the shape ratio from ``result``.
    """
    normals = np.asarray(normals, dtype=float)
    R = result.shape_ratio
    c1, c2, c3 = (normals @ result.principal_axes.T).T
    sigma_n = c1**2 + (1 - 2 * R) * c2**2 - c3**2
    tau = np.sqrt(np.maximum(c1**2 + (1 - 2 * R) ** 2 * c2**2 + c3**2 - sigma_n**2, 0.0))
    return (tau - friction * (sigma_n - 1)) / (friction + np.sqrt(1 + friction**2))

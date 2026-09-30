"""Vector geometry for fault planes and focal mechanisms.

Coordinate convention (Aki & Richards, 2002): x = North, y = East, z = Down.
Angles are in degrees:

* strike: 0-360, measured clockwise from North, with the dip to the right
  (right-hand rule)
* dip:    0-90, measured down from horizontal
* rake:   -180 to 180, direction the hanging wall slips, measured in the
  fault plane counter-clockwise from the strike direction
* trend/plunge: azimuth and downward inclination of a line (lower hemisphere)
"""

from __future__ import annotations

import numpy as np

_EPS = 1e-10


def plane_vectors(strike, dip, rake):
    """Return unit normal and slip vectors for fault planes.

    Parameters
    ----------
    strike, dip, rake : array-like, degrees, shape (n,)

    Returns
    -------
    normal, slip : ndarray, shape (n, 3), North-East-Down
        ``normal`` points up, out of the footwall; ``slip`` is the slip
        direction of the hanging wall.
    """
    phi, delta, lam = (np.radians(np.asarray(a, dtype=float)) for a in (strike, dip, rake))
    normal = np.column_stack([
        -np.sin(delta) * np.sin(phi),
        np.sin(delta) * np.cos(phi),
        -np.cos(delta),
    ])
    slip = np.column_stack([
        np.cos(lam) * np.cos(phi) + np.cos(delta) * np.sin(lam) * np.sin(phi),
        np.cos(lam) * np.sin(phi) - np.cos(delta) * np.sin(lam) * np.cos(phi),
        -np.sin(lam) * np.sin(delta),
    ])
    return normal, slip


def vectors_to_plane(normal, slip):
    """Return strike, dip, rake (degrees) from normal and slip vectors."""
    n = np.atleast_2d(np.asarray(normal, dtype=float)).copy()
    u = np.atleast_2d(np.asarray(slip, dtype=float)).copy()
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    u /= np.linalg.norm(u, axis=1, keepdims=True)

    # The normal must point up (z <= 0 in NED). Flipping it also flips the
    # slip vector so the pair still describes the same motion.
    down = n[:, 2] > 0
    n[down] *= -1
    u[down] *= -1

    dip = np.degrees(np.arccos(np.clip(-n[:, 2], -1.0, 1.0)))
    strike = np.mod(np.degrees(np.arctan2(-n[:, 0], n[:, 1])), 360.0)
    sin_dip = np.sin(np.radians(dip))
    horizontal = sin_dip < _EPS

    phi = np.radians(strike)
    along_strike = u[:, 0] * np.cos(phi) + u[:, 1] * np.sin(phi)
    with np.errstate(divide="ignore", invalid="ignore"):
        rake = np.degrees(np.arctan2(-u[:, 2] / sin_dip, along_strike))

    # Horizontal plane: strike is undefined, so fix it at 0 and measure rake
    # from North.
    if horizontal.any():
        strike[horizontal] = 0.0
        rake[horizontal] = np.degrees(np.arctan2(-u[horizontal, 1], u[horizontal, 0]))
    return strike, dip, rake


def to_trend_plunge(vectors):
    """Return trend and plunge (degrees, lower hemisphere) of NED vectors."""
    v = np.atleast_2d(np.asarray(vectors, dtype=float)).copy()
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    v[v[:, 2] < 0] *= -1
    plunge = np.degrees(np.arcsin(np.clip(v[:, 2], -1.0, 1.0)))
    trend = np.mod(np.degrees(np.arctan2(v[:, 1], v[:, 0])), 360.0)
    return trend, plunge


def from_trend_plunge(trend, plunge):
    """Return NED unit vectors from trend and plunge (degrees)."""
    t, p = np.radians(np.asarray(trend, dtype=float)), np.radians(np.asarray(plunge, dtype=float))
    return np.column_stack([np.cos(t) * np.cos(p), np.sin(t) * np.cos(p), np.sin(p)])

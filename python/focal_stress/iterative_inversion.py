"""Step 2: iterative stress inversion (Vavrycuk, 2014), repeated from random starts.

A focal mechanism cannot tell you which of its two nodal planes actually
slipped. Each run here:

1. picks plane 1 or plane 2 at random for every event,
2. inverts those planes for stress,
3. computes fault instability on both planes of every event under that stress,
4. switches each event to whichever plane is closer to failure,
5. repeats 2-4 until the set of chosen planes stops changing, starts
   repeating itself, or ``max_iterations`` is reached.

Because the answer can depend on the random start, many runs are made and
the results are tallied in step 3.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .catalog import MIN_EVENTS_FOR_INVERSION, RECOMMENDED_EVENTS
from .geometry import plane_vectors
from .linear_inversion import StressResult, fault_instability, invert_stress
from .nodal_planes import validate_nodal_planes

log = logging.getLogger(__name__)

RUN_COLUMNS = [
    "run_id", "plane_choice", "iterations", "status",
    "sigma1_trend", "sigma1_plunge", "sigma2_trend", "sigma2_plunge",
    "sigma3_trend", "sigma3_plunge", "shape_ratio",
    "misfit_mean_deg", "misfit_std_deg", "friction",
]


@dataclass
class PlaneSet:
    """Normals and slip vectors for both nodal planes of every event."""
    normals: np.ndarray  # (2, n, 3)
    slips: np.ndarray    # (2, n, 3)

    @classmethod
    def from_table(cls, mechanisms: pd.DataFrame) -> "PlaneSet":
        n1, s1 = plane_vectors(mechanisms["strike_1"], mechanisms["dip_1"], mechanisms["rake_1"])
        n2, s2 = plane_vectors(mechanisms["strike_2"], mechanisms["dip_2"], mechanisms["rake_2"])
        return cls(np.stack([n1, n2]), np.stack([s1, s2]))

    def select(self, choice: np.ndarray):
        """``choice`` holds 1 or 2 per event; returns (normals, slips) of the chosen planes."""
        idx = np.asarray(choice) - 1
        rows = np.arange(len(idx))
        return self.normals[idx, rows], self.slips[idx, rows]


def iterate_once(planes: PlaneSet, start: np.ndarray, friction: float, max_iterations: int):
    """Run one iterative inversion from a starting plane choice.

    Returns (choice, result, iterations, status). ``result`` is always the
    stress inverted from ``choice``. ``status`` is one of:

    * ``"converged"``: the chosen planes stopped changing.
    * ``"cycling"``: the choice flips back and forth between a few states,
      usually because some events have two planes almost equally close to
      failure. The state in the cycle with the lowest mean misfit is kept.
    * ``"max_iterations"``: stopped at the iteration cap.
    """
    choice = np.asarray(start, dtype=int).copy()
    history = []  # (choice, result) for every iteration, to detect cycles
    seen = {}
    for iteration in range(1, max_iterations + 1):
        result = invert_stress(*planes.select(choice))
        seen[encode_choice(choice)] = len(history)
        history.append((choice, result))

        inst_1 = fault_instability(planes.normals[0], result, friction)
        inst_2 = fault_instability(planes.normals[1], result, friction)
        new_choice = np.where(inst_1 >= inst_2, 1, 2)
        if np.array_equal(new_choice, choice):
            return choice, result, iteration, "converged"
        key = encode_choice(new_choice)
        if key in seen:
            cycle = history[seen[key]:]
            best_choice, best_result = min(cycle, key=lambda cr: cr[1].misfit_deg.mean())
            return best_choice, best_result, iteration, "cycling"
        choice = new_choice
    return history[-1][0], history[-1][1], max_iterations, "max_iterations"


def encode_choice(choice) -> str:
    """Plane choice as text, e.g. '1|2|2|1'. Text, not a number, so no tool mangles it."""
    return "|".join(str(int(c)) for c in choice)


def decode_choice(text: str) -> np.ndarray:
    return np.array([int(c) for c in str(text).split("|")], dtype=int)


def run_stress_inversions(
    mechanisms: pd.DataFrame,
    friction: float = 0.6,
    n_runs: int = 1000,
    max_iterations: int = 20,
    random_seed: int | None = None,
) -> pd.DataFrame:
    """Run ``n_runs`` iterative inversions from random starting plane choices.

    Parameters
    ----------
    mechanisms : output of step 1 (needs strike/dip/rake for planes 1 and 2).
    friction : fault friction coefficient used for instability (typically 0.4-0.8).
    n_runs : number of random starting choices.
    max_iterations : cap on iterations per run.
    random_seed : set for reproducible results.

    Returns
    -------
    DataFrame with one row per run and columns ``RUN_COLUMNS``.
    """
    if not 0 < friction <= 2:
        raise ValueError(f"friction should be between 0 and 2, got {friction}")
    if n_runs < 1 or max_iterations < 1:
        raise ValueError("n_runs and max_iterations must be at least 1")

    mech = validate_nodal_planes(mechanisms)
    n_events = len(mech)
    if n_events < MIN_EVENTS_FOR_INVERSION:
        raise ValueError(f"Stress inversion needs at least {MIN_EVENTS_FOR_INVERSION} events, got {n_events}")
    if n_events < RECOMMENDED_EVENTS:
        log.warning("Only %d events; results are more reliable with %d or more", n_events, RECOMMENDED_EVENTS)
    planes = PlaneSet.from_table(mech)
    rng = np.random.default_rng(random_seed)

    rows = []
    for run in range(1, n_runs + 1):
        start = rng.integers(1, 3, size=n_events)
        choice, result, iterations, status = iterate_once(planes, start, friction, max_iterations)
        rows.append(_run_row(run, choice, result, iterations, status, friction))

    runs = pd.DataFrame(rows, columns=RUN_COLUMNS)
    counts = runs["status"].value_counts()
    log.info("%d runs: %d converged, %d cycling, %d hit max_iterations; %d distinct plane choices",
             n_runs, counts.get("converged", 0), counts.get("cycling", 0),
             counts.get("max_iterations", 0), runs["plane_choice"].nunique())
    return runs


def _run_row(run, choice, result: StressResult, iterations, status, friction):
    axes = result.axes_trend_plunge()
    return [
        run, encode_choice(choice), iterations, status,
        *np.round(axes.ravel(), 2),
        round(result.shape_ratio, 4),
        round(float(result.misfit_deg.mean()), 2),
        round(float(result.misfit_deg.std(ddof=1)), 2),
        friction,
    ]

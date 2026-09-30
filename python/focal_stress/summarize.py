"""Step 3: tally the runs and check each stress solution against data uncertainty.

* Runs that end on the same set of planes give the same stress, so they are
  grouped into one **stress solution**. Solutions are ranked by how many runs
  reached them.
* For every event, count how often each nodal plane was chosen.
* If the catalog has ``uncertainty_deg``, each solution gets a variance test:
  misfit / uncertainty should have a variance of 1 if the misfits are no
  larger than the data uncertainty (the same test as MATLAB's ``vartest``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .iterative_inversion import PlaneSet, decode_choice
from .linear_inversion import fault_instability, invert_stress
from .nodal_planes import validate_nodal_planes


@dataclass
class Summary:
    stress_solutions: pd.DataFrame
    solution_planes: pd.DataFrame
    nodal_plane_frequency: pd.DataFrame


def variance_test(x, expected_variance: float = 1.0):
    """Two-sided chi-square test that ``x`` has variance ``expected_variance``.

    Returns (statistic, p_value). Matches MATLAB ``vartest(x, v)``.
    """
    x = np.asarray(x, dtype=float)
    dof = len(x) - 1
    statistic = dof * np.var(x, ddof=1) / expected_variance
    cdf = stats.chi2.cdf(statistic, dof)
    return float(statistic), float(min(1.0, 2 * min(cdf, 1 - cdf)))


def summarize_inversions(
    mechanisms: pd.DataFrame,
    runs: pd.DataFrame,
    friction: float | None = None,
    alpha: float = 0.05,
) -> Summary:
    """Build the three summary tables from step 1 and step 2 outputs."""
    mech = validate_nodal_planes(mechanisms)
    if friction is None:
        friction = float(runs["friction"].iloc[0])
    planes = PlaneSet.from_table(mech)
    has_uncertainty = "uncertainty_deg" in mech.columns
    n_runs = len(runs)

    grouped = (runs.groupby("plane_choice", sort=False)
                   .agg(n_runs=("run_id", "size"),
                        n_converged=("status", lambda x: int((x == "converged").sum())))
                   .reset_index()
                   .sort_values("n_runs", ascending=False, kind="stable")
                   .reset_index(drop=True))

    solution_rows, plane_rows = [], []
    for sid, row in enumerate(grouped.itertuples(index=False), start=1):
        choice = decode_choice(row.plane_choice)
        result = invert_stress(*planes.select(choice))
        axes = result.axes_trend_plunge()
        chosen_inst = fault_instability(planes.select(choice)[0], result, friction)

        sol = {
            "solution_id": sid,
            "n_runs": int(row.n_runs),
            "fraction_of_runs": round(row.n_runs / n_runs, 4),
            "n_runs_converged": int(row.n_converged),
            "sigma1_trend": axes[0, 0], "sigma1_plunge": axes[0, 1],
            "sigma2_trend": axes[1, 0], "sigma2_plunge": axes[1, 1],
            "sigma3_trend": axes[2, 0], "sigma3_plunge": axes[2, 1],
            "shape_ratio": result.shape_ratio,
            "misfit_mean_deg": result.misfit_deg.mean(),
            "misfit_std_deg": result.misfit_deg.std(ddof=1),
        }
        normalized = None
        if has_uncertainty:
            normalized = result.misfit_deg / mech["uncertainty_deg"].to_numpy(float)
            stat, p = variance_test(normalized)
            sol.update({
                "variance_test_statistic": stat,
                "variance_test_p_value": p,
                "fits_within_uncertainty": bool(p >= alpha),
            })
        sol["plane_choice"] = row.plane_choice
        solution_rows.append(sol)

        for i, event in enumerate(mech["event_id"]):
            k = choice[i]
            pr = {
                "solution_id": sid, "event_id": event, "chosen_plane": int(k),
                "strike": mech[f"strike_{k}"].iat[i], "dip": mech[f"dip_{k}"].iat[i],
                "rake": mech[f"rake_{k}"].iat[i],
                "misfit_deg": result.misfit_deg[i],
                "instability": chosen_inst[i],
            }
            if normalized is not None:
                pr["normalized_misfit"] = normalized[i]
            plane_rows.append(pr)

    stress_solutions = pd.DataFrame(solution_rows).round(4)
    solution_planes = pd.DataFrame(plane_rows).round(4)
    frequency = _plane_frequency(mech, runs)
    return Summary(stress_solutions, solution_planes, frequency)


def _plane_frequency(mech: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    choices = np.vstack([decode_choice(c) for c in runs["plane_choice"]])  # (runs, events)
    n_runs = len(runs)
    rows = []
    for i, event in enumerate(mech["event_id"]):
        counts = {1: int((choices[:, i] == 1).sum()), 2: int((choices[:, i] == 2).sum())}
        preferred = 1 if counts[1] >= counts[2] else 2
        for k in (1, 2):
            rows.append({
                "event_id": event, "plane": k,
                "strike": mech[f"strike_{k}"].iat[i], "dip": mech[f"dip_{k}"].iat[i],
                "rake": mech[f"rake_{k}"].iat[i],
                "times_chosen": counts[k],
                "fraction_chosen": round(counts[k] / n_runs, 4),
                "is_preferred": k == preferred,
            })
    return pd.DataFrame(rows)

import json

import numpy as np
import pandas as pd
import pytest

from focal_stress import (CatalogError, compute_nodal_planes, invert_stress, run_stress_inversions,
                          summarize_inversions, validate_catalog)
from focal_stress.geometry import plane_vectors, to_trend_plunge, vectors_to_plane
from focal_stress.iterative_inversion import decode_choice, encode_choice
from focal_stress.pipeline import Config, run_pipeline
from focal_stress.summarize import variance_test
from focal_stress.synthetic import make_synthetic_catalog


def angle_between_axes(a, b):
    """Angle in degrees between two lines given as [trend, plunge] (sign ignored)."""
    from focal_stress.geometry import from_trend_plunge
    va, vb = from_trend_plunge(*a)[0], from_trend_plunge(*b)[0]
    return np.degrees(np.arccos(np.clip(abs(va @ vb), -1, 1)))


def ang_diff(a, b):
    return np.abs((np.asarray(a) - np.asarray(b) + 180) % 360 - 180)


# ---------- geometry / step 1 ----------

def test_plane_round_trip():
    rng = np.random.default_rng(0)
    s, d, r = rng.uniform(0, 360, 200), rng.uniform(1, 89, 200), rng.uniform(-179, 179, 200)
    s2, d2, r2 = vectors_to_plane(*plane_vectors(s, d, r))
    assert ang_diff(s, s2).max() < 1e-8
    assert np.abs(d - d2).max() < 1e-8
    assert ang_diff(r, r2).max() < 1e-8


def test_auxiliary_plane_is_consistent():
    """Both nodal planes must give the same P and T axes, and plane 2 of plane 2 is plane 1."""
    rng = np.random.default_rng(1)
    cat = pd.DataFrame({"strike": rng.uniform(0, 360, 100), "dip": rng.uniform(1, 89, 100),
                        "rake": rng.uniform(-179, 179, 100)})
    m = compute_nodal_planes(cat, decimals=10)
    n2, u2 = plane_vectors(m.strike_2, m.dip_2, m.rake_2)
    p_trend, p_plunge = to_trend_plunge((n2 - u2) / np.sqrt(2))
    assert ang_diff(p_trend, m.p_trend).max() < 1e-6
    assert np.abs(p_plunge - m.p_plunge).max() < 1e-6
    back = compute_nodal_planes(m[["strike_2", "dip_2", "rake_2"]]
                                .rename(columns=lambda c: c[:-2]), decimals=10)
    assert ang_diff(back.strike_2, cat.strike).max() < 1e-6


def test_inversion_needs_three_events():
    mech = compute_nodal_planes(pd.DataFrame({"strike": [0, 90], "dip": [60, 60], "rake": [-90, -90]}))
    with pytest.raises(ValueError, match="at least 3 events"):
        run_stress_inversions(mech, n_runs=1)


def test_known_strike_slip_mechanism():
    # Vertical N-S right-lateral fault: auxiliary plane is vertical E-W, left-lateral.
    m = compute_nodal_planes(pd.DataFrame({"strike": [0], "dip": [90], "rake": [180]}), decimals=6)
    assert m.dip_2[0] == pytest.approx(90)
    assert m.strike_2[0] % 180 == pytest.approx(90)
    assert m.p_plunge[0] == pytest.approx(0, abs=1e-6)
    assert ang_diff(m.p_trend[0] % 180, 135) < 1e-6 or ang_diff(m.p_trend[0] % 180, 45) < 1e-6


# ---------- catalog validation ----------

def test_catalog_requires_columns():
    with pytest.raises(CatalogError, match="missing required column"):
        validate_catalog(pd.DataFrame({"strike": [1, 2, 3], "dip": [10, 20, 30]}))


def test_catalog_rejects_out_of_range():
    with pytest.raises(CatalogError, match="'dip' must be between"):
        validate_catalog(pd.DataFrame({"strike": [1, 2, 3], "dip": [10, 95, 30], "rake": [0, 0, 0]}))


def test_catalog_adds_event_ids_and_keeps_extra_columns():
    df = validate_catalog(pd.DataFrame({"Strike": [1, 2, 3], "dip": [10, 20, 30], "rake": [0, 0, 0],
                                        "magnitude": [1.1, 1.2, 1.3]}))
    assert list(df.event_id) == ["EV001", "EV002", "EV003"]
    assert "magnitude" in df.columns


# ---------- step 2 ----------

@pytest.fixture(scope="module")
def synthetic():
    return make_synthetic_catalog(n_events=40, sigma1=(30, 70), sigma2=(200, 20), shape_ratio=0.4, seed=3)


def test_linear_inversion_recovers_known_stress(synthetic):
    faults = synthetic.copy()
    # Use the true fault plane for every event.
    m = compute_nodal_planes(faults, decimals=10)
    k = faults.true_plane.to_numpy()
    s = np.where(k == 1, m.strike_1, m.strike_2)
    d = np.where(k == 1, m.dip_1, m.dip_2)
    r = np.where(k == 1, m.rake_1, m.rake_2)
    result = invert_stress(*plane_vectors(s, d, r))
    axes = result.axes_trend_plunge()
    # The linear method assumes equal shear stress on every fault, which the
    # synthetic faults do not have, so allow a few degrees of bias.
    assert angle_between_axes(axes[0], (30, 70)) < 8
    assert result.shape_ratio == pytest.approx(0.4, abs=0.2)
    assert np.median(result.misfit_deg) < 10


def test_matches_original_matlab_inversion():
    """Regression check against the original MATLAB inversion.m (Michael 1984 method).

    Reference values were produced by running inversion.m on the same planes.
    """
    from pathlib import Path
    catalog = pd.read_csv(Path(__file__).parents[2] / "examples/san_emidio_2016/catalog.csv")
    m = compute_nodal_planes(catalog, decimals=10)
    pick = np.random.default_rng(1).integers(1, 3, len(m))
    cols = [np.where(pick == 1, m[f"{c}_1"], m[f"{c}_2"]) for c in ("strike", "dip", "rake")]
    result = invert_stress(*plane_vectors(*cols))
    assert result.shape_ratio == pytest.approx(0.346933, abs=1e-6)
    expected_axes = [[18.098, 62.238], [160.959, 22.765], [257.428, 15.030]]
    assert np.allclose(result.axes_trend_plunge(), expected_axes, atol=1e-3)
    assert result.misfit_deg.mean() == pytest.approx(46.0660, abs=1e-4)


def test_iterative_inversion_finds_true_planes(synthetic):
    mech = compute_nodal_planes(synthetic)
    runs = run_stress_inversions(mech, friction=0.6, n_runs=50, random_seed=7)
    assert set(runs.status) <= {"converged", "cycling"}
    summary = summarize_inversions(mech, runs)
    best = summary.stress_solutions.iloc[0]
    assert angle_between_axes((best.sigma1_trend, best.sigma1_plunge), (30, 70)) < 5
    chosen = decode_choice(best.plane_choice)
    assert (chosen == synthetic.true_plane.to_numpy()).mean() > 0.9


def test_runs_are_reproducible_with_seed(synthetic):
    mech = compute_nodal_planes(synthetic)
    a = run_stress_inversions(mech, n_runs=20, random_seed=11)
    b = run_stress_inversions(mech, n_runs=20, random_seed=11)
    pd.testing.assert_frame_equal(a, b)


def test_choice_encoding_round_trip():
    choice = np.array([1, 2, 2, 1])
    assert encode_choice(choice) == "1|2|2|1"
    assert (decode_choice("1|2|2|1") == choice).all()


# ---------- step 3 ----------

def test_variance_test_matches_reference():
    # vartest([1 2 3 4 5], 1): chi-square statistic 10 on 4 degrees of freedom,
    # two-sided p = 2 * (1 - chi2cdf(10, 4)) = 0.080855
    stat, p = variance_test([1, 2, 3, 4, 5])
    assert stat == pytest.approx(10.0)
    assert p == pytest.approx(0.080855, abs=1e-5)


def test_frequency_table_sums_to_one(synthetic):
    mech = compute_nodal_planes(synthetic)
    runs = run_stress_inversions(mech, n_runs=30, random_seed=1)
    freq = summarize_inversions(mech, runs).nodal_plane_frequency
    totals = freq.groupby("event_id").fraction_chosen.sum()
    assert np.allclose(totals, 1.0)


# ---------- whole pipeline ----------

def test_pipeline_writes_all_outputs(tmp_path, synthetic):
    synthetic.drop(columns="true_plane").to_csv(tmp_path / "catalog.csv", index=False)
    (tmp_path / "config.json").write_text(json.dumps({
        "input_catalog": "catalog.csv", "output_dir": "out", "n_runs": 25,
        "random_seed": 5, "make_plots": True,
    }))
    written = run_pipeline(Config.from_file(tmp_path / "config.json"))
    for key in ("focal_mechanisms", "inversion_runs", "stress_solutions", "solution_planes",
                "nodal_plane_frequency", "metadata"):
        assert written[key].exists()
    assert (tmp_path / "out" / "plots" / "stress_solutions.png").exists()
    sols = pd.read_csv(written["stress_solutions"])
    assert "fits_within_uncertainty" in sols.columns
    assert sols.n_runs.sum() == 25

"""Create a synthetic catalog from a known stress state.

Run from the repository root:
    python examples/synthetic/make_catalog.py

Writes catalog.csv (the pipeline input) and true_answer.json (what the
pipeline should recover). The catalog lists the real fault plane for about
half the events and the auxiliary plane for the rest, like a real catalog.
"""

import json
from pathlib import Path

from focal_stress.synthetic import make_synthetic_catalog

HERE = Path(__file__).parent
TRUTH = {
    "sigma1_trend_plunge": [300, 80],
    "sigma2_trend_plunge": [30, 0],
    "shape_ratio": 0.5,
    "friction": 0.6,
    "rake_noise_deg": 10,
    "n_events": 40,
    "seed": 2024,
}

catalog = make_synthetic_catalog(
    n_events=TRUTH["n_events"],
    sigma1=TRUTH["sigma1_trend_plunge"],
    sigma2=TRUTH["sigma2_trend_plunge"],
    shape_ratio=TRUTH["shape_ratio"],
    friction=TRUTH["friction"],
    rake_noise_deg=TRUTH["rake_noise_deg"],
    seed=TRUTH["seed"],
)
TRUTH["true_plane_by_event"] = dict(zip(catalog.event_id, catalog.true_plane.astype(int).tolist()))
catalog.drop(columns="true_plane").to_csv(HERE / "catalog.csv", index=False)
(HERE / "true_answer.json").write_text(json.dumps(TRUTH, indent=2))
print(f"Wrote {len(catalog)} events to {HERE / 'catalog.csv'}")

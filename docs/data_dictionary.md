# Data dictionary

Every file the pipeline reads or writes, column by column. The Python and MATLAB versions read and write exactly the same files.

**True/false columns** are written as `True`/`False` by Python and `1`/`0` by MATLAB.

**Angle conventions** (Aki & Richards, 2002). All angles are in degrees.

| Term | Range | Meaning |
|---|---|---|
| strike | 0–360 | Azimuth of the plane's horizontal line, clockwise from North, with the plane dipping to the right (right-hand rule) |
| dip | 0–90 | Angle down from horizontal |
| rake | −180 to 180 | Slip direction of the hanging wall, measured in the plane counter-clockwise from strike. Positive = reverse component, negative = normal component |
| trend | 0–360 | Azimuth of a line, clockwise from North |
| plunge | 0–90 | Angle of a line below horizontal (lower hemisphere) |

---

## Input: catalog CSV

One row per earthquake. Comma-separated, with a header row. Column names are not case-sensitive.

| Column | Required? | Type | Rules | Used for |
|---|---|---|---|---|
| `strike` | **Required** | number | 0–360 | Step 1 |
| `dip` | **Required** | number | 0–90 | Step 1 |
| `rake` | **Required** | number | −180 to 180 | Step 1 |
| `event_id` | Optional | text | Must be unique. Created as `EV001`, `EV002`, … if missing | Joins every output back to the input |
| `uncertainty_deg` | Optional | number | > 0 | Step 3 variance test. Without it, the test is skipped |
| anything else | Optional | any | Not checked | Copied unchanged into `focal_mechanisms.csv` |

`strike`, `dip` and `rake` can describe **either** nodal plane; it doesn't matter which one your catalog lists.

**Minimum size.** Step 1 works on any number of events. Steps 2–3 need at least 3 and warn below 10; the results are more meaningful with 20 or more.

**Validation.** Before anything runs, the catalog is checked. A bad file stops with a message that names the column and the row numbers, for example:

```
catalog.csv: 'dip' must be between 0 and 90; rows 14, 22 are not
```

Example (`examples/san_emidio_2016/catalog.csv`):

```csv
event_id,origin_time,latitude,longitude,easting_m,northing_m,depth_km,strike,dip,rake,quality,uncertainty_deg
SE2016-01,2016-12-07T18:42:01.62,40.38267,-119.40282,296035.797,4473003.038,-0.527,260,17,155,C,45
```

---

## Config file (JSON)

Used by `focal-stress run` (Python) and `run_pipeline` (MATLAB). Relative paths are resolved from the config file's own folder.

| Key | Required? | Default | Meaning |
|---|---|---|---|
| `input_catalog` | **Required** | — | Path to the catalog CSV |
| `output_dir` | **Required** | — | Folder for all outputs (created if needed) |
| `friction` | Optional | `0.6` | Friction coefficient used to judge which plane is closer to failure. Typical rock values are 0.4–0.8 |
| `n_runs` | Optional | `1000` | Number of random starting plane choices |
| `max_iterations` | Optional | `20` | Cap on iterations per run (most runs finish in 2–5) |
| `random_seed` | Optional | none | Set it to get the same results every time |
| `significance_level` | Optional | `0.05` | Alpha for the variance test |
| `make_plots` | Optional | `true` | Write the PNG stereonets |

---

## Step 1 output: `focal_mechanisms.csv`

One row per event.

| Column | Meaning |
|---|---|
| `event_id` | From the input, or generated |
| `strike_1`, `dip_1`, `rake_1` | The plane given in the input |
| `strike_2`, `dip_2`, `rake_2` | The second (auxiliary) nodal plane, computed |
| `p_trend`, `p_plunge` | P (pressure) axis |
| `t_trend`, `t_plunge` | T (tension) axis |
| `n_trend`, `n_plunge` | N (null) axis, the line where the two planes meet |
| `uncertainty_deg` | Copied from the input, if present |
| *other input columns* | Copied unchanged |

## Step 2 output: `inversion_runs.csv`

One row per run. Each run starts from a random choice of plane for every event and iterates until the choice settles.

| Column | Meaning |
|---|---|
| `run_id` | 1, 2, 3, … |
| `plane_choice` | Final plane chosen for each event, in catalog order, e.g. `2\|1\|1\|2`. Stored as text so no tool reads it as a number |
| `iterations` | Iterations used |
| `status` | `converged`: the choice stopped changing. `cycling`: it flipped between a few states, usually because an event's two planes are almost equally close to failure; the best-fitting state in the cycle is kept. `max_iterations`: stopped at the cap |
| `sigma1_trend` … `sigma3_plunge` | Principal stress directions. σ1 is the most compressive, σ3 the least |
| `shape_ratio` | R = (σ1 − σ2) / (σ1 − σ3), from 0 to 1. Describes the relative size of σ2 |
| `misfit_mean_deg`, `misfit_std_deg` | Angle between each event's observed slip and the slip the stress predicts |
| `friction` | Friction used (kept here so step 3 can run from this file alone) |

## Step 3 outputs

### `stress_solutions.csv`

Runs that end on the same set of planes give the same stress, so they are grouped into one solution. One row per solution, most common first.

| Column | Meaning |
|---|---|
| `solution_id` | 1 = reached by the most runs |
| `n_runs`, `fraction_of_runs` | How many runs reached this solution. A high fraction means the answer does not depend on the starting guess |
| `n_runs_converged` | How many of those runs had `status = converged` |
| `sigma1_trend` … `sigma3_plunge`, `shape_ratio` | The stress state |
| `misfit_mean_deg`, `misfit_std_deg` | Fit to the data. Lower is better |
| `variance_test_statistic`, `variance_test_p_value` | *Only if `uncertainty_deg` was given.* Chi-square test that misfit ÷ uncertainty has a variance of 1 (same as MATLAB `vartest`) |
| `fits_within_uncertainty` | *Only if `uncertainty_deg` was given.* `True` when p ≥ `significance_level`, meaning one uniform stress state explains the data within their uncertainty |
| `plane_choice` | The set of planes behind this solution (links to `inversion_runs.csv`) |

### `solution_planes.csv`

The plane each solution picked as the one that slipped, for every event. Rows = solutions × events.

| Column | Meaning |
|---|---|
| `solution_id`, `event_id` | Keys |
| `chosen_plane` | 1 or 2 |
| `strike`, `dip`, `rake` | The chosen plane |
| `misfit_deg` | This event's misfit under this solution |
| `instability` | How close the plane is to failure: 1 = optimally oriented, 0 = most stable |
| `normalized_misfit` | *Only if `uncertainty_deg` was given.* `misfit_deg / uncertainty_deg` |

### `nodal_plane_frequency.csv`

How often each nodal plane was chosen across all runs. Two rows per event.

| Column | Meaning |
|---|---|
| `event_id`, `plane` | Keys (`plane` is 1 or 2) |
| `strike`, `dip`, `rake` | The plane |
| `times_chosen`, `fraction_chosen` | Across all runs. The two fractions for an event add to 1 |
| `is_preferred` | `True` for the plane chosen more often |

### `run_metadata.json`

Settings, software versions, run time, counts of run status and the list of files written, so every result can be traced back to how it was made.

### `plots/`

| File | Shows |
|---|---|
| `stress_solutions.png` | Stereonet for each of the top four solutions: σ1, σ2, σ3 and the poles of the chosen planes |
| `nodal_plane_frequency.png` | Poles of both planes of every event, shaded by how often each was chosen |

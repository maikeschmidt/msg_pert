# msg_pert — MSG Perturbation Analysis Toolbox

**Systematic perturbation analysis for Magnetospinography (MSG) forward models:
generating anatomically shifted geometry files and perturbed tissue
conductivities, then quantifying the sensitivity of BEM, FEM, Biot-Savart, and
single-sphere forward models to perturbations in source space, sensor array
position, and tissue conductivity — for MSG and ESG side by side, noise-free
and through realistic sensor noise.**

Developed by **Maike Schmidt** at the **Department of Imaging Neuroscience,
University College London**.

---

## Overview

This toolbox answers questions about MSG forward models:

> *Which forward model is most sensitive to imperfect anatomical knowledge?*
> *How large a shift in anatomy still produces consistent source reconstruction?*
> *And once realistic sensor noise is added, how much of the field survives?*

Starting from an original geometry file produced by `msg_coreg`, msg_pert
generates three sets of perturbed configurations and provides the tools to
analyse the resulting forward fields once they have been computed in `msg_fwd`.

Two **modalities** are supported side by side: **MSG** (triaxial magnetometer
array) and **ESG** (tangential/radial surface electrodes). Both are configured
in one place and the analysis pipeline loops over them, so nothing has to be
edited between modalities. A final combined stage compares the two directly.

### Perturbation types

| Type | What is perturbed | What stays fixed | Configurations |
|---|---|---|---|
| **Source space** | Spinal cord mesh, bone mesh, source positions | Torso, heart, lungs | 24: 3 bundles × 8 random shifts |
| **Sensor array** | Sensor coil / electrode positions | All anatomy meshes | 24: 3 bundles × 8 random shifts |
| **Conductivity** (BEM only) | Tissue compartment conductivities | Geometry (original) | 24: 3 bundles × 8 random scalings |

All three use the same bundle structure: three error bundles (small / medium /
large) with 8 random realisations each.

**Source space shifts** model uncertainty in anatomical understanding — the whole
cord–bone unit moves together. Bundles are ~2 mm (small), ~5 mm (medium),
~10 mm (large) per axis.

**Sensor array shifts** model registration error — e.g. the scanner cast placed a
few mm off its true position. Same ~2/5/10 mm bundle scales. For ESG (surface
electrodes) the shift is applied in X and Y only (Z held at 0).

**Conductivity perturbations** scale each BEM tissue compartment's conductivity
upward by a random factor (bundles: up to +5%, +10%, +50%). Because HBF bakes
conductivity into the BEM transfer matrices, each perturbation rebuilds the head
model (see `msg_fwd/run_conductivity_perturbation.m`).

---

## Workflow

```
Phase 1 — perturbation generation (once per study)

  msg_coreg geometry file
         │
         ├─ pt_generate_source_shifts   →  24 source-shift geometries
         ├─ pt_generate_sensor_shifts   →  24 sensor-shift geometries
         ▼
  msg_fwd: BEM (MSG + ESG), Biot-Savart (MSG), conductivity perturbation,
           optional FEM — on every geometry, for the analysed array
         │
         ▼
Phase 2 — staged analysis  (run_perturbation_analysis)

  0  pt_load_leadfields, pt_compute_metrics   every metric, every perturbation
  1  pt_baseline            noise-free fields: topoplots, amplitude, model type
  2  pt_within_modality     MSG: perturbation types and sizes, with statistics
  3  pt_within_modality     ESG: the same
  4  pt_compare_modalities  MSG vs ESG under identical perturbations
  5  pt_noise_simulate      stages 2–4 repeated through realistic sensor noise
     pt_noise_analyse       (SQUID MSG, OP-MSG, ESG)
  6  pt_summary_table       every source of error on one scale
```

### What each stage answers

| Stage | Question |
|---|---|
| 1 | What do Biot-Savart MSG, BEM MSG and BEM ESG predict before anything goes wrong, and how different are two reasonable model types? |
| 2, 3 | Within one modality: how large is each perturbation, where along the cord, is it an amplitude or a topography change, does it scale with shift size, which type matters most, and does the forward model type change the sensitivity? |
| 4 | Do MSG and ESG respond differently to the same perturbation, and do they rank the perturbation types differently? |
| 5 | Once realistic sensor noise is added, which perturbations are still visible, at what noise level do they disappear, and do the stage-4 differences survive? |
| 6 | Of model type, geometry, conductivity and noise, which dominates? |

Every stage reports the same four metrics (RE, r², RDM and gain from lnMAG)
and uses the same statistics, so numbers can be compared across stages.
**[RUN_ORDER.md](RUN_ORDER.md)** says what to run and what it costs;
**[INTERPRETATION.md](INTERPRETATION.md)** says how to read every output.

---

## Companion Repositories

This toolbox is part of the MSG toolbox family and depends on:

**msg_coreg** — MSG Coregistration Toolbox (produces the input geometry files)  
https://github.com/maikeschmidt/msg_coreg

**msg_fwd** — MSG Forward Modelling Toolbox (runs BEM/FEM/Biot-Savart/sphere on shifted geometries)  
https://github.com/maikeschmidt/msg_fwd

Both must be cloned as sibling directories to `msg_pert`.

---

## Directory Structure

```
msg_pert/
├── config_pert.m                  — every path and setting
├── run_perturbation_analysis.m    — PHASE 2 master: runs stages 0–6
├── RUN_ORDER.md                   — what to run, in what order, at what cost
├── INTERPRETATION.md              — how to read every output
│
├── pt_generate_source_shifts.m    — PHASE 1: 24 source-shift geometries
├── pt_generate_sensor_shifts.m    — PHASE 1: 24 sensor-shift geometries
│
├── pt_load_leadfields.m           — stage 0: load one array, scale each file, check
├── pt_compute_metrics.m           — stage 0: RE, r², RDM, lnMAG for every perturbation
├── pt_baseline.m                  — stage 1: noise-free baseline
├── pt_within_modality.m           — stages 2/3: within MSG / within ESG
├── pt_compare_modalities.m        — stage 4: MSG vs ESG
├── pt_noise_simulate.m            — stage 5a: perturbation + sensor noise
├── pt_noise_analyse.m             — stage 5b: figures, tables, statistics
├── pt_summary_table.m             — stage 6: every factor on one scale
│
├── pt_modality.m                  — holds the active modality between scripts
├── pt_run_step.m                  — runs a stage in an isolated workspace
├── pt_add_functions.m, pert_path.m — path setup
├── pt_diagnose_leadfields.m       — inspect a loaded set when results look wrong
│
├── pt_compute_rsq.m, pt_plot_curves.m, pt_plot_heatmaps.m,
│   pt_plot_displacement.m, pt_plot_slope_vs_position.m,
│   pt_compute_table.m, pt_compare_perturbations.m
│                                  — earlier r²-only outputs (run_legacy = true)
│
├── functions/                     — metric, statistics, table and plot helpers
│   ├── pt_metrics_block.m         — vectorised lf_metrics (identical numbers)
│   ├── pt_lf_matrix.m             — comparison vectors per orientation / axis set
│   ├── pt_unit_scale.m            — per-file unit detection
│   ├── pt_compare_two.m           — paired / unpaired test + effect size
│   ├── pt_stats_within.m, pt_stats_between.m, pt_fdr_by_group.m
│   ├── pt_perm_test2.m, pt_cliffs_delta.m, pt_spearman_perm.m
│   ├── pt_describe.m, pt_quantile.m, pt_shift_summary.m, pt_metric_info.m
│   └── pt_write_table.m, pt_long_row.m, pt_plot_band.m, pt_save_figure.m, ...
│
├── tests/
│   ├── test_pt_core.m             — metric and statistics core vs known answers
│   └── test_pt_pipeline.m         — every stage on synthetic lead fields
│
└── simulations/                   — standalone noise-simulation package
```

---

## Requirements

1. **MATLAB** (R2020a or later recommended)

2. **SPM** (developmental version)  
   https://www.fil.ion.ucl.ac.uk/spm/

3. **FieldTrip** (bundled with SPM — do not install standalone)

4. **Helsinki BEM Framework (HBF)** by Matti Stenroos  
   Clone into `msg_coreg/hbf_lc_p`:  
   https://github.com/MattiStenroos/hbf_lc_p

5. **msg_coreg** sibling repository  
   https://github.com/maikeschmidt/msg_coreg

6. **msg_fwd** sibling repository  
   https://github.com/maikeschmidt/msg_fwd

---

## Getting Started

### Step 1: Set up dependencies

```matlab
pt_add_functions;   % adds msg_pert to path; checks msg_coreg, msg_fwd, HBF
```

### Step 2: Configure

Everything is set in `config_pert.m`. Paths, sensor description and method
availability are declared **per modality**, in one block each:

```matlab
base_geom_name = 'original';   % short stem in file names, WITHOUT the
                               % leading 'geometries_' prefix

mods_cfg.msg.geoms_path           = '';   % original geometry .mat
mods_cfg.msg.perturbed_geoms_path = '';   % output for shifted geometry files
mods_cfg.msg.forward_fields_base  = '';   % msg_fwd leadfield output
mods_cfg.msg.save_base_dir        = '';   % figures and tables
mods_cfg.msg.sensor_n_axes        = 3;    % 3 = triaxial MSG
mods_cfg.msg.sensor_is_meg        = true;
mods_cfg.msg.have_bem             = true; % which forward models you computed
% ... and the same block again for mods_cfg.esg (sensor_n_axes = 2)

combined_results_dir = '';                % MSG-vs-ESG comparison output
pert_modalities      = {'msg', 'esg'};    % e.g. {'msg'} to run MSG only
```

`sensor_n_axes` and `sensor_is_meg` are **declared, not inferred**: an ESG
electrode count can also be divisible by 3, so guessing would mis-split the
leadfield.

### Step 3: Generate source-shift geometries (Phase 1a)

```matlab
pt_generate_source_shifts;
```

This produces 24 shifted geometry files (3 bundles × 8 random shifts) in
`perturbed_geoms_path` and prints the shift vectors to paste into `config_pert.m`
plus a filename list for `msg_fwd`.

### Step 4: Generate sensor-shift geometries (Phase 1b)

```matlab
pt_generate_sensor_shifts;
```

This produces 24 shifted geometry files and prints:
1. The shift vectors — paste into `config_pert.m` under `sensor_shift_vectors`
2. A filename list for `msg_fwd`

### Step 5: Run forward models in msg_fwd

Copy the printed filename lists into `msg_fwd` and run the forward models you need
for both front and back arrays:

- `run_bem_leadfields.m` — BEM (required)
- `run_fem_leadfields.m`, Biot-Savart, single sphere — optional extra methods
- `run_conductivity_perturbation.m` — BEM conductivity leadfields (for the
  conductivity mode)

The per-geometry leadfield `.mat` files are the input to Phase 2. (msg_pert
organises them itself via `pt_load_leadfields` in the next step — you do **not**
run msg_fwd's own `load_and_organise_leadfields` here.)

### Step 6: Run the analysis (Phase 2)

Set the `have_<method>` flags per modality in `config_pert.m` to match the
forward models you computed, set `stage_results_dir`, then:

```matlab
run_perturbation_analysis;
```

Flags at the top of that script run a subset of stages. Any stage can be
re-run alone once its inputs exist, see [RUN_ORDER.md](RUN_ORDER.md).
Outputs go to numbered folders under `stage_results_dir`:

| Folder | Contents |
|---|---|
| `1_baseline` | topoplots, amplitude and focality along the cord, model-type decomposition |
| `2_within_msg`, `3_within_esg` | decomposition per bundle, dose response, types and methods compared, statistics tables |
| `4_msg_vs_esg` | headline comparison, effect-size maps, interaction, slopes, along-cord tests |
| `5_noise` | noise alone, perturbation + noise, detectability, critical noise levels, error budget, systems compared per level |
| `6_summary` | one ranking of every source of error per modality |

Every table is written as `.csv` (every orientation and axis set), `.txt`
(readable) and `.tex` (booktabs; headline orientation, whole array).

---

## Script Reference

| Script | Stage | Description |
|---|---|---|
| `pt_generate_source_shifts` | Phase 1 | 24 geometries: 3 bundles × 8 random source-space shifts |
| `pt_generate_sensor_shifts` | Phase 1 | 24 geometries: 3 bundles × 8 random sensor-array shifts |
| `run_perturbation_analysis` | all | Master script: stages 0–6, both modalities |
| `pt_load_leadfields` | 0 | Load the analysed array for every method and perturbation; per-file unit scaling; magnitude check against each reference |
| `pt_compute_metrics` | 0 | RE, r², RDM and lnMAG per source, method, orientation and axis set → `pert_metrics.mat` |
| `pt_baseline` | 1 | Unperturbed topoplots, amplitude and focality along the cord, forward-model-type decomposition |
| `pt_within_modality` | 2, 3 | Within one modality: descriptive table, decomposition per bundle, dose response, type and method comparisons, along-cord source-vs-sensor test |
| `pt_compare_modalities` | 4 | MSG vs ESG: descriptive, tests, interaction, dose-response slopes, along-cord tests, effect-size maps |
| `pt_noise_simulate` | 5 | Every perturbed field measured through trial-averaged sensor noise, scored against the true field |
| `pt_noise_analyse` | 5 | Noise alone, perturbation + noise, detectability, critical noise level, error budget, system comparisons per level, within-modality tests under noise |
| `pt_summary_table` | 6 | Every source of error ranked on one scale |
| `pt_modality` | — | Get/set the active modality; survives each script's `clearvars` |
| `pt_run_step` | — | Run one script in an isolated workspace |
| `pt_diagnose_leadfields` | — | Inspect a loaded lead-field set |

---

## Standalone noise simulation (`simulations/`)

The `simulations/` package predates stage 5 and remains usable on its own. It
scores the noisy sensor-by-time data matrix with r², for a single
representative shift per bundle. Stage 5 of the main pipeline asks the same
question for every realisation with all four metrics and statistics, so
use stage 5 for results that sit alongside stages 1–4. The noise floors and
source waveform in `config_pert.m` are copied from `simulations/config_sim.m`;
keep the two in step if you change either. See `simulations/README.md`.

---

## Metrics

All metrics come from msg_fwd's shared definitions (`lf_metrics`), computed
here by the vectorised `pt_metrics_block`, which gives identical numbers. The
unperturbed geometry is always the reference.

| Metric | Measures |
|---|---|
| RE (%) | `‖L_pert − L_orig‖₂ / ‖L_orig‖₂ × 100`, magnitude and shape |
| r² | squared Pearson correlation, shape only |
| RDM | distance between unit-normalised fields, shape only |
| gain (%) | `(exp(lnMAG) − 1) × 100`, magnitude only |

`RE ≈ √(gain² + (RDM×100)²)`, so every result can be split into an
amplitude part and a topography part. Statistics treat one perturbation
realisation as one observation (its median over the cord), use permutation
tests with rank-biserial or Cliff's δ effect sizes, and are FDR-corrected
within each family and metric. See [INTERPRETATION.md](INTERPRETATION.md).

---

## Coordinate Convention

All geometries follow the **msg_coreg scanner-cast frame** (mm):

| Axis | Direction |
|---|---|
| X | Left → Right |
| Y | Posterior → Anterior (Rostral → Caudal along cord) |
| Z | Inferior → Superior (Ventral → Dorsal for cord cross-section) |

Source shifts are applied in this frame: a +Y shift moves the cord anterior
relative to the sensor array; a +Z shift moves the cord superior.

---

## Citation

If you use this toolbox, please cite it along with the companion toolboxes you
used:

> msg_coreg: https://github.com/maikeschmidt/msg_coreg  
> msg_fwd:   https://github.com/maikeschmidt/msg_fwd

---

## Contact

For questions, issues, or contributions, open an issue or pull request on GitHub.  
Contact: maike.schmidt.23@ucl.ac.uk

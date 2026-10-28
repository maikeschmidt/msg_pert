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
generates three sets of perturbed configurations, computes their forward
fields with `msg_fwd`'s solvers, and analyses them.

### Forward models

| Model | Modalities | Role |
|---|---|---|
| **BEM** (MRI-derived) | MSG, ESG | the reference: every model comparison is measured against it |
| **FEM** | MSG | the second volume-conductor model; DUNEuro as msg_fwd calls it is MEG only |
| **Biot-Savart** | MSG | anatomy-free: the purely geometric sensitivity, and the upper bound through noise |

Every perturbation is solved with every model that applies, so one shift can
be read through all three.

Two **modalities** are supported side by side: **MSG** (triaxial magnetometer
array) and **ESG** (one set of surface electrodes, common-average referenced).
Both are configured
in one place and the analysis pipeline loops over them, so nothing has to be
edited between modalities. A final combined stage compares the two directly.

### Perturbation types

| Type | What is perturbed | What stays fixed | Configurations |
|---|---|---|---|
| **Source space** | Spinal cord mesh, bone mesh, source positions | Torso, heart, lungs | 24: 3 bundles × 8 random shifts |
| **Sensor array** | Sensor coil / electrode positions | All anatomy meshes | 24: 3 bundles × 8 random shifts |
| **Conductivity** (BEM, FEM) | Tissue compartment conductivities | Geometry (original) | 24: 3 bundles × 8 random scalings |

All three use the same bundle structure: three error bundles (small / medium /
large) with 8 random realisations each.

**Source space shifts** model uncertainty in anatomical understanding — the whole
cord–bone unit moves together. Bundles are ~2 mm (small), ~5 mm (medium),
~10 mm (large) per axis.

**Sensor array shifts** model registration error — e.g. the scanner cast placed a
few mm off its true position. Same ~2/5/10 mm bundle scales. For ESG (surface
electrodes) the shift is applied in X and Y only (Z held at 0).

**Conductivity perturbations** scale each tissue compartment's conductivity
upward by a random factor (bundles: up to +5%, +10%, +50%). BEM and FEM get the
same draws, so the two solvers are compared on identical changes. Because HBF
bakes conductivity into the BEM transfer matrices, each BEM draw rebuilds the
head model; the FEM reuses its mesh and changes the labels' conductivities.

**FEM re-meshing.** TetGen places its region seeds at random, so two FEM meshes
of one geometry give slightly different lead fields. A sensor shift or
conductivity draw reuses the unshifted mesh and contains no re-meshing; a
source shift moves the cord and bone and needs a new mesh. The unshifted
geometry is therefore also meshed with 8 other seeds, and that spread — the
re-meshing floor — is reported beside every FEM source-shift result.

---

## Workflow

Two scripts run everything, with every folder set in `config_pert_paths.m`:

```
Phase 1 — geometries and lead fields  (run_all_perturbations, once per study)

  msg_coreg geometry file  (<data_root>/geometries/<base_geom_name>.mat)
         │
         ├─ pt_generate_source_shifts   →  24 source-shift geometries
         ├─ pt_generate_sensor_shifts   →  24 sensor-shift geometries
         ▼
  pt_compute_leadfields, with msg_fwd's solvers:
         BEM (MSG + ESG), FEM (MSG), Biot-Savart (MSG),
         conductivity perturbation (BEM + FEM), FEM re-meshing control
         — on every geometry, for the analysed array
         │
         ▼
Phase 2 — staged analysis  (run_all_pert_analysis)

  0  pt_load_leadfields, pt_compute_metrics   every metric, every perturbation
  1  pt_baseline            noise-free fields: topoplots, amplitude, model type
  2  pt_within_modality     MSG: perturbation types and sizes, with intervals
  3  pt_within_modality     ESG: the same
  4  pt_compare_modalities  MSG vs ESG under identical perturbations
  5  pt_noise_simulate      an N13-like evoked response at C7 through realistic
     pt_noise_analyse       sensor noise (SQUID MSG, OP-MSG, ESG through BEM;
                            the MSG systems also through FEM and through the
                            Biot-Savart upper bound), stages 2–4 repeated
  6  pt_summary_table       every source of error on one scale
  7  pt_cord_activity       N11 volley + N13 as a recording, N13 extent swept:
                            BEM vs FEM divergence, detectability, perturbations
```

### What each stage answers

| Stage | Question |
|---|---|
| 1 | What do Biot-Savart, BEM and FEM MSG and BEM ESG predict before anything goes wrong, and how different are the model types (MRI-derived BEM as the reference)? |
| 2, 3 | Within one modality: how large is each perturbation, where along the cord, is it an amplitude or a topography change, does it scale with shift size, which type matters most, do BEM and FEM respond alike, how much of the sensitivity does the anatomy add over the anatomy-free Biot-Savart model, and is a FEM source-shift effect larger than re-meshing alone? |
| 4 | Do MSG and ESG respond differently to the same perturbation, and do they rank the perturbation types differently? |
| 5 | For an N13-like evoked response at C7: how well is it recovered through each system's noise with perfect coregistration, and once the geometry is wrong, which perturbations are still visible, at what noise level do they disappear, do the stage-4 differences survive, and how far below the Biot-Savart upper bound does each MSG system sit? |
| 6 | Of model type, geometry, conductivity and noise, which dominates — over the cord and at C7? |
| 7 | With the N11 volley and the N13 together as a recording, and the N13 spread from one source to a wide patch: do BEM and FEM diverge more for spread sources, is a spread N13 more or less detectable through each system's noise under each model, and does spreading change how robust each model is to the perturbations? |

Every stage reports the same four metrics (RE, r², RDM and gain from lnMAG)
and the same kind of summary — medians with bootstrap intervals, as msg_fwd
and msg_models do — so numbers can be compared across stages and across the
toolboxes.
**[RUN_ORDER.md](RUN_ORDER.md)** says what to run and what it costs;
**[INTERPRETATION.md](INTERPRETATION.md)** says how to read every output.

---

## Companion Repositories

This toolbox is part of the MSG toolbox family and depends on:

**msg_coreg** — MSG Coregistration Toolbox (produces the input geometry files)  
https://github.com/maikeschmidt/msg_coreg

**msg_fwd** — MSG Forward Modelling Toolbox (runs BEM/FEM/Biot-Savart/sphere
on shifted geometries, and defines the metrics, percentiles and intervals
every msg toolbox reports)  
https://github.com/maikeschmidt/msg_fwd

Both must be cloned as sibling directories to `msg_pert`. `pt_add_functions`
warns if msg_fwd lacks a function msg_pert calls.

Related, not required:

**msg_analysis** — the real-data counterpart: concurrent spinal cord and
brain OPM recordings after median nerve stimulation, with source analysis on
msg_fwd forward models. The factors
msg_pert ranks in simulation (forward model, source and sensor error,
conductivity, sensor noise) are the ones to test there on measured N13
responses.  
https://github.com/maikeschmidt/msg_analysis

**msg_models** — numerical and mesh factors (convergence, bone conductivity
sweep, mesh repeatability). Its numbers use the same metrics and intervals,
so its mesh and bone-conductivity results can be read on the same scale as
the perturbations here. Its FEM seed-pair study measures the same re-meshing
effect as msg_pert's re-meshing floor, across mesh resolutions. Its bone-conductivity sweep varies one compartment
over a wide range; the conductivity perturbations here scale every
compartment by a few percent to 50% at once.  
https://github.com/maikeschmidt/msg_models

---

## Directory Structure

```
msg_pert/
├── config_pert_paths.m            — every folder (set data_root here)
├── config_pert.m                  — every other setting
├── run_all_perturbations.m        — PHASE 1 master: geometries + lead fields
├── run_all_pert_analysis.m        — PHASE 2 master: runs stages 0–6
├── RUN_ORDER.md                   — what to run, in what order, at what cost
├── INTERPRETATION.md              — how to read every output
│
├── pt_generate_source_shifts.m    — PHASE 1: 24 source-shift geometries
├── pt_generate_sensor_shifts.m    — PHASE 1: 24 sensor-shift geometries
│   (functions/pt_compute_leadfields.m — PHASE 1: BEM, Biot-Savart and
│    conductivity lead fields with msg_fwd's solvers)
│
├── pt_load_leadfields.m           — stage 0: load one array, scale each file, check
├── pt_compute_metrics.m           — stage 0: RE, r², RDM, lnMAG for every perturbation
├── pt_baseline.m                  — stage 1: noise-free baseline
├── pt_within_modality.m           — stages 2/3: within MSG / within ESG
├── pt_compare_modalities.m        — stage 4: MSG vs ESG
├── pt_noise_simulate.m            — stage 5a: perturbation + sensor noise
├── pt_noise_analyse.m             — stage 5b: figures, tables, intervals
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
│   ├── pt_compare_two.m           — paired / unpaired difference + bootstrap interval
│   ├── pt_stats_within.m, pt_stats_between.m, pt_compare_row.m
│   ├── pt_rank_corr.m             — Spearman rho + bootstrap interval
│   ├── pt_describe.m              — median, IQR, 95th percentile, intervals
│   ├── pt_quantile.m, pt_shift_summary.m, pt_metric_info.m, pt_sep_label.m
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

### Step 2: Set the paths

Every folder is in `config_pert_paths.m`. Usually only the first line
changes:

```matlab
data_root    = 'D:\Simulations\Pertubations';   % SET THIS
results_root = fullfile(data_root, 'results');
```

Put the original geometry from msg_coreg at
`<data_root>/geometries/<base_geom_name>.mat` (MSG) and
`<data_root>/geoms_elec/<base_geom_name>.mat` (ESG). Everything else is
created by the two run scripts. The file is not called `config_paths`
because msg_fwd's file of that name is on the path too.

### Step 3: Configure

Everything else is in `config_pert.m`. Sensor description and method
availability are declared **per modality**:

```matlab
base_geom_name = 'original';   % stem of the msg_coreg geometry file

mods_cfg.msg.sensor_n_axes = 3;    % 3 = triaxial MSG
mods_cfg.msg.sensor_is_meg = true;
mods_cfg.msg.have_bem      = true; % which forward models to compute and analyse
% ... and the same block for mods_cfg.esg (sensor_n_axes = 2)

pert_modalities = {'msg', 'esg'};  % e.g. {'msg'} to run MSG only
analysis_array  = 'back';          % the array computed and analysed
evoked_src_mm   = 60;              % C7, distance along the source space (mm)
```

`sensor_n_axes` and `sensor_is_meg` are **declared, not inferred**: an ESG
electrode count can also be divisible by 3, so guessing would mis-split the
leadfield.

### Step 4: Geometries and lead fields (Phase 1)

```matlab
run_all_perturbations;
```

For each modality this generates the 24 source-shift and 24 sensor-shift
geometries (3 bundles × 8 random shifts, plus the unshifted references),
saves the shift vectors next to them (`config_pert` reads them, so nothing
is pasted by hand), and computes with msg_fwd's solvers, in msg_fwd's file
layout:

- BEM on every geometry, and on 24 conductivity draws
- Biot-Savart on every geometry (MSG)
- FEM on every geometry, on the same 24 conductivity draws, and on 8
  re-meshes of the unshifted geometry (MSG; Windows with DUNEuro, as for
  msg_fwd's `run_fem_leadfields`; `duneuro_binpath` in `config_pert_paths`)

Existing files are skipped, so an interrupted run resumes. This is the slow
part: hours for the BEM, days for the FEM. Set `do_fem = false` at the top
to get the BEM and Biot-Savart results first.

Sphere lead fields are not part of Phase 1. To add them, run msg_fwd's
runner on the printed geometry list, write them to `msg_sphere_dir` and set
`have_sphere`.

### Step 5: Run the analysis (Phase 2)

```matlab
run_all_pert_analysis;
```

Flags at the top of that script run a subset of stages. Any stage can be
re-run alone once its inputs exist, see [RUN_ORDER.md](RUN_ORDER.md).
Outputs go to numbered folders under `stage_results_dir`:

| Folder | Contents |
|---|---|
| `1_baseline` | topoplots, amplitude and focality along the cord, model-type decomposition |
| `2_within_msg`, `3_within_esg` | decomposition per bundle, dose response, types and methods compared, anatomy contribution and FEM re-meshing floor (MSG), comparison tables |
| `4_msg_vs_esg` | headline comparison, difference maps, interaction, slopes, along-cord differences |
| `5_noise/c7`, `5_noise/cord` | noise alone, perturbation + noise, detectability, critical noise levels, error budget, systems compared per level — at the C7 source and over the cord |
| `6_summary` | one ranking of every source of error per modality, over the cord and at C7 |
| `7_cord_activity` | simulated N11 + N13, recordings per model, BEM–FEM divergence, critical noise levels and perturbation effects, all against N13 extent |

Every table is written as `.csv` (every orientation and axis set), `.txt`
(readable) and `.tex` (booktabs; headline orientation, whole array).

---

## Script Reference

| Script | Stage | Description |
|---|---|---|
| `run_all_perturbations` | Phase 1 | Master script: geometries and lead fields, both modalities |
| `pt_generate_source_shifts` | Phase 1 | 24 geometries: 3 bundles × 8 random source-space shifts; saves the shift vectors |
| `pt_generate_sensor_shifts` | Phase 1 | 24 geometries: 3 bundles × 8 random sensor-array shifts; saves the shift vectors |
| `pt_compute_leadfields` | Phase 1 | BEM, FEM, Biot-Savart, conductivity-perturbed BEM/FEM and FEM re-meshing lead fields with msg_fwd's solvers, in msg_fwd's file layout |
| `run_all_pert_analysis` | Phase 2 | Master script: stages 0–6, both modalities |
| `pt_load_leadfields` | 0 | Load the analysed array for every method and perturbation; per-file unit scaling; magnitude check against each reference |
| `pt_compute_metrics` | 0 | RE, r², RDM and lnMAG per source, method, orientation and axis set → `pert_metrics.mat` |
| `pt_baseline` | 1 | Unperturbed topoplots, amplitude and focality along the cord, forward-model-type decomposition |
| `pt_within_modality` | 2, 3 | Within one modality: descriptive table, decomposition per bundle, dose response, type and method comparisons, along-cord source-vs-sensor differences, anatomy contribution (BEM and FEM), FEM re-meshing floor |
| `pt_compare_modalities` | 4 | MSG vs ESG: descriptive, differences with intervals, interaction, dose-response slopes, along-cord differences, difference maps |
| `pt_noise_simulate` | 5 | An N13-like evoked response at C7 (and at every cord position) through every perturbed field and trial-averaged sensor noise, scored against the true field |
| `pt_noise_analyse` | 5 | Noise alone, perturbation + noise, detectability, critical noise level, error budget, system comparisons per level (including the Biot-Savart bound), within-modality comparisons under noise — once at C7, once over the cord |
| `pt_summary_table` | 6 | Every source of error ranked on one scale, over the cord and at C7 |
| `pt_cord_activity` | 7 | The N11 volley and N13 as a sensors × time recording through every model and system, N13 extent swept; BEM vs FEM divergence, detectability and perturbation sensitivity |
| `pt_modality` | — | Get/set the active modality; survives each script's `clearvars` |
| `pt_run_step` | — | Run one script in an isolated workspace |
| `pt_diagnose_leadfields` | — | Inspect a loaded lead-field set |

---

## Standalone noise simulation (`simulations/`)

The `simulations/` package predates stage 5 and remains usable on its own. It
scores the noisy sensor-by-time data matrix with r², for a single
representative shift per bundle. Stage 5 of the main pipeline asks the same
question for every realisation with all four metrics and intervals, so
use stage 5 for results that sit alongside stages 1–4. The noise floors in
`config_pert.m` are copied from `simulations/config_sim.m`; keep the two in
step if you change either. The source waveform differs on purpose: stage 5
uses an N13-like response at 13 ms, the package its own 70 Hz burst. See
`simulations/README.md`.

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
amplitude part and a topography part. One perturbation realisation is one
observation (its median over the cord, or its value at C7 in stage 5).
Because every value is computed exactly from a simulated geometry, results
are reported as msg_fwd reports its replicate geometries: medians and
differences with percentile bootstrap intervals, the IQR, and the 95th
percentile with its interval — not significance tests. Percentiles and
intervals come from msg_fwd's `pctl` and `st_boot_ci_median`. See
[INTERPRETATION.md](INTERPRETATION.md).

**Reference field.** Forward models (Biot-Savart, FEM) are compared against
the unperturbed MRI-derived BEM, as in msg_fwd and msg_models. Perturbations are measured
against each model's own unperturbed field. Biot-Savart, which has no
anatomy, is the anatomy-independent comparison: the geometric part of the
sensitivity, and the upper bound on recoverability through noise.

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

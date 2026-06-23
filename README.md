# msg_pert — MSG Perturbation Analysis Toolbox

**Systematic perturbation analysis for Magnetospinography (MSG) forward models:
generating anatomically shifted geometry files and quantifying the sensitivity
of BEM, FEM, Biot-Savart, and single-sphere forward models to perturbations
in source space and sensor array position.**

Developed by **Maike Schmidt** at the **Department of Imaging Neuroscience,
University College London**.

---

## Overview

This toolbox answers two questions about MSG forward models:

> *Which forward model is most sensitive to imperfect anatomical knowledge?*
> *How large a shift in anatomy still produces consistent source reconstruction?*

Starting from an original geometry file produced by `msg_coreg`, msg_pert
generates two sets of systematically shifted geometry files and provides the
tools to analyse the resulting forward fields once they have been computed
in `msg_fwd`.

### Perturbation types

| Type | What shifts | What stays fixed | Configurations |
|---|---|---|---|
| **Source space** | Spinal cord mesh, bone mesh, source positions | Torso, heart, lungs | 18: ±2, ±4, ±6 mm × X, Y, Z |
| **Sensor array** | Sensor coil positions | All anatomy meshes | 24: 3 bundles × 8 random shifts |

**Source space shifts** model uncertainty in anatomical understanding — e.g. the cord
sits 4 mm more anterior than segmented. The entire cord-bone unit moves together.

**Sensor array shifts** model registration error — e.g. the scanner cast was placed
2–10 mm off its true position. Three error bundles represent small (~2 mm),
medium (~5 mm), and large (~10 mm) registration uncertainty.

---

## Two-phase workflow

```
Phase 1 — Geometry generation

  msg_coreg geometry file
         │
         ▼
  pt_generate_source_shifts   →  18 geometry .mat files
  pt_generate_sensor_shifts   →  24 geometry .mat files  +  shift vectors (paste into config_pert)
         │
         │  (Copy printed filename list into msg_fwd)
         ▼
  msg_fwd: run BEM, FEM, Biot-Savart, single sphere
  on each shifted geometry (front + back arrays)
         │
         ▼
  msg_fwd: load_and_organise_leadfields
  → leadfields_organised.mat


Phase 2 — Analysis  (run_perturbation_analysis)

  leadfields_organised.mat
         │
         ├─ pt_compute_rsq          r² per source, per orientation, per model
         ├─ pt_plot_curves          r² vs cord distance  (source and sensor mode)
         ├─ pt_plot_displacement    displacement vs r²   (sensor mode)
         └─ pt_compute_table        summary tables (.txt and .csv)
```

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
├── pert_path.m                    — path function (locates repository root)
├── pt_add_functions.m             — dependency setup (msg_coreg, msg_fwd, HBF)
├── config_pert.m                  — all paths and perturbation parameters
│
├── pt_generate_source_shifts.m    — PHASE 1: generate 18 source-shift geometries
├── pt_generate_sensor_shifts.m    — PHASE 1: generate 24 sensor-shift geometries
│
├── run_perturbation_analysis.m    — PHASE 2: master analysis script
├── pt_compute_rsq.m               — compute per-source r² (source and sensor modes)
├── pt_plot_curves.m               — r² vs cord distance figures
├── pt_plot_displacement.m         — median displacement vs r² (sensor mode)
├── pt_compute_table.m             — summary tables (.txt and .csv)
│
└── README.md
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

### Step 2: Configure paths

Edit `config_pert.m` and set the four path variables and `base_geom_name`:

```matlab
geoms_path           = 'D:\my_study\geometries';       % original geometry .mat
perturbed_geoms_path = 'D:\my_study\pert_geometries';  % output for shifted files
forward_fields_base  = 'D:\my_study\leadfields';       % msg_fwd leadfield output
save_base_dir        = 'D:\my_study\figures';          % figures and tables
base_geom_name       = 'geometries_sub001_experimental';
```

### Step 3: Generate source-shift geometries (Phase 1a)

```matlab
pt_generate_source_shifts;
```

This produces 18 shifted geometry files in `perturbed_geoms_path` and prints a
filename list to paste into `msg_fwd`.

### Step 4: Generate sensor-shift geometries (Phase 1b)

```matlab
pt_generate_sensor_shifts;
```

This produces 24 shifted geometry files and prints:
1. The shift vectors — paste into `config_pert.m` under `sensor_shift_vectors`
2. A filename list for `msg_fwd`

### Step 5: Run forward models in msg_fwd

Copy the printed filename lists into `msg_fwd` (into `run_bem_leadfields.m` and the
Biot-Savart / single-sphere scripts) and run all four forward model types for both
front and back arrays. Then run `load_and_organise_leadfields.m` in msg_fwd to
produce `leadfields_organised.mat`.

### Step 6: Run perturbation analysis (Phase 2)

```matlab
run_perturbation_analysis;
```

Or run individual analysis scripts standalone:

```matlab
pt_compute_rsq;            % compute r² (run first)
pt_plot_curves;            % r² vs cord distance
pt_plot_displacement;      % displacement vs r² (sensor mode)
pt_compute_table;          % summary tables
```

---

## Script Reference

| Script | Phase | Description |
|---|---|---|
| `pert_path` | — | Returns the absolute path to the msg_pert root directory |
| `pt_add_functions` | — | Adds msg_coreg, msg_fwd, HBF, and FieldTrip wrappers to the MATLAB path |
| `config_pert` | — | Shared configuration: paths, shift parameters, naming conventions, plot styling |
| `pt_generate_source_shifts` | 1 | Generate 18 geometry files for ±2/4/6 mm source-space shifts along X/Y/Z |
| `pt_generate_sensor_shifts` | 1 | Generate 24 geometry files for 3 bundles × 8 random sensor-array shifts |
| `run_perturbation_analysis` | 2 | Master script: runs all analysis steps in order |
| `pt_compute_rsq` | 2 | Compute per-source r² between each shifted model and the original |
| `pt_plot_curves` | 2 | r² vs cord distance figures (individual shift axes / bundles + overviews) |
| `pt_plot_displacement` | 2 | Median displacement vs r² scatter plots (sensor mode; requires shift vectors in config) |
| `pt_compute_table` | 2 | Write median r², min r², and first-drop thresholds as .txt and .csv |

---

## Metrics

**r² (squared Pearson correlation)**

Computed per source position by comparing the full leadfield vector of the
shifted model against the unshifted original:

```
r² = (Pearson r)^2 between shifted and original leadfield at each source
r² = 1.0 — identical leadfields (no effect of perturbation)
r² = 0.0 — no correlation
```

Computed separately per dipole orientation (VD / RC / LR) and sensor axis.
Edges (first and last source) are excluded.

**Threshold conventions:**
- r² < 0.99 — first position where the perturbation has a measurable effect
- r² < 0.95 — first position where the effect is practically significant

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

If you use this toolbox, please cite:

> Schmidt, M. et al. (2026). *Forward model sensitivity in Magnetospinography.*
> [Journal TBC] [DOI TBC]

Please also cite the companion toolboxes:

> msg_coreg: https://github.com/maikeschmidt/msg_coreg  
> msg_fwd:   https://github.com/maikeschmidt/msg_fwd

---

## Contact

For questions, issues, or contributions, open an issue or pull request on GitHub.  
Contact: maike.schmidt.23@ucl.ac.uk

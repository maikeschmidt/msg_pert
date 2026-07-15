# msg_pert — MSG Perturbation Analysis Toolbox

**Systematic perturbation analysis for Magnetospinography (MSG) forward models:
generating anatomically shifted geometry files and perturbed tissue
conductivities, then quantifying the sensitivity of BEM, FEM, Biot-Savart, and
single-sphere forward models to perturbations in source space, sensor array
position, and tissue conductivity. Includes a self-contained `simulations/`
package that adds realistic sensor noise to evoked responses.**

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

## Two-phase workflow

```
Phase 1 — Perturbation generation

  msg_coreg geometry file
         │
         ▼
  pt_generate_source_shifts        →  24 source-shift geometry .mat files
  pt_generate_sensor_shifts        →  24 sensor-shift geometry .mat files
                                       (+ shift vectors to paste into config_pert)
         │
         │  (run forward models in msg_fwd)
         ▼
  msg_fwd: run_bem_leadfields / run_fem_leadfields / Biot-Savart / sphere
           on each shifted geometry (front + back arrays)
  msg_fwd: run_conductivity_perturbation  →  24 BEM conductivity leadfields
         │
         ▼
  leadfield .mat files (per geometry, per method, per array)


Phase 2 — Analysis  (run_perturbation_analysis)

         ├─ pt_load_leadfields      load BEM/FEM/BS/sphere + BEM-cond leadfields
         │                          for all perturbations → leadfields_organised.mat
  leadfields_organised.mat
         │
         ├─ pt_compute_rsq          r² per source/orientation/method (source, sensor, cond)
         ├─ pt_plot_curves          r² vs cord distance (detail, summary, cross-model)
         ├─ pt_plot_heatmaps        pairwise RE and r² heatmaps
         ├─ pt_plot_displacement    displacement / % change vs r²  (cervical + full cord)
         ├─ pt_plot_slope_vs_position   slope of r² change vs cord position
         └─ pt_compute_table        summary tables (.txt and .csv)


Optional — Realistic-measurement simulation  (simulations/, see below)

         run_simulation_analysis    evoked response + sensor noise across
                                    systems (SQUID/OP-MSG/ESG) and geometries
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
├── pt_generate_source_shifts.m    — PHASE 1: generate 24 source-shift geometries
├── pt_generate_sensor_shifts.m    — PHASE 1: generate 24 sensor-shift geometries
│                                    (conductivity leadfields come from
│                                     msg_fwd/run_conductivity_perturbation.m)
│
├── run_perturbation_analysis.m    — PHASE 2: master analysis script
├── pt_load_leadfields.m           — load + organise all leadfields (incl. BEM-cond)
├── pt_compute_rsq.m               — per-source r² (source, sensor, conductivity)
├── pt_plot_curves.m               — r² vs cord distance figures
├── pt_plot_heatmaps.m             — pairwise RE and r² heatmaps
├── pt_plot_displacement.m         — displacement / % change vs r²
├── pt_plot_slope_vs_position.m    — slope of r² change vs cord position
├── pt_compute_table.m             — summary tables (.txt and .csv)
├── pt_compare_perturbations.m     — cross-perturbation comparison utilities
│
├── simulations/                   — self-contained realistic-measurement package
│   ├── run_simulation_analysis.m  — master script (6 steps)
│   ├── config_sim.m               — models, systems, noise floors, geometry list
│   ├── sim_run_geometries.m       — evoked response + noise, looped over geometries
│   ├── sim_plot_original.m        — base noise curves for the original geometry
│   ├── sim_plot_comparison.m      — noise curves overlaid across geometry variants
│   ├── sim_plot_worstcase.m       — systems compared under the largest shift
│   ├── sim_plot_topoplots.m       — perfect-field topoplots per model
│   ├── sim_plot_noise_topoplot.m  — measured (noisy) topoplot at a chosen source
│   ├── functions/                 — sim helpers (load, positions, closed-form r²)
│   └── README.md
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

### Step 6: Run perturbation analysis (Phase 2)

Before running, open `pt_load_leadfields.m` and set the `have_<method>` flags
to match the forward models you computed in msg_fwd:

```matlab
have_bem      = true;    % BEM via Helsinki BEM Framework
have_fem      = false;   % FEM via DUNEuro
have_bslaw    = false;   % Biot-Savart (infinite space)
have_sphere   = false;   % Single sphere (Sarvas analytical)
have_bem_cond = false;   % BEM with perturbed tissue conductivities
```

Also set `sensor_n_axes` / `sensor_is_meg` (3/true for MSG, 2/false for ESG) —
these are declared, not inferred, because an ESG electrode count can be
divisible by 3 and would otherwise be mis-split.

Then run the full pipeline:

```matlab
run_perturbation_analysis;
```

Or run individual steps standalone:

```matlab
pt_load_leadfields;          % load and organise leadfields (run first)
pt_compute_rsq;              % r² for source, sensor, conductivity
pt_plot_curves;              % r² vs cord distance
pt_plot_heatmaps;           % pairwise RE / r² heatmaps
pt_plot_displacement;       % displacement / % change vs r²
pt_plot_slope_vs_position;  % slope of r² vs cord position
pt_compute_table;           % summary tables
```

---

## Script Reference

| Script | Phase | Description |
|---|---|---|
| `pert_path` | — | Returns the absolute path to the msg_pert root directory |
| `pt_add_functions` | — | Adds msg_pert/functions, msg_coreg, msg_fwd, HBF, and FieldTrip wrappers to the MATLAB path |
| `config_pert` | — | Shared configuration: paths, source/sensor/conductivity parameters, naming, plot styling |
| `pt_generate_source_shifts` | 1 | Generate 24 geometry files for 3 bundles × 8 random source-space shifts (~2/5/10 mm) |
| `pt_generate_sensor_shifts` | 1 | Generate 24 geometry files for 3 bundles × 8 random sensor-array shifts |
| `run_perturbation_analysis` | 2 | Master script: runs all analysis steps in order |
| `pt_load_leadfields` | 2 | Load and organise BEM/FEM/BS/sphere + BEM-conductivity leadfields; saves `leadfields_organised.mat` |
| `pt_compute_rsq` | 2 | Compute per-source r² for source, sensor, and conductivity perturbations vs the original |
| `pt_plot_curves` | 2 | r² vs cord distance figures (detail, summary, cross-model) for all three modes |
| `pt_plot_heatmaps` | 2 | Pairwise RE and r² heatmaps (within- and cross-method) for all three modes |
| `pt_plot_displacement` | 2 | Displacement (mm) or % conductivity change vs r² (individual: cervical; combined + trend table: full cord) |
| `pt_plot_slope_vs_position` | 2 | Slope of r² change vs cord position, from the displacement trend tables |
| `pt_compute_table` | 2 | Write median r², min r², and first-drop thresholds as .txt and .csv |
| `pt_compare_perturbations` | 2 | Cross-perturbation comparison utilities |

---

## Realistic-measurement simulation (`simulations/`)

The perturbation pipeline above asks how the *noise-free* forward field changes
when the model is wrong. The self-contained `simulations/` package asks the
complementary question: given the forward field, how much does realistic
**sensor noise** degrade a measured evoked response, and how does that differ
between sensor systems and between geometry variants?

It simulates a Gaussian-windowed evoked burst at every source on the cord,
projects it through a chosen leadfield, adds trial-averaged white sensor noise
across a sweep of levels, and scores r² against the noise-free field. Three
systems are compared — **SQUID MSG**, **OP-MSG**, and **ESG** — each with its own
published noise floor and measurement bandwidth. Because MSG (fT) and ESG (µV)
noise floors are not comparable in absolute terms, everything is expressed as a
multiple of each system's own baseline.

Run `simulations/run_simulation_analysis` (configure `simulations/config_sim.m`
first). Outputs: perfect-field topoplots, per-geometry noise curves, cross-variant
comparison curves (mean across cord + IQR band), a worst-case system comparison
(largest shift, scored against the original field), and measured noisy topoplots.
See `simulations/README.md` for full detail.

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

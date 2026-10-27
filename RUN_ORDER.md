# Run order

What to run, in what order, and roughly what each stage costs.

See [README.md](README.md) for what the toolbox does and
[INTERPRETATION.md](INTERPRETATION.md) for how to read the outputs.

---

## Before you start

### Platform

Geometry generation (Phase 1) needs msg_coreg with the Helsinki BEM
Framework, and the lead fields come from msg_fwd (BEM, and FEM on Windows via
DUNEuro). The analysis (Phase 2) needs only MATLAB and msg_fwd's `functions/`
folder, so it runs on any machine once the lead fields exist.

Two test suites run anywhere MATLAB does, with no toolbox dependencies:

```bash
matlab -batch "cd msg_pert/tests; test_pt_core"
```

```bash
matlab -batch "cd msg_pert/tests; test_pt_pipeline"
```

`test_pt_core` checks the metric and statistics core against known answers.
`test_pt_pipeline` builds a small synthetic data set and runs every stage on
it (about 15 minutes), so a broken stage shows up before a long real run.
Run both before and after any change to `functions/` or a stage script.

### Configuration

Two files, split as msg_fwd splits `config_paths` from `config_models`:

- **`config_pert_paths.m`** — every folder. Set `data_root` (and
  `results_root` if results go elsewhere); the layout under it is listed at
  the top of the file.
- **`config_pert.m`** — everything else:

| Block | Holds |
|---|---|
| per-modality (`mods_cfg.msg`, `mods_cfg.esg`) | which methods to compute and analyse, sensor axes and how they correspond between modalities (folders come from `config_pert_paths`) |
| `analysis_array`, `unit_scale_mode` | which array is analysed; how lead fields are brought to fT/nAm and µV/nAm |
| source / sensor / conductivity | the perturbation definitions (shift vectors read from the files Phase 1 saves) |
| staged analysis | output folder, orientations, interval settings (`stats_n_boot`, `stats_ci_level`, `stats_cover_pct`), pairing |
| stage 1 | which models are shown and which are compared |
| FEM meshing | `fem_maxvol_mm3`, `fem_mesh_seed`, `fem_remesh_seeds`, `fem_cond_nominal` (msg_fwd's production values) |
| stage 5 | sensor systems (BEM, FEM and Biot-Savart), C7 position (`evoked_src_mm`), noise floors, bandwidths, N13 waveform, noise sweep |

Check these before a first run:

- **`msg_esg_paired`** says, per perturbation type, whether shift *k* in MSG
  is the same geometric change as shift *k* in ESG. It decides whether the
  MSG-vs-ESG tests are paired. ESG sensor shifts hold Z at 0, so they are
  not the same shift.
- **The shift vectors** give the perturbation sizes on every dose-response
  axis. Phase 1 saves the ones it applied (`pert_source_shift_vectors.mat`,
  `pert_sensor_shift_vectors.mat` next to the geometries) and `config_pert`
  reads them; the pasted block in `config_pert` is only the fallback for
  geometries generated before that.
- **`evoked_src_mm`** must be C7's distance along the source space (source
  index × `src_spacing_mm`). Stage 5 prints the source it picked; check it
  against the vertebral levels of your geometry.

---

## Phase 1 — perturbed geometries and lead fields (once per study)

```matlab
run_all_perturbations
```

For each modality in `pert_modalities`:

| Step | Does | Output (`config_pert_paths`) | Needed for |
|---|---|---|---|
| 1 | `pt_generate_source_shifts`: 24 source-space geometries + reference, shift vectors saved | `msg_geoms_dir` / `esg_geoms_dir` | everything |
| 2 | `pt_generate_sensor_shifts`: 24 sensor-array geometries + reference | same | everything |
| 3 | BEM on all 50 geometries | `*_bem_dir` | every stage |
| 4 | Biot-Savart on all 50 (MSG only) | `msg_bslaw_dir` | stage 1 model-type comparison, stage 2 method comparison and anatomy contribution, stage 5 upper-bound systems |
| 5 | BEM with 24 conductivity draws on the original geometry | `*_bem_cond_dir` | conductivity rows everywhere |
| 6 | FEM on all 50 (MSG only) | `msg_fem_dir` | stage 1 FEM vs BEM, FEM rows in stages 2 and 6, stage 5 FEM systems |
| 7 | FEM with the same 24 conductivity draws (MSG only) | `msg_fem_dir` | FEM conductivity rows, BEM vs FEM conductivity sensitivity |
| 8 | FEM on the original geometry re-meshed with `fem_remesh_seeds` (MSG only) | `msg_fem_dir` | the re-meshing floor (stage 2 section I, a stage-6 row) |

The lead fields come from `pt_compute_leadfields`, which performs the steps
of msg_fwd's `run_bem_leadfields`, `run_biot_savart_leadfields`,
`run_conductivity_perturbation` and `run_fem_leadfields` with msg_fwd's own
functions (`bem_build_boundaries`, `bs_leadfield`, `fem_build_volume_mesh`,
`convert_duneuro_to_fieldtrip`, `log_compute_time`) and writes msg_fwd's
file layout. Those runners are scripts with their own folders and
geometry lists; this does the same work with msg_pert's. Only
`analysis_array` is computed. A method whose `have_*` flag is false is not
computed. Every step skips files that already exist, so an interrupted run
resumes, and flags at the top run a subset.

**FEM meshes.** Volume meshes are cached in `msg_fem_dir/meshes`, one per
anatomy and seed. Sensor shifts and conductivity draws are solved on the
mesh of the unshifted geometry, so their FEM differences contain no
re-meshing; source shifts get their own mesh, all with `fem_mesh_seed`.
Step 8 measures what re-meshing alone does. FEM needs Windows, the DUNEuro
binary (`duneuro_binpath` in `config_pert_paths`), ISO2Mesh and
`fem_calc_fwds` on the path, exactly as msg_fwd's `run_fem_leadfields`; the
FEM steps stop with a message if any is missing, and the others still run.
FEM is MSG only: the DUNEuro call msg_fwd uses has no EEG route.

Cost: 74 BEM head models per modality (hours); 33 FEM meshes and 82 DUNEuro
solves for MSG (days). Biot-Savart is seconds. `do_fem = false` at the top
runs everything else first.

Sphere is not part of Phase 1. To add it, run msg_fwd's runner on the
geometry list Phase 1 prints, put the files in `msg_sphere_dir`, and set
`have_sphere`. Do **not** run msg_fwd's `load_and_organise_leadfields`:
msg_pert organises the files itself in stage 0.

---

## Phase 2 — analysis

```matlab
run_all_pert_analysis
```

runs every stage in order. The flags at its top run a subset; each stage
reads only files written by earlier ones, so any stage can be re-run alone
once its inputs exist.

| Stage | Script(s) | Reads | Time* | Output folder |
|---|---|---|---|---|
| 0 | `pt_load_leadfields`, `pt_compute_metrics` (per modality) | lead-field files | minutes | `<forward_fields_base>` |
| 1 | `pt_baseline` | `leadfields_organised.mat` | ~1 min | `1_baseline` |
| 2 | `pt_within_modality` (MSG) | `pert_metrics.mat` | ~5 min | `2_within_msg` |
| 3 | `pt_within_modality` (ESG) | `pert_metrics.mat` | ~5 min | `3_within_esg` |
| 4 | `pt_compare_modalities` | both `pert_metrics.mat` | ~5 min | `4_msg_vs_esg` |
| 5 | `pt_noise_simulate`, then `pt_noise_analyse` | lead fields + metrics | 20–90 min | `5_noise/c7`, `5_noise/cord` |
| 6 | `pt_summary_table` | the stage `*_long.csv` files | seconds | `6_summary` |
| 7 | `pt_cord_activity` | lead fields + metrics | 10–30 min | `7_cord_activity` |

\*Rough times for 110 sources and ~100 sensors per axis. Stage 5 scales with
`noise_n_real` × number of noise levels × number of systems (seven by
default: three BEM, two FEM, two Biot-Savart); drop `noise_n_real` to 5, or
remove systems from `noise_systems`, for a quick look.

### Running one stage by hand

Stages 0, 2 and 3 work on one modality at a time, chosen with
`pt_modality`:

```matlab
pt_modality('set', 'esg');
pt_compute_metrics
pt_within_modality
```

Stages 1, 4, 5 and 6 use both modalities and need no selection.

### After changing anything

| You changed | Re-run from |
|---|---|
| `data_root` or any folder in `config_pert_paths` | nothing, if the files moved with it; otherwise Phase 1 |
| lead-field files, `analysis_array`, `unit_scale_mode` | stage 0 |
| shift vectors, `metric_orientations`, `metric_defaults` (msg_fwd) | stage 0 (`pt_compute_metrics`) |
| interval settings, pairing flags, figure options | the stage concerned, then 6 |
| noise systems, `evoked_src_mm`, floors, waveform, sweep | stage 5, then 6, then 7 |
| `cord_*` (N11 volley, N13 extents) | stage 7 |

### Stage 0 — check the scale log

`pt_load_leadfields` writes `leadfield_scale_log.csv` next to
`leadfields_organised.mat`: one row per file, with the factor applied and
why. It also warns if any perturbed lead field is more than 10× larger or
smaller than its reference. Read the log once after a first load. A
perturbation cannot change field strength by an order of magnitude, so such
a warning always means a units problem, and every RE and gain value would be
wrong while r² still looked normal.

---

## Dependency map

```
run_all_perturbations: geometries ──> lead fields (msg_fwd solvers)
                              │
                              v
        stage 0: pt_load_leadfields ──> pt_compute_metrics      (per modality)
                   │                          │
                   ├──> stage 1 pt_baseline   ├──> stage 2/3 pt_within_modality
                   │                          ├──> stage 4 pt_compare_modalities
                   └──────────────┬───────────┘
                                  v
                   stage 5 pt_noise_simulate ──> pt_noise_analyse
                                  │
   stage 1, 2, 3, 5 *_long.csv ───┴──────────> stage 6 pt_summary_table
```

---

## Earlier r²-only outputs

The scripts from before the staged pipeline (`pt_compute_rsq`,
`pt_plot_curves`, `pt_plot_heatmaps`, `pt_plot_displacement`,
`pt_plot_slope_vs_position`, `pt_compute_table`, `pt_compare_perturbations`)
still work and read the same loaded lead fields. Set `run_legacy = true` in
`run_all_pert_analysis` to produce them. They report r² only, and
`pt_compare_perturbations` still uses permutation tests and needs the
Statistics Toolbox; the staged pipeline replaces them with all four metrics
and the interval-based summaries shared with msg_fwd.

---

## If time is tight

1. **Stages 0, 2, 3, 4** — the noise-free comparison, which is the core.
2. **Stage 1** — needed to read everything else against (the model-type
   difference and the field amplitude along the cord).
3. **Stage 5** with `noise_n_real = 5` — the qualitative noise picture.
4. **Stage 5** at full `noise_n_real`, then **stage 6**.

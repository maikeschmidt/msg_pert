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

Everything is set in one file, `config_pert.m`:

| Block | Holds |
|---|---|
| per-modality (`mods_cfg.msg`, `mods_cfg.esg`) | geometry and lead-field paths, which methods exist, sensor axes and how they correspond between modalities |
| `analysis_array`, `unit_scale_mode` | which array is analysed; how lead fields are brought to fT/nAm and µV/nAm |
| source / sensor / conductivity | the perturbation definitions (shift vectors pasted from Phase 1) |
| staged analysis | output folder, orientations, statistics settings, pairing |
| stage 1 | which models are shown and which are compared |
| stage 5 | sensor systems, noise floors, bandwidths, source waveform, noise sweep |

Check two things before a first run:

- **`msg_esg_paired`** says, per perturbation type, whether shift *k* in MSG
  is the same geometric change as shift *k* in ESG. It decides whether the
  MSG-vs-ESG tests are paired. ESG sensor shifts hold Z at 0, so they are
  not the same shift.
- **The shift vectors** in the source and sensor blocks must be the ones
  Phase 1 actually applied: they give the perturbation sizes on every
  dose-response axis.

---

## Phase 1 — perturbed geometries and lead fields (once per study)

```matlab
pt_add_functions
pt_generate_source_shifts     % 24 source-space geometries
pt_generate_sensor_shifts     % 24 sensor-array geometries
```

Paste the printed shift vectors into `config_pert.m`. Then, in msg_fwd, run
the forward models on every printed geometry for the array set in
`analysis_array`:

| msg_fwd script | Gives | Needed for |
|---|---|---|
| `run_bem_leadfields` | BEM, MSG and ESG | every stage |
| Biot-Savart | MSG only | stage 1 model-type comparison, stage 2 method comparison |
| `run_conductivity_perturbation` | BEM with 24 conductivity draws, MSG and ESG | conductivity rows everywhere |
| `run_fem_leadfields` (optional) | FEM | an extra method in stages 2–3 |

Set the `have_<method>` flags in `config_pert.m` to match what you computed.
msg_pert organises the files itself: do **not** run msg_fwd's
`load_and_organise_leadfields` here.

---

## Phase 2 — analysis

```matlab
run_perturbation_analysis
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
| 5 | `pt_noise_simulate`, then `pt_noise_analyse` | lead fields + metrics | 15–60 min | `5_noise` |
| 6 | `pt_summary_table` | the stage `*_long.csv` files | seconds | `6_summary` |

\*Rough times for 110 sources and ~100 sensors per axis. Stage 5 scales with
`noise_n_real` × number of noise levels × number of systems; drop
`noise_n_real` to 5 for a quick look.

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
| lead-field files, paths, `analysis_array`, `unit_scale_mode` | stage 0 |
| shift vectors, `metric_orientations`, `metric_defaults` (msg_fwd) | stage 0 (`pt_compute_metrics`) |
| statistics settings, pairing flags, figure options | the stage concerned, then 6 |
| noise systems, floors, waveform, sweep | stage 5, then 6 |

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
Phase 1 geometries ──> msg_fwd lead fields
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
`run_perturbation_analysis` to produce them. They report r² only; the staged
pipeline replaces them with all four metrics and one consistent set of
statistics.

---

## If time is tight

1. **Stages 0, 2, 3, 4** — the noise-free comparison, which is the core.
2. **Stage 1** — needed to read everything else against (the model-type
   difference and the field amplitude along the cord).
3. **Stage 5** with `noise_n_real = 5` — the qualitative noise picture.
4. **Stage 5** at full `noise_n_real`, then **stage 6**.

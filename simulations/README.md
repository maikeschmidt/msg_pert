# msg_pert / simulations

Self-contained simulation package: how much does realistic **sensor noise**
degrade a measured evoked response, across sensor systems (SQUID-MSG, OP-MSG,
ESG) and across geometry variants (the unperturbed model vs representative
source / sensor / conductivity perturbations)?

This is separate from the `pt_*` perturbation pipeline in the parent folder,
which asks how the *noise-free* forward field changes when the model is wrong.
Here the model is taken as given and **noise** is the variable.

## Pipeline

Run `run_simulation_analysis` (or the steps individually):

| Step | Script | Produces |
|------|--------|----------|
| 1 | `sim_plot_topoplots` | Noise-free forward fields (Biot-Savart MSG, BEM MSG, BEM ESG) for the chosen geometry variant |
| 2 | `sim_run_geometries` | For every variant in `sim_geometries` × every system: evoked response + trial-averaged noise, scored by r² vs the noise-free field. One output subfolder per variant. |
| 3 | `sim_plot_original` | Base figures for the **original** geometry, all systems overlaid: median-across-cord (IQR band), single chosen source (±1 s.d.), and r² along the cord at each noise level |
| 4 | `sim_plot_comparison` | Overlaid r²-vs-noise curves across variants, one figure per system (rows = source/sensor/cond, lines = baseline/small/medium/large) |
| 5 | `sim_plot_noise_topoplot` | What one system actually measures at a chosen source + noise level |

## Configure `config_sim.m`

- **Leadfield roots + geometry stems** — where the `.mat` files live.
- **`sim_geometries`** — the list looped over. Pick one representative shift per
  bundle via `src_reps` / `sen_reps` / `cond_reps`. To compare noise levels for a
  single geometry only, leave one entry in the list.
- **`bem_patched` / unit scales** — a wrong leadfield scale flattens every curve.
  `sim_run_geometries` prints `peak|g|` vs `sigma@1x` each run; if they differ by
  orders of magnitude, fix the scale, not the noise sweep.
- **Noise floors, bandwidths, trials, evoked burst** — the physical model.

## Modality coverage

Source, sensor, and conductivity variants all exist for both MSG and ESG. Source
and sensor shifts load from each modality's `model.root`; conductivity loads from
each modality's `model.cond_root` (`msg_cond_root` / `esg_cond_root`). The loop
skips any (variant, system) whose leadfield file is absent, so a partially
generated set still runs without special-casing.

## Notes

- MSG (fT) and ESG (µV) noise floors are not comparable in absolute terms;
  everything is expressed as a multiple of each system's own baseline.
- The ESG floor (~1 µV/√Hz) is an amplifier-noise estimate; real ESG is
  dominated by much larger cardiac artefact, so ESG curves are an **optimistic**
  bound.
- Trial averaging is applied analytically (σ → σ/√N); for independent Gaussian
  noise this is exact, not an approximation.

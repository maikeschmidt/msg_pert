# Interpreting the results

How to read every output of the staged pipeline: what each number means, what
a sensible result looks like, and what points to a problem rather than a
finding.

See [README.md](README.md) for what the toolbox does and
[RUN_ORDER.md](RUN_ORDER.md) for what to run in what order.

---

## 1. The metrics

Every comparison in every stage reports the same four numbers. They are
defined once, in msg_fwd's `functions/lf_metrics.m`; msg_pert's vectorised
`pt_metrics_block` reproduces them exactly (checked by `test_pt_core`), so a
msg_pert number and a msg_fwd number mean the same thing.

| Metric | Formula | Sensitive to | Range |
|---|---|---|---|
| **RE** | `‖L₂−L₁‖₂ / ‖L₁‖₂ × 100` | magnitude **and** shape | 0 → ∞ % |
| **r²** | `(Pearson r)²` | shape only | 0 → 1 |
| **RDM** | `‖L̂₂−L̂₁‖` on unit-normalised fields | shape only | 0 → 2 |
| **gain** | `(e^lnMAG − 1) × 100`, lnMAG = `ln(‖L₂‖/‖L₁‖)` | magnitude only | −100 → ∞ % |
| **\|gain\|** | absolute value of gain | size of the magnitude change | 0 → ∞ % |

**The reference is always the unperturbed geometry** (`L₁`), and RE is
asymmetric: it is measured relative to the reference field.

**The decomposition** `RE ≈ √(gain² + (RDM×100)²)` says *why* a field
changed. RE carried mostly by gain means the field got stronger or weaker
but kept its shape. RE carried by RDM means the pattern moved. The first
biases amplitude estimates; the second also biases localisation.

**Signed gain is descriptive only.** A source shift towards the array raises
the field and one away lowers it, so across random shifts the signed gain can
average to about zero even when every shift changes the amplitude a lot. All
tests on amplitude therefore use |gain|.

**Orientations and axis sets.** Every metric is computed per dipole
orientation (VD, RC, LR) and for `ALL`, the three stacked into one vector,
which is the single number to quote. It is also computed per sensor axis and
for the whole array, all axes stacked. The headline figures and LaTeX tables
use `headline_orientation` (default `ALL`) on the whole array; the CSVs hold
everything.

---

## 2. The statistics

### The unit of observation

Every test treats **one perturbation realisation** (one random shift or one
random conductivity draw) as one observation, summarised by its **median
over cord positions**. So n = 8 per bundle and 24 pooled. Two consequences:

- Bundle-level tests have little power by design. A non-significant result
  at n = 8 is not evidence of no effect; read the effect size and the
  pooled row.
- The cord positions are not independent observations and are never
  pooled as if they were. Where positions are tested individually (the
  along-cord figures), each position is its own test across realisations,
  FDR-corrected along the cord.

### Tests and effect sizes

| Situation | Test | Effect size |
|---|---|---|
| paired (same realisation, two conditions) | sign-flip permutation on the paired differences (exact for n ≤ 20) | matched-pairs rank-biserial |
| unpaired | label-shuffling permutation on the mean difference (exact for 8 vs 8) | Cliff's δ |
| metric vs perturbation size | Spearman ρ, permutation p | ρ |

Both effect sizes run from −1 to +1 and are **positive when the first-named
condition shows the larger value** (for r², larger means *more* similar, so
the sign reads the other way from RE, RDM and |gain|). Rough guides for
|Cliff's δ|: < 0.15 negligible, < 0.33 small, < 0.47 medium, otherwise large.
A p-value says whether a difference is distinguishable from none; the effect
size says how consistently one condition exceeds the other; the medians say
by how much in % RE. Report all three.

### What is paired with what

| Comparison | Paired? | Why |
|---|---|---|
| source vs sensor shift, one modality | yes (`source_sensor_paired`) | shift *k* of both uses the same shift vector |
| anything vs conductivity | no | a different kind of draw |
| bundle vs bundle | no | different random shifts |
| Biot-Savart vs BEM | yes | one geometry, solved twice |
| MSG vs ESG | per type, `msg_esg_paired` | only if the same geometric change was applied to both |
| SQUID MSG vs OP-MSG | yes | same lead fields, different noise floor |

### Multiple comparisons

p-values are Benjamini–Hochberg adjusted (`p_fdr`) within one **family ×
metric**, across methods, orientations, axis sets and subsets. So "every
MSG-vs-ESG test on RE" is one family. Along-cord tests are adjusted across
positions within each panel. Stars in the figures are always FDR-adjusted.

### Confidence intervals

95% percentile bootstrap intervals of the median. What they resample is
stated in every table: **realisations** in stages 2–5 ("how much would the
median move for a different set of random shifts of this size"), **source
positions** for the stage-1 model-type comparison. Neither is a
between-subject interval. This is one anatomy, perturbed.

---

## 3. Stage 1 — baseline

**Topoplots.** Unperturbed fields for Biot-Savart MSG, BEM MSG and BEM ESG,
one row per model, with a colour scale per row (the units differ). One
figure per axis slot: X/tangential, Y, Z/radial. ESG has no Y panel, and that
is shown explicitly rather than dropped. Look for how spatially concentrated
each map is: a tight, structured map has more to lose from both a geometric
error and from noise than a broad, smooth one.

**Amplitude along the cord.** RMS and peak field per source, and focality
(fraction of sensors reaching half the peak; smaller = more concentrated).
Amplitude sets the SNR in stage 5, and focality usually predicts which
modality is more sensitive to sensor shifts. Read these before reading the
perturbation results.

**Model-type decomposition.** Biot-Savart against BEM on the unperturbed
geometry. This is the scale of "how different are two reasonable forward
models". Every perturbation effect in stages 2–4 can be judged against it.
It also appears as a row in the stage-6 ranking. A difference that is mostly
gain (the volume conductor changes how much field escapes) is physically
expected; a large RDM means the two models disagree about the field pattern.

---

## 4. Stages 2 and 3 — within one modality

Read in this order.

1. **`within_descriptive`** — how big each perturbation is. The worst-case
   columns (maximum RE, minimum r² over every realisation and position) show
   whether a small median hides a bad region; the `% cord r²<0.95` column
   says how much of the cord is affected.
2. **`decomposition_<type>_<method>_<axes>`** — where along the cord the
   error sits, bundle by bundle, and whether it is gain or topography.
   Expect RE to peak where the field is weakest (RE is relative, so a small
   reference field inflates it): check against the stage-1 amplitude profile
   before calling a peak anatomical.
3. **Dose response** — larger shifts should give larger errors (ρ > 0). A
   flat dose response means the metric is saturated or that direction
   matters more than size. Check the per-realisation scatter: two clouds at
   the same |shift| usually mean the error depends on *which way* the shift
   went.
4. **Types compared** — which perturbation matters most at matched size.
   Source and sensor shifts of the same vector move the cord relative to the
   sensors by the same amount. What differs is whether the anatomy moves
   through the volume conductor (source shift) or the sensors move past a
   fixed volume conductor (sensor shift). So a significant source-vs-sensor
   difference is the effect of the volume conductor moving with the cord.
5. **Gain vs RDM** — each point is a realisation. Points hugging the
   horizontal axis are pure amplitude changes; points rising away from it
   moved the pattern.
6. **Methods compared** (MSG) — whether Biot-Savart and BEM are equally
   sensitive to the same perturbation. If Biot-Savart shows a similar
   perturbation effect to BEM, the sensitivity is geometric (source–sensor
   distance). If they differ, the volume conductor changes how perturbations
   propagate.
7. **Source minus sensor along the cord** — the paired difference at each
   position, with FDR-significant positions marked. Use it to say *where* the
   two perturbation types diverge.

---

## 5. Stage 4 — MSG versus ESG

The core comparison. Both modalities use the same forward model type
(`modality_compare_method`) and see the same perturbations. All metrics are
dimensionless, so comparing fT/nAm with µV/nAm here is legitimate. What is
compared is how much each modality's field changes relative to itself.

**Axis sets.** The whole array (all axes stacked) is the headline: it needs
no correspondence between MSG and ESG axes. Matched-axis rows (ESG
tangential with MSG X, ESG radial with MSG Z) add detail but compare
physically different channels, so read them as secondary.

- **`headline_<axes>`** — every metric × perturbation type, MSG and ESG side
  by side per bundle, with FDR stars.
- **`effect_size_map`** — the same tests as a map; blue = ESG larger,
  red = MSG larger. The orientation version shows whether a modality
  difference holds for all dipole orientations or only some.
- **`modality_interaction`** — whether the source-minus-sensor difference
  itself differs between modalities. A significant row means the two
  modalities *rank the perturbation types differently*: one is more exposed
  to anatomical error, the other to registration error. That is a stronger
  and more useful statement than "MSG is more sensitive".
- **`modality_slopes`** — error per mm of shift (or per % conductivity) in
  each modality, with a CI of the difference. A CI excluding zero means one
  modality degrades faster as registration gets worse.
- **`along_cord_<type>`** — where along the cord the modalities differ.

When a MSG-vs-ESG test is unpaired (sensor shifts by default), it compares
two different sets of random shifts. It is valid, but has less power than a
paired test, and a difference could partly reflect the different shift
directions. The ESG sensor shift has no Z component, so it is also a smaller
shift at the same bundle.

---

## 6. Stage 5 — through sensor noise

### What is simulated

An evoked response is generated through each lead field, trial-averaged
sensor noise is added at the system's own floor and bandwidth, and the field
pattern is recovered from the data with a matched filter
(`ĝ = Y w / wᵀw`). The recovered pattern is scored against the **true,
unperturbed, noise-free** field with the same four metrics. At noise level 0
this reproduces stages 2–4 exactly (`test_pt_pipeline` checks this).

Noise levels are multiples of each system's own floor, so "1×" means
realistic for every system even though the absolute floors (fT/√Hz, µV/√Hz)
cannot be compared.

### Reading the outputs

- **`noise_only_<system>`** — the error noise causes with a perfect model.
  This is the floor any perturbation must rise above. Noise raises RE and
  RDM and lowers r², and it biases **gain upwards**: noise adds energy, so
  `‖ĝ‖ > ‖g‖` on average. A positive gain at high noise is therefore not an
  amplitude effect of the geometry.
- **`noise_pert_<system>_<type>`** — perturbed + noise (colour) against
  unperturbed + noise (grey). Where the two converge, the noise hides the
  perturbation.
- **`detectability_<system>`** — the fraction of realisations whose error
  exceeds the 95th percentile of the noise-only error. **`noise_critical_levels`**
  gives the noise level (× floor) at which that fraction drops below
  `noise_detect_rate`. A critical level above 1× means a perturbation of that
  size would be visible in a real recording with that system; below 1× means
  it would be hidden by the noise.
- **`error_budget`** — the perturbation's share of the squared RE. Because
  the noise is zero-mean and independent of the field,
  `E‖ĝ − a‖² = ‖g − a‖² + E‖n‖²` exactly, so squared RE from the two sources
  adds. 100% means the noise is negligible; near 0% means the error is
  almost all noise.
- **`noise_system_stats` / `systems_<a>_vs_<b>`** — the stage-4 comparisons
  repeated at every noise level. The 0× row is the noise-free result. The
  question is whether a noise-free MSG-vs-ESG difference **survives**, is
  **masked** (the effect size shrinks towards 0 as noise grows) or **reverses**.
- **`noise_within_stats`** — the stage 2/3 magnitude and type tests at 0× and
  at the reference level. A type difference that disappears at 1× is real in
  the model but not observable with that system.

### Caveats to state in any write-up

- White noise only. Real recordings add environmental interference, and for
  ESG cardiac artefact, which is far larger than the amplifier floor used
  here. **ESG results are an optimistic bound.**
- The matched filter assumes the source waveform is known exactly. Real
  analyses estimate it, so the simulated recovery is a best case for every
  system equally.
- One source is active at a time. Crosstalk between sources is not
  modelled.

---

## 7. Stage 6 — every factor on one scale

`summary_<modality>` ranks the forward-model-type difference, every
perturbation bundle, noise alone and perturbation + noise at the reference
level by median RE. It answers "which source of error matters most" for each
modality. The `CI over` column states what each interval resamples, and it
differs between rows. Compare magnitudes across rows, not CI overlap.
`summary_forest` shows the same ranking graphically on a log scale.

---

## 8. When something looks wrong

| Symptom | Likely cause | Check |
|---|---|---|
| RE ≈ 100%, r² ≈ 1, RDM ≈ 0 | units mismatch between a perturbed file and its reference | `leadfield_scale_log.csv`; the stage-0 magnitude warning |
| RE in the 1e9 % range for conductivity only | conductivity files scaled with the standard factor | `unit_scale_mode = 'auto'` |
| every perturbation gives RE = 0 | the perturbed file is identical to the reference (msg_fwd was run on the wrong geometry) | file dates and contents |
| MSG results change after re-loading | front and back arrays were mixed | `analysis_array`; `analysis_array_loaded` in `leadfields_organised.mat` |
| stage-5 curves flat at r² = 1 for every level | signal far above noise, typically a scale error | the SNR column of `noise_only` |
| dose-response ρ ≈ 0 with large RE | error dominated by shift direction, not size | the gain-vs-RDM and dose-response scatters |

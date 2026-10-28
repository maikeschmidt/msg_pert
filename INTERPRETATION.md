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
comparisons of amplitude therefore use |gain|.

**Orientations and axis sets.** Every metric is computed per dipole
orientation (VD, RC, LR) and for `ALL`, the three stacked into one vector,
which is the single number to quote. It is also computed per sensor axis and
for the whole array, all axes stacked. The headline figures and LaTeX tables
use `headline_orientation` (default `ALL`) on the whole array; the CSVs hold
everything.

### Which field is the reference

Wherever two **forward models** are compared (stage 1), the reference (the RE
denominator) is the unperturbed **MRI-derived BEM**, as in msg_fwd and
msg_models. Wherever a **perturbation** is measured, each model is compared
with its own unperturbed field, so the number is that model's sensitivity.
**Biot-Savart** has no volume conductor: it depends only on where the
sources and sensors are. That makes it the anatomy-independent comparison in
three places: its response to a shift is the purely geometric part of the
sensitivity (stage 2, section H), it is unaffected by conductivity by
construction, and in stage 5 it gives the upper bound on what an MSG system
could recover from the same sensors.

**FEM** (MSG only) is the second volume-conductor model, compared against
the BEM in stage 1 and solved for every perturbation, so each result can be
read for both solvers. Conductivity draws are identical for the two, so
"does the FEM respond to conductivity error like the BEM" is a paired
comparison.

### FEM re-meshing

TetGen places its region seeds at random, so two meshes of one geometry give
different FEM lead fields; msg_models' seed-pair study found this can be as
large as the BEM–FEM difference itself. Which FEM results contain it:

| FEM result | Mesh | Contains re-meshing? |
|---|---|---|
| sensor shift | the unshifted geometry's | no |
| conductivity draw | the unshifted geometry's | no |
| source shift | its own (cord and bone moved), `fem_mesh_seed` | **yes** |

So the unshifted geometry is also meshed with `fem_remesh_seeds` (8 seeds by
default) and each re-mesh is scored against the production mesh exactly as
a perturbation is. That spread is the **re-meshing floor**
(`within_fem_remesh`, `fem_remesh_floor`, a row in stage 6). A FEM
source-shift bundle inside the floor cannot be told apart from re-meshing;
only the part above it is the shift. The same caveat applies to every FEM
quantity built from source shifts: the source-vs-sensor difference, the
FEM anatomy contribution, and the FEM source rows of stage 5. The BEM has
no such floor: its boundaries are the surfaces themselves.

---

## 2. The statistics

### Intervals, not tests

Every value in this pipeline is computed exactly from a simulated geometry.
There is no measurement noise in stages 1–4, so there is no null hypothesis
for a significance test to test: a difference between two perturbations is
real by construction. What *is* uncertain is how well the random draws
represent the perturbations that could occur, and that is a question about
spread, which an interval answers directly. This is the convention msg_fwd
uses for its replicate geometries (`st_warp_summary`) and msg_models for its
sweeps, so the three toolboxes report the same kind of number:

| Quantity | How it is reported |
|---|---|
| one condition | median, 95% bootstrap interval of the median, IQR, 95th percentile with its own interval (`pt_describe`) |
| two conditions, paired | median of the per-realisation difference, 95% bootstrap interval (`pt_compare_two`) |
| two conditions, unpaired | difference of medians, 95% bootstrap interval resampling each condition |
| metric vs perturbation size | Spearman ρ with a 95% bootstrap interval over realisations (`pt_rank_corr`) |
| along the cord | the same, at each position |

A comparison is marked **separated** (column `separated`, ◇ in the figures)
when the interval of the difference excludes zero. That describes the
interval; it is not a significance test, and no multiple-comparison
correction is applied because nothing is being tested. Percentiles
everywhere come from msg_fwd's `pctl`, and median intervals from
`st_boot_ci_median`, so a msg_pert interval and a msg_fwd interval are
computed identically.

**The coverage statement**, as in msg_fwd, is worded descriptively:

> In 95% of the 8 realisations, the cord-median RE was at or below X%
> (95% CI …).

Not "a new perturbation has a 95% chance of falling below X". With 8
realisations per bundle the 95th percentile rests on the largest one or two
values; quote its interval with it.

### The unit of observation

Every comparison treats **one perturbation realisation** (one random shift
or one random conductivity draw) as one observation, summarised by its
**median over cord positions** (or, in stage 5, by its value at the C7
source). So n = 8 per bundle and 24 pooled. Two consequences:

- Bundle-level intervals are wide by design. An interval that includes zero
  at n = 8 is not evidence of no difference; read the pooled row and the
  medians.
- The cord positions are not independent observations and are never
  pooled as if they were. Where positions are compared individually (the
  along-cord figures), each position has its own interval across
  realisations.

Signs: every difference is **first-named minus second-named**, so positive
means the first condition shows the larger value. For r², larger means
*more* similar, so the sign reads the other way from RE, RDM and |gain|. The
stage-4 difference map converts r² to 1 − r² so that + always means "more
affected".

### What is paired with what

| Comparison | Paired? | Why |
|---|---|---|
| source vs sensor shift, one modality | yes (`source_sensor_paired`) | shift *k* of both uses the same shift vector |
| anything vs conductivity | no | a different kind of draw |
| bundle vs bundle | no | different random shifts |
| Biot-Savart vs BEM | yes | one geometry, solved twice |
| MSG vs ESG | per type, `msg_esg_paired` | only if the same geometric change was applied to both |
| SQUID MSG vs OP-MSG | yes | same lead fields, different noise floor |
| MSG system vs its Biot-Savart bound | yes | same shifts, same noise draws |

### What the intervals resample

Stated in every table: **realisations** in stages 2–5 ("how much would the
median move for a different set of random shifts of this size"), **source
positions** for the stage-1 model-type comparison. Neither is a
between-subject interval. This is one anatomy, perturbed.

---

## 3. Stage 1 — baseline

**Topoplots.** Unperturbed fields for Biot-Savart, BEM and FEM MSG and BEM ESG,
one row per model, with a colour scale per row (the units differ). One
figure per MSG axis: X, Y, Z. The ESG potential has no direction, so the same
ESG map appears beside every MSG axis. Look for how spatially concentrated
each map is: a tight, structured map has more to lose from both a geometric
error and from noise than a broad, smooth one.

**Amplitude along the cord.** RMS and peak field per source, and focality
(fraction of sensors reaching half the peak; smaller = more concentrated).
Amplitude sets the SNR in stage 5, and focality usually predicts which
modality is more sensitive to sensor shifts. Read these before reading the
perturbation results.

**Model-type decomposition.** Biot-Savart and FEM, each against the
MRI-derived BEM on the unperturbed geometry. FEM vs BEM is the scale of "how
different are two reasonable volume-conductor models" (msg_fwd's central
comparison, here on this geometry); Biot-Savart vs BEM is what the volume
conductor does at all. Every perturbation effect in stages 2–4 can be judged
against them, and both appear as rows in the stage-6 ranking. A difference that is mostly
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
3. **Dose response** — larger shifts should give larger errors (ρ > 0, its
   interval above zero). A flat dose response means the metric is saturated
   or that direction matters more than size. Check the per-realisation scatter: two clouds at
   the same |shift| usually mean the error depends on *which way* the shift
   went.
4. **Types compared** — which perturbation matters most at matched size.
   Source and sensor shifts of the same vector move the cord relative to the
   sensors by the same amount. What differs is whether the anatomy moves
   through the volume conductor (source shift) or the sensors move past a
   fixed volume conductor (sensor shift). So a separated source-vs-sensor
   difference is the effect of the volume conductor moving with the cord.
5. **Gain vs RDM** — each point is a realisation. Points hugging the
   horizontal axis are pure amplitude changes; points rising away from it
   moved the pattern.
6. **Methods compared** (MSG) — whether Biot-Savart, BEM and FEM are
   equally sensitive to the same perturbation, every pair as a paired
   difference per realisation. BEM vs FEM on sensor shifts and conductivity
   is a clean solver comparison; on source shifts it includes the FEM's
   re-meshing (item 9).
7. **Source minus sensor along the cord** — the paired difference at each
   position, marked where its interval excludes zero. Use it to say *where*
   the two perturbation types diverge.
8. **Anatomy contribution** (MSG; `within_anatomy`,
   `anatomy_contribution_<bem|fem>`, `anatomy_along_cord_<bem|fem>_<type>`) —
   Biot-Savart's response to a shift is the purely geometric sensitivity,
   because it has no anatomy. BEM (or FEM) minus Biot-Savart for the same
   shift is therefore what the volume conductor adds (positive) or removes
   (negative); the model/BS column gives it as a ratio of the error, where 1
   means the anatomy adds nothing. BEM and FEM agreeing here means the
   anatomy contribution is a property of the anatomy, not of the solver. Near zero:
   the perturbation is a source–sensor distance effect and an anatomy-free
   model predicts it. Clearly positive: the volume conductor amplifies the
   error, so the anatomy has to be known to predict it. Conductivity has no
   row, since Biot-Savart does not depend on it; there the whole effect is
   anatomical.
9. **FEM re-meshing floor** (MSG; `within_fem_remesh`, `fem_remesh_floor`) —
   see section 1. Read every FEM source-shift number against it.

---

## 5. Stage 4 — MSG versus ESG

The core comparison. Both modalities use the same forward model type
(`modality_compare_method`) and see the same perturbations. All metrics are
dimensionless, so comparing fT/nAm with µV/nAm here is legitimate. What is
compared is how much each modality's field changes relative to itself.

**Axis sets.** MSG and ESG are compared on the whole array: every MSG axis
stacked, against every electrode. The ESG is one set of electrodes measuring
a potential, which has no direction, so there is no ESG channel to match to a
single MSG axis.

**The ESG electrodes.** msg_coreg builds the ESG array from the grid points of
the magnetometer array, which a triaxial array holds three times; the stored
electrode set therefore holds every electrode twice. The copies measure the
identical potential, so msg_pert keeps one copy on load
(`pt_dedup_electrodes`, after checking they are identical). The potential is
common-average referenced over the array (msg_coreg's reference matrix,
applied by FieldTrip): r² does not depend on that choice, RE, RDM and gain
do, and a different reference (a single electrode, as often used in ESG
recordings) would change them.

- **`headline_<axes>`** — every metric × perturbation type, MSG and ESG side
  by side per bundle; ◇ where the interval of MSG − ESG excludes zero.
- **`difference_map`** — the same comparisons as a map of the relative
  difference `(MSG − ESG) / mean`, in error form (1 − r² for r²), so
  red = MSG more affected and blue = ESG more affected on every metric. The
  orientation version shows whether a modality difference holds for all
  dipole orientations or only some.
- **`modality_interaction`** — whether the source-minus-sensor difference
  itself differs between modalities. A separated row means the two
  modalities *rank the perturbation types differently*: one is more exposed
  to anatomical error, the other to registration error. That is a stronger
  and more useful statement than "MSG is more sensitive".
- **`modality_slopes`** — error per mm of shift (or per % conductivity) in
  each modality, with a CI of the difference. A CI excluding zero means one
  modality degrades faster as registration gets worse.
- **`along_cord_<type>`** — where along the cord the modalities differ.

When a MSG-vs-ESG comparison is unpaired (sensor shifts by default), it
compares two different sets of random shifts. It is valid, but its interval
is wider than a paired one, and a difference could partly reflect the
different shift directions. The ESG sensor shift has no Z component, so it is also a smaller
shift at the same bundle.

---

## 6. Stage 5 — through sensor noise

### What is simulated

The evoked response is the spinal **N13** after median nerve stimulation: a
stationary generator in the dorsal horn, simulated as a single source at C7
(`evoked_src_mm`) with a response peaking at 13 ms (`noise_latency`,
`noise_env_sd`, `noise_freq`; inside the 8–16 ms spinal window used by
msg_analysis). It is generated through each lead field, trial-averaged
sensor noise is added at the system's own floor and bandwidth, and the field
pattern is recovered from the data with a matched filter (`ĝ = Y w / wᵀw`).
The recovered pattern is scored with the same four metrics against the
**noise-free, unperturbed field of the same forward model**. At noise level 0
this reproduces stages 2–4 exactly (`test_pt_pipeline` checks this, at C7
and for the cord median).

Seven systems are simulated:

| System | Forward model | Role |
|---|---|---|
| SQUID MSG, OP-MSG | MRI-derived BEM | realistic MSG |
| ESG | MRI-derived BEM | realistic ESG |
| SQUID MSG (FEM), OP-MSG (FEM) | FEM | realistic MSG through the second solver |
| SQUID MSG (Biot-Savart), OP-MSG (Biot-Savart) | Biot-Savart | anatomy-free upper bound |

The FEM and Biot-Savart systems use the same sensors and the same noise
draws as their BEM counterparts. A conclusion about recoverability that
holds for both the BEM and FEM systems does not depend on the solver
(`systems_squid_msg_fem_vs_squid_msg`, `systems_op_msg_fem_vs_op_msg`). Without a volume conductor nothing attenuates or
redistributes the field, so they show the best an MSG system with those
sensors could do. The gap between the two (`systems_squid_msg_bs_vs_squid_msg`,
`systems_op_msg_bs_vs_op_msg`) is what the anatomy costs in recoverability.

Noise levels are multiples of each system's own floor, so "1×" means
realistic for every system even though the absolute floors (fT/√Hz, µV/√Hz)
cannot be compared.

### Two summaries

Every stage-5 output exists twice:

- **`5_noise/c7/`** — the evoked source alone. This is the headline: the
  signal a median nerve experiment records. The baseline noise question
  (`noise_only`: how well is the N13 recovered with perfect coregistration,
  as noise grows from 0× to 8× each system's floor) and the recoverability
  question (`noise_pert`, `detectability`, `noise_critical_levels`,
  `error_budget`: does the N13 survive once the geometry is wrong) are both
  answered here.
- **`5_noise/cord/`** — the same response placed at every cord position,
  summarised by the cord median as in stages 2–4. Use it to say whether C7
  is typical of the cord or a favourable or unfavourable level.

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
  repeated at every noise level, as a difference with its interval. The 0×
  row is the noise-free result. The question is whether a noise-free
  MSG-vs-ESG difference **survives**, is **masked** (the difference shrinks
  towards 0 and its interval comes to include it as noise grows) or
  **reverses**. The Biot-Savart-vs-BEM pairs show how far the realistic
  model sits below the anatomy-free bound at each level.
- **`noise_within_stats`** — the stage 2/3 magnitude and type comparisons at
  0× and at the reference level. A type difference that disappears at 1× is
  real in the model but not observable with that system.

### Caveats to state in any write-up

- White noise only. Real recordings add environmental interference, and for
  ESG cardiac artefact, which is far larger than the amplifier floor used
  here. **ESG results are an optimistic bound.**
- **SQUID-MSG and OP-MSG share one set of sensor positions**, the msg_coreg
  array about 10 mm from the skin, where OPMs sit. Real SQUID sensors sit
  several centimetres further out, inside the dewar, and the field falls
  steeply with distance. The two systems therefore differ here only in noise
  floor and bandwidth, so **SQUID results are an optimistic bound relative
  to OP-MSG.**
- The matched filter assumes the source waveform is known exactly. Real
  analyses estimate it, so the simulated recovery is a best case for every
  system equally.
- One source is active at a time. Crosstalk between sources, and the
  cortical N20 that follows the N13, are not modelled.
- `evoked_src_mm` places C7 by distance along the source space. Check it
  against the vertebral levels of the geometry before quoting C7 results.

---

## 7. Stage 6 — every factor on one scale

`summary_<modality>` ranks the forward-model-type differences, every
perturbation bundle for every model, the FEM re-meshing floor, noise alone
and perturbation + noise at the reference level by median RE, every row as a
cord median. It answers "which source of
error matters most" for each modality. `summary_<modality>_c7` is the same
ranking at the evoked C7 source, from stage 5: each perturbation noise-free,
noise alone, and perturbation + noise. The `CI over` column states what each
interval resamples, and it differs between rows. Compare magnitudes across
rows, not CI overlap. `summary_forest` and `summary_forest_c7` show the
rankings on a log scale.

---

## 8. Stage 7 — cord activity: N11 volley and N13 extent

### What is simulated

Stages 5 and 6 test one source at a time. Stage 7 simulates the spinal
response as a recording, every active source summed over time:

- **N11** — the afferent volley entering the cord at C7 at 10.5 ms and
  ascending in the dorsal columns at 60 m/s. A propagating volley is
  modelled as a quadrupole: a leading and a trailing longitudinal dipole of
  opposite sign, 30 mm apart (conduction velocity × action-potential
  duration).
- **N13** — the stationary dorsal-horn response at C7, transverse, peaking
  at 13 ms. Its spatial extent is swept from one source to a Gaussian patch
  40 mm wide (full width at half maximum) along the cord, with the same
  total moment at every extent.

`source_model.png` shows both on a cord × time map. `recordings_<modality>`
shows each recording the way an evoked response is read in real data: the
RMS over all channels against time (a global field power), noise-free (black)
and in three noisy realisations at 1× the floor, against the RMS of the noise
alone (dashed), one row per sensor technology and one column per forward
model, with the N11 and N13 windows shaded. It is the figure to set beside an
empirical average with the same number of trials. Every system is analysed in
at most `cord_max_bandwidth_hz` (1 kHz), as a real analysis band would be;
OP-MSG keeps its own 150 Hz. Every number is computed in two windows: the N11 window (the
volley's travel) and the N13 window (its peak ± 2 s.d.). Each system's
bandwidth is applied, so a narrow-band system smears the fast N11 as it
would in a real recording.

### What it can and cannot say

Every recording here is generated by the model being tested, so the
simulation cannot say whether BEM or FEM is *right* for spread sources. It
can say three things, one per output:

| Output | Question | How to read it |
|---|---|---|
| `spread_effect`, `cord_spread_effect` | how much does spreading change the recording at all, within one model? | Read this first. It bounds every other extent effect: if a 40 mm N13 changes the recording by only a few percent against a focal one, no measure below can change by much more. That is the physics of a deep source blurred along the cord by much less than its depth. |
| `divergence`, `cord_divergence` | do BEM and FEM disagree more as the N13 spreads? | RE/RDM rising with extent: the solvers' difference grows for spread activity, so the choice of model matters more there. Flat: it does not depend on extent. |
| `snr_by_system` | which system, through which model, gives the strongest N11 and N13? | SNR at 1× the floor at the headline extent, bars grouped by sensor technology, one bar (and shape) per model. Exact, no noise draws. |
| `detection_vs_noise` | how far is each system from detecting each component? | Fraction of noise realisations detected against noise level. Where a curve crosses the 80% line is its critical level; where it sits at 1× is what a real recording with these trials would show. |
| `detectability_snr`, `cord_snr` | is a spread N13 more or less visible through each system's noise, under each model? | SNR against extent, as % of the focal N13 (top) and relative to the BEM system with the same sensors (bottom, 1 = identical). Exact, so changes of a few percent are visible. |
| `detectability_<window>`, `cord_critical_levels`, `cord_bem_vs_fem_noise` | up to what noise level is each component still detected? | The coarser summary. Built from noise draws, so it moves in steps: a change of a few percent in SNR usually leaves it unchanged. Equal critical levels across extents with a slightly falling SNR mean the extent effect is real but small. Above 1× = visible in a real recording. |
| `perturbation_<modality>`, `cord_perturbation_bem_vs_fem` | does spreading make either model more or less robust to the source, sensor and conductivity errors? | FEM minus BEM per extent, paired by realisation; FEM source shifts against the grey re-meshing band. |

Whether the FEM is the better model for wider sources is then tested on
real data in msg_analysis: fit the measured N13 with BEM and with FEM lead
fields, focal and spread, and compare how much of the data each explains.
These outputs say where in that comparison a difference should be expected.

### Detection

A component counts as detected in one noise realisation when the mean
square of the measured recording in its window exceeds the 95th percentile
(`noise_detect_pct`) of the same statistic for noise alone. The critical
noise level is the highest level at which at least `noise_detect_rate` of
realisations are detected, by the same rule as stage 5
(`pt_critical_level`). Noise draws are shared across extents and across
BEM, FEM and Biot-Savart for the same sensors, so differences between
models are not differences between noise draws.

---

## 9. When something looks wrong

| Symptom | Likely cause | Check |
|---|---|---|
| RE ≈ 100%, r² ≈ 1, RDM ≈ 0 | units mismatch between a perturbed file and its reference | `leadfield_scale_log.csv`; the stage-0 magnitude warning |
| RE in the 1e9 % range for conductivity only | conductivity files scaled with the standard factor | `unit_scale_mode = 'auto'` |
| every perturbation gives RE = 0 | the perturbed file is identical to the reference (msg_fwd was run on the wrong geometry) | file dates and contents |
| MSG results change after re-loading | front and back arrays were mixed | `analysis_array`; `analysis_array_loaded` in `leadfields_organised.mat` |
| stage-5 curves flat at r² = 1 for every level | signal far above noise, typically a scale error | the SNR column of `noise_only` |
| dose-response ρ ≈ 0 with large RE | error dominated by shift direction, not size | the gain-vs-RDM and dose-response scatters |

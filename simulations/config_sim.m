% config_sim - Shared configuration for the msg_pert noise-simulation module
%
% Defines the models, sensor systems, noise floors, and signal parameters
% used by the sim_* scripts. Run as a script at the top of each one.
%
% SCOPE
%   The perturbation pipeline (pt_*) asks how the forward field changes when
%   the *model* is wrong (source position, sensor position, conductivity).
%   This module asks a different question: given a PERFECT forward model, how
%   much does *sensor noise* degrade the measurable field, and does that
%   degradation differ between SQUID-MSG, OP-MSG, and ESG?
%
%   Motivating hypothesis: BEM MSG produces highly individualised, spatially
%   sharp field maps, whereas Biot-Savart MSG and BEM ESG produce smoother,
%   more diffuse maps. Smoother maps should be more robust to both model error
%   and sensor noise. sim_plot_topoplots visualises that claim directly.
%
% VARIABLES DEFINED:
%   Paths:
%     msg_geoms_path     - Folder holding the MSG geometry .mat files
%     esg_geoms_path     - Folder holding the ESG geometry .mat files
%                          (MSG and ESG geometries live in separate folders;
%                           set both to the same path if yours do not)
%     sim_save_dir       - Output folder for sim figures and tables
%     sim_out_dir        - Output folder for sim .mat files
%
%   Models (sim_models struct array) — one entry per forward model to compare:
%     .label      - Display name, e.g. 'MSG — BEM'
%     .short      - Filename-safe stem, e.g. 'msg_bem'
%     .front      - Full path to the FRONT array leadfield .mat
%     .back       - Full path to the BACK array leadfield .mat
%     .var        - Variable name inside the .mat
%     .scale      - Unit scale factor applied on load
%     .is_meg     - true = magnetic (fT/nAm), false = electric (uV/nAm)
%     .geom_file  - Full path to the geometry .mat holding sensor positions
%     .axis_names - Display name of each sensor axis, in leadfield channel
%                   order. MSG is triaxial (X/Y/Z field components); ESG has
%                   two electrode sets (tangential and radial).
%     .axis_slot  - Which COMPARISON SLOT each of this model's axes belongs in.
%                   The topoplot figures are organised one-per-slot, so this is
%                   what decides which ESG panel sits underneath which MSG panel.
%                   Slots follow the MSG triaxial convention:
%                     slot 1 = X   slot 2 = Y   slot 3 = Z
%                   ESG measures potentials, not field components, but its two
%                   electrode sets correspond to MSG axes:
%                     tangential -> slot 1 (X),  radial -> slot 3 (Z)
%                   so ESG is axis_slot = [1 3] — it has no slot-2 panel, and
%                   the slot-2 figure shows a placeholder for the ESG row.
%
%   Sensor systems (sim_systems struct array) — one entry per real system:
%     .label          - Display name, e.g. 'SQUID MSG'
%     .short          - Filename-safe stem
%     .model          - Index into sim_models: which forward model it measures
%     .noise_baseline - White noise floor, in units/sqrt(Hz)
%     .noise_unit     - Unit string for FIGURE text (TeX: renders \surd as a
%                       radical). Do NOT fprintf this — the console prints the
%                       raw escape, e.g. "fT/\surdHz".
%     .noise_unit_txt - Plain-ASCII unit string for console output.
%     .color          - [1x3] RGB
%
%   Signal:
%     sim_fs              - Sampling rate (Hz)
%     sim_duration        - Epoch length (s)
%     sim_freq            - Source oscillation frequency (Hz)
%     sim_dipole_nAm      - Source strength (nA*m)
%
%   Noise sweep:
%     sim_noise_factors   - Multipliers applied to each system's baseline.
%                           The SAME factors are used for every system, so the
%                           x-axis is comparable even though the absolute noise
%                           floors are in different physical units.
%     sim_n_realisations  - Noise realisations averaged per (source, level)
%     sim_noise_seed      - RNG seed for reproducibility
%
% UNITS WARNING:
%   MSG leadfields are fT/nAm and ESG leadfields are uV/nAm. Absolute noise
%   levels are therefore NOT comparable between MSG and ESG. Everything in
%   this module is compared in units of "x baseline noise floor", which IS
%   comparable — it asks how each system degrades relative to its own
%   real-world operating point.
%
% NOISE FLOOR REFERENCES:
%   SQUID MSG  2-5   fT/sqrt(Hz)   Brookes et al. 2022; MSG systems <4-5
%   OP-MSG     7-20  fT/sqrt(Hz)   Brookes et al. 2022; O'Neill et al. 2024
%   ESG        ~1    uV/sqrt(Hz)   Not formally specified; EEG-grade amplifier
%                                  noise assuming good reference + impedance.
%                                  In practice cardiac artefact (~29 uV
%                                  cervical, ~657 uV lumbar) dominates and is
%                                  far larger than this — so the ESG curve here
%                                  is an OPTIMISTIC bound, not a realistic one.
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_pert
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk


% =========================================================================
% PACKAGE PATH BOOTSTRAP (self-contained)
% =========================================================================
% This module lives in msg_pert/simulations and is self-contained. Add its own
% folder + functions/, and the msg_pert root (for pt_add_functions, which pulls
% in msg_fwd for plot_topoplot_publication).
sim_dir = fileparts(mfilename('fullpath'));
if isempty(sim_dir); sim_dir = fileparts(which('config_sim')); end
addpath(sim_dir);
addpath(fullfile(sim_dir, 'functions'));
addpath(fileparts(sim_dir));   % msg_pert root (pt_add_functions -> msg_fwd)


% =========================================================================
% USER CONFIGURATION — paths
% =========================================================================

% MSG and ESG geometries live in separate folders. If yours happen to share a
% folder, just set both to the same path.
msg_geoms_path = 'D:\Simulations\Pertubations\geometries';        % SET THIS
esg_geoms_path = 'D:\Simulations\Pertubations\geoms_elec';    % SET THIS

sim_save_dir   = 'D:\Simulations\Pertubations\results\simulation';% SET THIS
sim_out_dir    = 'D:\Simulations\Pertubations\fields\simulation'; % SET THIS

% Roots of the leadfield sets (perfect + perturbed forward fields).
% These hold BOTH the original and the source/sensor-shift geometries, each in
% its own geometries_<short> subfolder (standard msg_fwd output layout).
msg_bem_root   = 'D:\Simulations\Pertubations\fields\mag\bem';        % SET THIS
msg_bslaw_root = 'D:\Simulations\Pertubations\fields\mag\bs_law';     % SET THIS
esg_bem_root   = 'D:\Simulations\Pertubations\fields\elec\bem_elec';  % SET THIS

% Roots holding the BEM-conductivity leadfields (produced by
% run_conductivity_perturbation) — one per modality. Cond files live under the
% root in geometries_<cond_short>/leadfield_<cond_short>_bem_cond_bundleB_shiftS_<array>.mat
msg_cond_root  = 'D:\Simulations\Pertubations\fields\mag\bem_cond_msg';   % SET THIS
esg_cond_root  = 'D:\Simulations\Pertubations\fields\elec\bem_cond_esg';  % SET THIS
cond_short     = 'original_source_original';   % geometry the conductivity change was applied to

if ~exist(sim_save_dir, 'dir'); mkdir(sim_save_dir); end
if ~exist(sim_out_dir,  'dir'); mkdir(sim_out_dir);  end


% =========================================================================
% FORWARD MODELS TO COMPARE
% =========================================================================
% File naming follows msg_fwd output (see pt_load_leadfields):
%   BEM:    <root>/geometries_<short>/leadfield_<short>_bem_<array>.mat
%           variable leadfield_cord   | scale 1e15 (T/nAm -> fT/nAm)
%   BS law: <root>/leadfield_geometries_<short>_bslaw_<array>.mat
%           variable leadfield_bs     | scale 1 (already fT/nAm)
%
% UNIT SCALES — depend on whether the ft_prepare_leadfield dipoleunit patch
% is active in your FieldTrip.
%
% run_bem_leadfields sets cfg.dipoleunit = 'nA*m', but that option only takes
% effect with a patched ft_prepare_leadfield (see the note at
% run_sphere_leadfields.m:325). UNPATCHED, FieldTrip returns the leadfield per
% AMPERE-metre, not per nanoampere-metre — a factor of 1e9 larger than intended.
%
%   Patched   (per nA*m):  BEM MSG x 1e15 (T->fT),   BEM ESG x 1e6  (V->uV)
%   Unpatched (per A*m):   BEM MSG x 1e6,            BEM ESG x 1e-3
%                          (i.e. the patched value x 1e-9)
%
% Symptom of getting this wrong: BEM peaks come out ~1e9 too large, and BEM
% disagrees wildly with Biot-Savart even though msg_fwd shows the two agreeing
% closely for MSG. Biot-Savart is immune because run_biot_savart_leadfields
% computes its own fT/nAm scale directly and never goes through FieldTrip.
%
% The amplitude diagnostic printed by sim_plot_topoplots checks this for you.

bem_patched = false;   % SET THIS: true if your ft_prepare_leadfield honours cfg.dipoleunit

if bem_patched
    msg_bem_scale = 1e15;   % T/nAm -> fT/nAm
    esg_scale     = 1e6;    % V/nAm -> uV/nAm
else
    msg_bem_scale = 1e6;    % T/(A*m) -> fT/nAm
    esg_scale     = 1e-3;   % V/(A*m) -> uV/nAm
end

% CONDUCTIVITY leadfields need a DIFFERENT scale than the standard BEM ones.
% run_bem_leadfields saves the raw leadfield, but run_conductivity_perturbation
% bakes an extra *1e15 into the file it saves (and labels it 'fT/nAm') for BOTH
% modalities. So a cond file is 1e15x larger than the matching standard file;
% loading it with the standard scale makes the signal ~1e15 too big, the noise
% negligible, and every cond r^2 curve pins flat at 1.
%
% The correction is therefore cond_scale = standard_scale / 1e15. If you
% regenerate the cond leadfields WITHOUT that bake (i.e. matching
% run_bem_leadfields), set cond_extra_scale = 1 so cond and standard agree.
cond_extra_scale = 1e15;   % SET TO 1 if run_conductivity_perturbation no longer bakes *1e15

% Models are now defined by COMPONENTS (root, method, scale, ...) rather than
% fixed file paths, so the same model can be loaded for ANY geometry stem. The
% actual .mat path is built at load time by sim_lf_path(model, short, array),
% and the geometry .mat (sensor positions) by sim_geom_file(model, short).
sim_models = struct('label', {}, 'id', {}, 'method', {}, 'root', {}, ...
                    'cond_root', {}, 'cond_scale', {}, 'geoms_path', {}, 'var', {}, 'scale', {}, ...
                    'is_meg', {}, 'n_axes', {}, 'axis_names', {}, 'axis_slot', {});

% --- Model 1: MSG, Biot-Savart (infinite homogeneous space — smooth fields)
sim_models(1).label      = 'MSG — Biot-Savart';
sim_models(1).id         = 'msg_bslaw';
sim_models(1).method     = 'bslaw';       % flat folder, leadfield_geometries_<short>_bslaw_<array>.mat
sim_models(1).root       = msg_bslaw_root;
sim_models(1).cond_root  = '';            % Biot-Savart has no conductivity variant
sim_models(1).cond_scale = [];           % (unused)
sim_models(1).geoms_path = msg_geoms_path;
sim_models(1).var        = 'leadfield_bs';
sim_models(1).scale      = 1;
sim_models(1).is_meg     = true;
sim_models(1).n_axes     = 3;
sim_models(1).axis_names = {'X-axis', 'Y-axis', 'Z-axis'};
sim_models(1).axis_slot  = [1 2 3];

% --- Model 2: MSG, BEM (individualised anatomy — sharp fields)
sim_models(2).label      = 'MSG — BEM';
sim_models(2).id         = 'msg_bem';
sim_models(2).method     = 'bem';         % <root>/geometries_<short>/leadfield_<short>_bem_<array>.mat
sim_models(2).root       = msg_bem_root;
sim_models(2).cond_root  = msg_cond_root;
sim_models(2).cond_scale = msg_bem_scale / cond_extra_scale;   % cond files carry an extra *1e15
sim_models(2).geoms_path = msg_geoms_path;
sim_models(2).var        = 'leadfield_cord';
sim_models(2).scale      = msg_bem_scale;
sim_models(2).is_meg     = true;
sim_models(2).n_axes     = 3;
sim_models(2).axis_names = {'X-axis', 'Y-axis', 'Z-axis'};
sim_models(2).axis_slot  = [1 2 3];

% --- Model 3: ESG, BEM (surface potentials — smooth fields)
sim_models(3).label      = 'ESG — BEM';
sim_models(3).id         = 'esg_bem';
sim_models(3).method     = 'bem';
sim_models(3).root       = esg_bem_root;
sim_models(3).cond_root  = esg_cond_root;
sim_models(3).cond_scale = esg_scale / cond_extra_scale;       % cond files carry an extra *1e15
sim_models(3).geoms_path = esg_geoms_path;
sim_models(3).var        = 'leadfield_cord';
sim_models(3).scale      = esg_scale;
sim_models(3).is_meg     = false;
% ESG electrodes are not a Cartesian triad: two electrode sets, tangential and
% radial. n_axes MUST be stated — the ESG channel count (e.g. 342 = 2 x 171)
% can be divisible by 3, so any code that guesses the axis count from
% divisibility will wrongly read it as 3 x 114 and mis-slice the leadfield.
sim_models(3).n_axes     = 2;
sim_models(3).axis_names = {'Tangential', 'Radial'};
% ESG tangential lines up with the MSG X-axis, ESG radial with the MSG Z-axis.
sim_models(3).axis_slot  = [1 3];


% =========================================================================
% GEOMETRY VARIANTS TO LOOP OVER
% =========================================================================
% Each entry is one leadfield geometry the noise simulation runs on. The loop
% (sim_run_geometries) creates a per-variant output subfolder, simulates the
% evoked response + noise for every applicable system, and saves the r^2
% curves; sim_plot_comparison then overlays all variants.
%
% Fields:
%   .name   subfolder + legend label (filename-safe)
%   .short  geometry stem (WITHOUT 'geometries_' prefix) — how the leadfield
%           file is named. For source/sensor shifts pick ONE representative
%           shift per bundle.
%   .kind   'standard'  -> leadfield_<short>_bem_<array>.mat  (source/sensor/original)
%           'cond'      -> leadfield_<cond_short>_bem_cond_bundleB_shiftS_<array>.mat
%   .group  free-text grouping for the comparison plot ('baseline'/'source'/
%           'sensor'/'cond'); variants in a group are drawn in one panel.
%   .bundle small/medium/large label used for colour + line ordering.
%   For 'cond' only: .cond_bundle, .cond_shift  and (optional) .root override.
%
% NOTE ON MODALITY: source, sensor, AND conductivity variants all exist for
% both MSG and ESG, each under its own root (model.root for source/sensor,
% model.cond_root for conductivity). Any (variant, system) whose leadfield file
% is missing is skipped automatically, so this stays robust if some combination
% was not generated.

sim_geometries = struct('name', {}, 'short', {}, 'kind', {}, 'group', {}, ...
                        'bundle', {}, 'cond_bundle', {}, 'cond_shift', {}, ...
                        'root', {});

gi = 0;
% -- Baseline (unperturbed) --------------------------------------------------
gi = gi+1;
sim_geometries(gi) = struct('name','original', 'short','original_source_original', ...
    'kind','standard', 'group','baseline', 'bundle','none', ...
    'cond_bundle',[], 'cond_shift',[], 'root','');

% -- Source shifts: one representative shift per bundle ----------------------
src_reps = {1, 1, 1};   % SET THIS: which shift index to use for bundle 1/2/3
for b = 1:3
    gi = gi+1;
    lbl = {'small','medium','large'};
    sim_geometries(gi) = struct( ...
        'name', sprintf('source_%s', lbl{b}), ...
        'short', sprintf('original_source_bundle%d_shift%d', b, src_reps{b}), ...
        'kind','standard', 'group','source', 'bundle',lbl{b}, ...
        'cond_bundle',[], 'cond_shift',[], 'root','');
end

% -- Sensor shifts: one representative shift per bundle ----------------------
sen_reps = {1, 1, 1};   % SET THIS: which shift index to use for bundle 1/2/3
for b = 1:3
    gi = gi+1;
    lbl = {'small','medium','large'};
    sim_geometries(gi) = struct( ...
        'name', sprintf('sensor_%s', lbl{b}), ...
        'short', sprintf('original_sensor_bundle%d_shift%d', b, sen_reps{b}), ...
        'kind','standard', 'group','sensor', 'bundle',lbl{b}, ...
        'cond_bundle',[], 'cond_shift',[], 'root','');
end

% -- Conductivity: one representative shift per bundle (MSG and ESG) ---------
% Cond files all live under the single cond_short geometry, named with
% bem_cond_bundleB_shiftS, under each modality's own cond root (model.cond_root).
cond_reps = {1, 1, 1};   % SET THIS: which shift index to use for bundle 1/2/3
for b = 1:3
    gi = gi+1;
    lbl = {'small','medium','large'};
    sim_geometries(gi) = struct( ...
        'name', sprintf('cond_%s', lbl{b}), ...
        'short', cond_short, ...
        'kind','cond', 'group','cond', 'bundle',lbl{b}, ...
        'cond_bundle', b, 'cond_shift', cond_reps{b}, 'root', '');
end

% Colour per bundle magnitude (shared across groups, light -> dark)
sim_bundle_colors = containers.Map( ...
    {'none','small','medium','large'}, ...
    {[0.2 0.2 0.2], [0.30 0.65 0.90], [0.10 0.40 0.70], [0.03 0.20 0.45]});


% =========================================================================
% SENSOR SYSTEMS (noise floors)
% =========================================================================
% SQUID and OP-MSG measure the SAME field (same array, same BEM leadfield) —
% they differ ONLY in their noise floor. That is the intended comparison:
% it isolates the effect of sensor noise from the effect of geometry.
%
% Set .model to the index of the forward model each system measures.
%
% .bandwidth_hz — the system's measurement bandwidth (upper cutoff, Hz).
%   This is NOT a cosmetic detail. Total noise POWER scales with bandwidth, so
%   the time-domain noise s.d. is
%       sigma = noise_density * sqrt(bandwidth)
%   A system that only measures out to 150 Hz therefore admits sqrt(500/150)
%   = 1.8x less noise than one measuring out to 500 Hz, at the same noise
%   density. Ignoring bandwidth would systematically penalise the narrow-band
%   systems.
%
%   The EFFECTIVE bandwidth is min(.bandwidth_hz, sim_fs/2): a system cannot
%   measure above the Nyquist frequency of the simulation's sampling rate. With
%   sim_fs = 1000 Hz, Nyquist is 500 Hz, so SQUID's 1 kHz spec is clipped to
%   500 Hz here. Raise sim_fs if you want its full bandwidth represented.
%
%   Set .bandwidth_hz = Inf to mean "whatever Nyquist allows".
%
% CAVEAT: the bandwidth is applied to the NOISE only. The 90 Hz evoked burst
% sits inside every system's passband, so a real anti-alias filter would leave
% the signal essentially untouched — but this does mean the simulation does not
% model any signal attenuation near a system's cutoff.

sim_systems = struct('label', {}, 'short', {}, 'model', {}, ...
                     'noise_baseline', {}, 'noise_unit', {}, ...
                     'noise_unit_txt', {}, 'bandwidth_hz', {}, 'color', {});

sim_systems(1).label          = 'SQUID MSG';
sim_systems(1).short          = 'squid_msg';
sim_systems(1).model          = 2;      % MSG BEM
sim_systems(1).noise_baseline = 5;      % fT/sqrt(Hz)  (2-5 typical)
sim_systems(1).noise_unit     = 'fT/\surdHz';
sim_systems(1).noise_unit_txt = 'fT/sqrt(Hz)';
sim_systems(1).bandwidth_hz   = 1000;   % up to 1 kHz (clipped to Nyquist)
sim_systems(1).color          = [0.10, 0.30, 0.80];

sim_systems(2).label          = 'OP-MSG';
sim_systems(2).short          = 'op_msg';
sim_systems(2).model          = 2;      % MSG BEM
sim_systems(2).noise_baseline = 20;     % fT/sqrt(Hz)  (7-20 typical)
sim_systems(2).noise_unit     = 'fT/\surdHz';
sim_systems(2).noise_unit_txt = 'fT/sqrt(Hz)';
sim_systems(2).bandwidth_hz   = 750;    % OPM dynamic range rolls off ~150 Hz (but for comparison we clip to Nyquist)
sim_systems(2).color          = [0.10, 0.60, 0.20];

sim_systems(3).label          = 'ESG';
sim_systems(3).short          = 'esg';
sim_systems(3).model          = 3;      % ESG BEM
sim_systems(3).noise_baseline = 1;      % uV/sqrt(Hz)  (amplifier-noise estimate)
sim_systems(3).noise_unit     = '\muV/\surdHz';
sim_systems(3).noise_unit_txt = 'uV/sqrt(Hz)';
sim_systems(3).bandwidth_hz   = Inf;    % EEG amps (e.g. Brain Products): 1 kHz
                                        % or Nyquist, whichever is lower
sim_systems(3).color          = [0.80, 0.15, 0.10];


% =========================================================================
% SIGNAL PARAMETERS — evoked response
% =========================================================================
% Mimics an evoked potential/field: a brief oscillatory burst at a fixed
% latency after stimulus onset, repeated over many trials and averaged.
%
% The source is a Gaussian-windowed sinusoid peaking at sim_evoked_latency.
% A windowed burst (rather than a sinusoid running for the whole epoch) is what
% makes this an EVOKED response: the signal is confined to a short post-stimulus
% window, and the pre- and post-burst samples carry noise only — which is
% exactly the regime a real evoked recording is scored in.

sim_fs             = 1500;    % sampling rate (Hz)
sim_duration       = 0.100;   % epoch length (s) — 100 ms post-stimulus
sim_freq           = 70;      % carrier frequency of the burst (Hz)
sim_dipole_nAm     = 3;       % PEAK source strength (nA*m)

sim_evoked_latency = 0.025;   % s — burst peaks 25 ms after trial start
sim_evoked_sigma   = 0.005;   % s — Gaussian envelope s.d. (~2 cycles at 90 Hz)

sim_time = (0 : 1/sim_fs : sim_duration - 1/sim_fs);

% Gaussian-windowed sine burst, peak amplitude sim_dipole_nAm at the latency
sim_envelope = exp(-((sim_time - sim_evoked_latency).^2) ...
                   / (2 * sim_evoked_sigma^2));
sim_waveform = sim_dipole_nAm ...
             * sin(2*pi*sim_freq*(sim_time - sim_evoked_latency)) ...
             .* sim_envelope;                                  % [1 x n_time]


% =========================================================================
% TRIAL AVERAGING
% =========================================================================
% A single 1 nA*m trial is far below every system's noise floor — r-squared
% comes out at essentially zero, which is physically correct and is precisely
% why real evoked recordings average over many trials.
%
% Averaging n independent trials leaves the (time-locked) evoked signal
% unchanged while reducing the standard deviation of the (independent) noise by
% sqrt(n). So the trial-averaged data has an effective noise level of
%
%   sigma_eff = sigma_single_trial / sqrt(sim_n_trials)
%
% 8000 trials therefore buys a sqrt(8000) ~ 89.4x reduction in noise.
%
% NOTE: the simulation applies this analytically rather than generating and
% averaging sim_n_trials separate noisy datasets. For independent Gaussian noise the
% two are exactly equivalent — the average of n iid N(0, sigma^2) draws is
% distributed as N(0, sigma^2/n) — so this is an identity, not an
% approximation, and it avoids sim_n_trials-fold computation for the same answer.

sim_n_trials = 8000;   % SET THIS: trials averaged per condition
                       % (8000 -> sqrt(8000) = 89.4x noise reduction)


% =========================================================================
% NOISE SWEEP
% =========================================================================
% Each system's noise floor is a spectral density (units/sqrt(Hz)). The noise
% power a system actually admits depends on ITS OWN measurement bandwidth
% (sim_systems(k).bandwidth_hz), capped at the simulation's Nyquist frequency:
%
%   bw_eff = min(bandwidth_hz, sim_fs/2)
%   sigma  = density * sqrt(bw_eff)         (single trial)
%
% and after averaging sim_n_trials trials:
%
%   sigma_eff = sigma / sqrt(sim_n_trials)
%
% No band-pass filter is applied. Filtering to a narrow band around sim_freq
% would shrink the noise by sqrt(bandwidth_ratio) and inflate every r-squared
% value; leaving the noise broadband keeps this a conservative, honest
% worst-case and avoids baking a filter choice into the result.

% Noise levels as multiples of each system's baseline floor. This range is
% calibrated so r^2 transitions from ~1 to ~0 at the current trial count and
% (correctly scaled) leadfield amplitudes. If curves come out flat, the usual
% cause is NOT this range but a leadfield unit-scale mismatch (sim_run_geometries
% prints peak|g| vs sigma@1x per run): if peak|signal| is orders of magnitude
% away from sigma, fix the scale rather than widening this sweep.
sim_noise_factors  = [0.125, 0.25, 0.5, 1, 2, 4, 8];   % x baseline
sim_n_realisations = 20;
sim_noise_seed     = 2026;

sim_array    = 'back';   % SET THIS: 'front' or 'back' — array used for the noise sim

% Which dipole orientations to simulate (a separate r-squared curve each)
sim_orientations = {'VD', 'RC', 'LR'};
sim_ori_display  = {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right'};


% =========================================================================
% SOURCE SELECTION
% =========================================================================

src_spacing_mm  = 5;     % mm between adjacent sources along the cord

sim_topo_src_mm = 75;   % SET THIS: cord distance (mm) for the topoplot figures
sim_focus_src_mm = 75;  % SET THIS: cord distance (mm) for the single-source r-sq curve
sim_focus_noise_factor = 0.5;   % SET THIS: noise level for the noisy topoplot (x baseline)

% Which geometry VARIANT the topoplot scripts visualise (a .name from
% sim_geometries). Leave '' to default to the unperturbed baseline.
sim_topo_geom = '';   % e.g. 'source_large' to topoplot a large source shift


% =========================================================================
% PLOT STYLING
% =========================================================================

pub_line_width  = 2.0;
pub_marker_size = 7;

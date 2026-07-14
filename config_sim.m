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
% USER CONFIGURATION — paths
% =========================================================================

% MSG and ESG geometries live in separate folders. If yours happen to share a
% folder, just set both to the same path.
msg_geoms_path = 'D:\Simulations\Pertubations\geometries';        % SET THIS
esg_geoms_path = 'D:\Simulations\Pertubations\geoms_elec';    % SET THIS

sim_save_dir   = 'D:\Simulations\Pertubations\results\simulation';% SET THIS
sim_out_dir    = 'D:\Simulations\Pertubations\fields\simulation'; % SET THIS

% Roots of the three leadfield sets (perfect / unperturbed forward fields)
msg_bem_root   = 'D:\Simulations\Pertubations\fields\mag\bem';        % SET THIS
msg_bslaw_root = 'D:\Simulations\Pertubations\fields\mag\bs_law';     % SET THIS
esg_bem_root   = 'D:\Simulations\Pertubations\fields\elec\bem_elec';    % SET THIS

% Geometry stems (WITHOUT the leading 'geometries_' prefix)
msg_geom_short = 'original_source_original';   % SET THIS
esg_geom_short = 'original_source_original';   % SET THIS

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

sim_models = struct('label', {}, 'short', {}, 'front', {}, 'back', {}, ...
                    'var', {}, 'scale', {}, 'is_meg', {}, 'geom_file', {}, ...
                    'n_axes', {}, 'axis_names', {}, 'axis_slot', {});

% --- Model 1: MSG, Biot-Savart (infinite homogeneous space — smooth fields)
sim_models(1).label     = 'MSG — Biot-Savart';
sim_models(1).short     = 'msg_bslaw';
sim_models(1).front     = fullfile(msg_bslaw_root, ...
    sprintf('leadfield_geometries_%s_bslaw_front.mat', msg_geom_short));
sim_models(1).back      = fullfile(msg_bslaw_root, ...
    sprintf('leadfield_geometries_%s_bslaw_back.mat',  msg_geom_short));
sim_models(1).var       = 'leadfield_bs';
sim_models(1).scale     = 1;
sim_models(1).is_meg    = true;
sim_models(1).geom_file = fullfile(msg_geoms_path, ...
    ['geometries_' msg_geom_short '.mat']);
sim_models(1).n_axes     = 3;
sim_models(1).axis_names = {'X-axis', 'Y-axis', 'Z-axis'};
sim_models(1).axis_slot  = [1 2 3];

% --- Model 2: MSG, BEM (individualised anatomy — sharp fields)
sim_models(2).label     = 'MSG — BEM';
sim_models(2).short     = 'msg_bem';
sim_models(2).front     = fullfile(msg_bem_root, ['geometries_' msg_geom_short], ...
    sprintf('leadfield_%s_bem_front.mat', msg_geom_short));
sim_models(2).back      = fullfile(msg_bem_root, ['geometries_' msg_geom_short], ...
    sprintf('leadfield_%s_bem_back.mat',  msg_geom_short));
sim_models(2).var       = 'leadfield_cord';
sim_models(2).scale     = msg_bem_scale;
sim_models(2).is_meg    = true;
sim_models(2).geom_file = fullfile(msg_geoms_path, ...
    ['geometries_' msg_geom_short '.mat']);
sim_models(2).n_axes     = 3;
sim_models(2).axis_names = {'X-axis', 'Y-axis', 'Z-axis'};
sim_models(2).axis_slot  = [1 2 3];

% --- Model 3: ESG, BEM (surface potentials — smooth fields)
sim_models(3).label     = 'ESG — BEM';
sim_models(3).short     = 'esg_bem';
sim_models(3).front     = fullfile(esg_bem_root, ['geometries_' esg_geom_short], ...
    sprintf('leadfield_%s_bem_front.mat', esg_geom_short));
sim_models(3).back      = fullfile(esg_bem_root, ['geometries_' esg_geom_short], ...
    sprintf('leadfield_%s_bem_back.mat',  esg_geom_short));
sim_models(3).var       = 'leadfield_cord';
sim_models(3).scale     = esg_scale;
sim_models(3).is_meg    = false;
sim_models(3).geom_file = fullfile(esg_geoms_path, ...
    ['geometries_' esg_geom_short '.mat']);
% ESG electrodes are not a Cartesian triad: two electrode sets, tangential and
% radial. n_axes MUST be stated — the ESG channel count (e.g. 342 = 2 x 171)
% can be divisible by 3, so any code that guesses the axis count from
% divisibility will wrongly read it as 3 x 114 and mis-slice the leadfield.
% SET axis_names to match the channel order your ESG leadfield was built in.
sim_models(3).n_axes     = 2;
sim_models(3).axis_names = {'Tangential', 'Radial'};
% ESG tangential lines up with the MSG X-axis, ESG radial with the MSG Z-axis.
% There is no ESG counterpart to the MSG Y-axis, so slot 2 has no ESG panel.
sim_models(3).axis_slot  = [1 3];


% =========================================================================
% SENSOR SYSTEMS (noise floors)
% =========================================================================
% SQUID and OP-MSG measure the SAME field (same array, same BEM leadfield) —
% they differ ONLY in their noise floor. That is the intended comparison:
% it isolates the effect of sensor noise from the effect of geometry.
%
% Set .model to the index of the forward model each system measures.

sim_systems = struct('label', {}, 'short', {}, 'model', {}, ...
                     'noise_baseline', {}, 'noise_unit', {}, ...
                     'noise_unit_txt', {}, 'color', {});

sim_systems(1).label          = 'SQUID MSG';
sim_systems(1).short          = 'squid_msg';
sim_systems(1).model          = 2;      % MSG BEM
sim_systems(1).noise_baseline = 5;      % fT/sqrt(Hz)  (2-5 typical)
sim_systems(1).noise_unit     = 'fT/\surdHz';
sim_systems(1).noise_unit_txt = 'fT/sqrt(Hz)';
sim_systems(1).color          = [0.10, 0.30, 0.80];

sim_systems(2).label          = 'OP-MSG';
sim_systems(2).short          = 'op_msg';
sim_systems(2).model          = 2;      % MSG BEM
sim_systems(2).noise_baseline = 20;     % fT/sqrt(Hz)  (7-20 typical)
sim_systems(2).noise_unit     = 'fT/\surdHz';
sim_systems(2).noise_unit_txt = 'fT/sqrt(Hz)';
sim_systems(2).color          = [0.10, 0.60, 0.20];

sim_systems(3).label          = 'ESG';
sim_systems(3).short          = 'esg';
sim_systems(3).model          = 3;      % ESG BEM
sim_systems(3).noise_baseline = 1;      % uV/sqrt(Hz)  (amplifier-noise estimate)
sim_systems(3).noise_unit     = '\muV/\surdHz';
sim_systems(3).noise_unit_txt = 'uV/sqrt(Hz)';
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

sim_fs             = 1000;    % sampling rate (Hz)
sim_duration       = 0.100;   % epoch length (s) — 100 ms post-stimulus
sim_freq           = 90;      % carrier frequency of the burst (Hz)
sim_dipole_nAm     = 1;       % PEAK source strength (nA*m)

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
% 2000 trials therefore buys a sqrt(2000) ~ 44.7x reduction in noise.
%
% NOTE: the simulation applies this analytically rather than generating and
% averaging 2000 separate noisy datasets. For independent Gaussian noise the
% two are exactly equivalent — the average of n iid N(0, sigma^2) draws is
% distributed as N(0, sigma^2/n) — so this is an identity, not an
% approximation, and it avoids 2000x the computation for the same answer.

sim_n_trials = 2000;   % SET THIS: trials averaged per condition


% =========================================================================
% NOISE SWEEP
% =========================================================================
% Each system's noise floor is a spectral density (units/sqrt(Hz)). Broadband
% white noise sampled at sim_fs occupies a bandwidth of sim_fs/2, so the
% SINGLE-TRIAL time-domain standard deviation is:
%
%   sigma = density * sqrt(sim_fs/2)
%
% and after averaging sim_n_trials trials:
%
%   sigma_eff = sigma / sqrt(sim_n_trials)
%
% No band-pass filter is applied. Filtering to a narrow band around sim_freq
% would shrink the noise by sqrt(bandwidth_ratio) and inflate every r-squared
% value; leaving the noise broadband keeps this a conservative, honest
% worst-case and avoids baking a filter choice into the result.

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

sim_topo_src_mm = 275;   % SET THIS: cord distance (mm) for the topoplot figures
sim_focus_src_mm = 275;  % SET THIS: cord distance (mm) for the single-source r-sq curve
sim_focus_noise_factor = 8;   % SET THIS: noise level for the noisy topoplot (x baseline)


% =========================================================================
% PLOT STYLING
% =========================================================================

pub_line_width  = 2.0;
pub_marker_size = 7;

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
%                   two electrode sets (tangential and radial). These are NOT
%                   the same quantity, so each model carries its own names and
%                   they are never implicitly paired across models.
%
%   Sensor systems (sim_systems struct array) — one entry per real system:
%     .label          - Display name, e.g. 'SQUID MSG'
%     .short          - Filename-safe stem
%     .model          - Index into sim_models: which forward model it measures
%     .noise_baseline - White noise floor, in units/sqrt(Hz)
%     .noise_unit     - Display string for the noise floor
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
% ESG BEM files are assumed to follow the BEM convention but hold potentials
% in volts, so scale 1e6 converts V/nAm -> uV/nAm. If your ESG pipeline
% already saved microvolts, set esg_scale = 1 instead.

esg_scale = 1e6;   % SET THIS: 1e6 if ESG leadfields are in V/nAm, 1 if already uV/nAm

sim_models = struct('label', {}, 'short', {}, 'front', {}, 'back', {}, ...
                    'var', {}, 'scale', {}, 'is_meg', {}, 'geom_file', {}, ...
                    'n_axes', {}, 'axis_names', {});

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

% --- Model 2: MSG, BEM (individualised anatomy — sharp fields)
sim_models(2).label     = 'MSG — BEM';
sim_models(2).short     = 'msg_bem';
sim_models(2).front     = fullfile(msg_bem_root, ['geometries_' msg_geom_short], ...
    sprintf('leadfield_%s_bem_front.mat', msg_geom_short));
sim_models(2).back      = fullfile(msg_bem_root, ['geometries_' msg_geom_short], ...
    sprintf('leadfield_%s_bem_back.mat',  msg_geom_short));
sim_models(2).var       = 'leadfield_cord';
sim_models(2).scale     = 1e15;
sim_models(2).is_meg    = true;
sim_models(2).geom_file = fullfile(msg_geoms_path, ...
    ['geometries_' msg_geom_short '.mat']);
sim_models(2).n_axes     = 3;
sim_models(2).axis_names = {'X-axis', 'Y-axis', 'Z-axis'};

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


% =========================================================================
% SENSOR SYSTEMS (noise floors)
% =========================================================================
% SQUID and OP-MSG measure the SAME field (same array, same BEM leadfield) —
% they differ ONLY in their noise floor. That is the intended comparison:
% it isolates the effect of sensor noise from the effect of geometry.
%
% Set .model to the index of the forward model each system measures.

sim_systems = struct('label', {}, 'short', {}, 'model', {}, ...
                     'noise_baseline', {}, 'noise_unit', {}, 'color', {});

sim_systems(1).label          = 'SQUID MSG';
sim_systems(1).short          = 'squid_msg';
sim_systems(1).model          = 2;      % MSG BEM
sim_systems(1).noise_baseline = 5;      % fT/sqrt(Hz)  (2-5 typical)
sim_systems(1).noise_unit     = 'fT/\surdHz';
sim_systems(1).color          = [0.10, 0.30, 0.80];

sim_systems(2).label          = 'OP-MSG';
sim_systems(2).short          = 'op_msg';
sim_systems(2).model          = 2;      % MSG BEM
sim_systems(2).noise_baseline = 20;     % fT/sqrt(Hz)  (7-20 typical)
sim_systems(2).noise_unit     = 'fT/\surdHz';
sim_systems(2).color          = [0.10, 0.60, 0.20];

sim_systems(3).label          = 'ESG';
sim_systems(3).short          = 'esg';
sim_systems(3).model          = 3;      % ESG BEM
sim_systems(3).noise_baseline = 1;      % uV/sqrt(Hz)  (amplifier-noise estimate)
sim_systems(3).noise_unit     = '\muV/\surdHz';
sim_systems(3).color          = [0.80, 0.15, 0.10];


% =========================================================================
% SIGNAL PARAMETERS
% =========================================================================

sim_fs         = 1000;   % sampling rate (Hz)
sim_duration   = 1.0;    % epoch length (s)
sim_freq       = 90;     % source oscillation frequency (Hz)
sim_dipole_nAm = 1;      % source strength (nA*m)

sim_time       = (0 : 1/sim_fs : sim_duration - 1/sim_fs);
sim_waveform   = sim_dipole_nAm * sin(2*pi*sim_freq*sim_time);   % [1 x n_time]


% =========================================================================
% NOISE SWEEP
% =========================================================================
% Each system's noise floor is a spectral density (units/sqrt(Hz)). Broadband
% white noise sampled at sim_fs occupies a bandwidth of sim_fs/2, so the
% time-domain standard deviation is:
%
%   sigma = density * sqrt(sim_fs/2)
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

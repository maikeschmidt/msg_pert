% config_pert - Shared configuration for all msg_pert perturbation scripts
%
% Defines paths, perturbation parameters, naming conventions, and plot
% styling used across the perturbation generation and analysis pipeline.
% Run as a script at the top of each msg_pert script.
%
% USAGE:
%   config_pert;
%
% VARIABLES DEFINED:
%   Paths:
%     geoms_path           - Path to original geometry .mat from msg_coreg
%     perturbed_geoms_path - Output path for perturbed geometry .mat files
%     forward_fields_base  - Path to leadfield .mat files (from msg_fwd runs)
%     save_base_dir        - Base path for saving figures and tables
%
%   Base geometry:
%     base_geom_name       - Short stem used in leadfield file names,
%                            WITHOUT the leading 'geometries_' prefix.
%                            e.g. if files are named
%                              leadfield_geometries_sub001_source_original_bslaw_back.mat
%                            set base_geom_name = 'sub001'
%
%   Source shift parameters (24 configs: 3 bundles × 8 random shifts):
%     source_shift_vectors  - {3×1} cell: each is [8×3] [dx,dy,dz] in mm.
%                             PASTE from pt_generate_source_shifts output.
%     n_source_bundles      - 3
%     n_source_shifts       - 8 (shifts per bundle)
%     sensitivity_ref_key   - Key for the original (unshifted) source model
%     sensitivity_keys      - [1 x 24] cell array of shifted model keys
%     sensitivity_labels    - [1 x 24] display labels
%     source_bundle_idx     - [1 x 24] bundle index (1, 2, or 3)
%     source_shift_idx      - [1 x 24] shift index within bundle (1–8)
%     source_bundle_display - {'~2mm (small)', '~5mm (medium)', '~10mm (large)'}
%     source_bundle_colors  - [3 x 3] RGB colours per bundle (orange family)
%
%   Sensor shift parameters (24 configs: 3 bundles × 8 random shifts):
%     sensor_shift_vectors          - {3 x 1} cell: each is [8 x 3] matrix of
%                                     [dx,dy,dz] shifts in mm. SET THESE after
%                                     generating shifts in pt_generate_sensor_shifts.m
%                                     (or paste from the printed output).
%     n_sensor_bundles              - 3
%     n_sensor_shifts               - 8 (shifts per bundle)
%     sensor_sensitivity_ref_key    - Key for original (unshifted) sensor model
%     sensor_sensitivity_keys       - [1 x 24] cell array of shifted model keys
%     sensor_sensitivity_labels     - [1 x 24] display labels
%     sensor_sensitivity_bundle_idx - [1 x 24] bundle index (1, 2, or 3)
%     sensor_sensitivity_shift_idx  - [1 x 24] shift index within bundle (1-8)
%     sensor_bundle_display         - {'~2mm (small)', '~5mm (medium)', '~10mm (large)'}
%     sensor_bundle_colors          - [3 x 3] RGB colours per bundle
%
%   Leadfield orientation labels (shared with msg_fwd):
%     orientation_labels   - {'VD', 'RC', 'LR'}
%     orientation_display  - {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right'}
%
%   Source spacing:
%     src_spacing_mm       - Source spacing along cord in mm (default: 5)
%
%   Forward model methods:
%     fwd_methods          - cell array of method names to compare
%     fwd_method_labels    - display labels (same order as fwd_methods)
%     fwd_method_colors    - [N x 3] RGB colours per method
%     fwd_method_styles    - line styles per method
%
%   Plot styling:
%     pub_line_width       - Line width for publication figures (default: 2.0)
%     pub_marker_size      - Marker size for publication figures (default: 7)
%
% NOTES:
%   - Set the four path variables and base_geom_name before running any script
%   - base_geom_name must NOT include the 'geometries_' prefix — that prefix
%     is attached automatically in file search patterns
%   - Paste sensor_shift_vectors from the output of pt_generate_sensor_shifts
%     so that displacement vs r² plots can use actual mm values on the x-axis
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
% USER CONFIGURATION — both modalities in one place (no switching)
% =========================================================================
% MSG and ESG are configured side by side here. The pipeline loops over both
% (run_perturbation_analysis), so you do NOT edit paths or flags between
% modalities. Each modality carries its own geometry + leadfield paths, method
% flags, sensor-axis count, and results (figures) folder. A third path holds the
% combined MSG-vs-ESG comparison.
%
% Make sure config_pert's own folder is on the path so pt_modality is found.
cfg_dir = fileparts(mfilename('fullpath'));
if ~isempty(cfg_dir); addpath(cfg_dir); end

base_geom_name = 'original';   % SET THIS: short stem in file names, WITHOUT the
                               % leading 'geometries_' prefix (shared by both
                               % modalities). Default dataset: 'original'.

% ---- MSG (triaxial magnetometer array) ----------------------------------
mods_cfg.msg.geoms_path           = 'D:\Simulations\Pertubations\geometries';         % SET THIS
mods_cfg.msg.perturbed_geoms_path = 'D:\Simulations\Pertubations\geometries';         % SET THIS
mods_cfg.msg.forward_fields_base  = 'D:\Simulations\Pertubations\fields\mag';         % SET THIS (holds pert_*_rsq.mat + leadfields_organised.mat)
mods_cfg.msg.bem_path             = 'D:\Simulations\Pertubations\fields\mag\bem';      % SET THIS
mods_cfg.msg.bslaw_path           = 'D:\Simulations\Pertubations\fields\mag\bs_law';  % SET THIS
mods_cfg.msg.fem_path             = '';                                                 % SET THIS
mods_cfg.msg.sphere_path          = '';                                               % SET THIS (optional)
mods_cfg.msg.bem_cond_path        = 'D:\Simulations\Pertubations\fields\mag\bem_cond_msg'; % SET THIS
mods_cfg.msg.save_base_dir        = 'D:\Simulations\Pertubations\results\msg';        % SET THIS: MSG figures/tables
mods_cfg.msg.sensor_n_axes        = 3;
mods_cfg.msg.sensor_is_meg        = true;
mods_cfg.msg.have_bem             = true;
mods_cfg.msg.have_fem             = false;
mods_cfg.msg.have_bslaw           = true;
mods_cfg.msg.have_sphere          = false;
mods_cfg.msg.have_bem_cond        = true;

% ---- ESG (tangential/radial surface electrodes) -------------------------
mods_cfg.esg.geoms_path           = 'D:\Simulations\Pertubations\geoms_elec';         % SET THIS
mods_cfg.esg.perturbed_geoms_path = 'D:\Simulations\Pertubations\geoms_elec';         % SET THIS
mods_cfg.esg.forward_fields_base  = 'D:\Simulations\Pertubations\fields\elec';        % SET THIS
mods_cfg.esg.bem_path             = 'D:\Simulations\Pertubations\fields\elec\bem_elec';    % SET THIS
mods_cfg.esg.bslaw_path           = '';                                               % Biot-Savart is magnetic — N/A for ESG
mods_cfg.esg.fem_path             = '';                                                  % SET THIS
mods_cfg.esg.sphere_path          = '';
mods_cfg.esg.bem_cond_path        = 'D:\Simulations\Pertubations\fields\elec\bem_cond_esg'; % SET THIS
mods_cfg.esg.save_base_dir        = 'D:\Simulations\Pertubations\results\esg';        % SET THIS: ESG figures/tables
mods_cfg.esg.sensor_n_axes        = 2;
mods_cfg.esg.sensor_is_meg        = false;
mods_cfg.esg.have_bem             = true;
mods_cfg.esg.have_fem             = false;
mods_cfg.esg.have_bslaw           = false;   % no Biot-Savart for ESG
mods_cfg.esg.have_sphere          = false;
mods_cfg.esg.have_bem_cond        = true;

% ---- Combined (MSG vs ESG comparison output) ----------------------------
combined_results_dir = 'D:\Simulations\Pertubations\results\combined';   % SET THIS

% Which modalities the master loop should run
pert_modalities = {'msg', 'esg'};   % SET THIS: e.g. {'msg'} to run MSG only


% =========================================================================
% ACTIVATE THE SELECTED MODALITY
% =========================================================================
% The active modality is chosen by pt_modality (set by run_perturbation_analysis
% as it loops). For standalone runs of a single script it falls back to 'msg' —
% or force it here / call pt_modality('set','esg') before running.
active_modality = pt_modality('get');
if isempty(active_modality); active_modality = 'msg'; end
if ~isfield(mods_cfg, active_modality)
    error('config_pert: no configuration for modality ''%s''.', active_modality);
end
M = mods_cfg.(active_modality);

geoms_path           = M.geoms_path;
perturbed_geoms_path = M.perturbed_geoms_path;
forward_fields_base  = M.forward_fields_base;
save_base_dir        = M.save_base_dir;
bem_path             = M.bem_path;
bslaw_path           = M.bslaw_path;
fem_path             = M.fem_path;
sphere_path          = M.sphere_path;
bem_cond_path        = M.bem_cond_path;
sensor_n_axes        = M.sensor_n_axes;
sensor_is_meg        = M.sensor_is_meg;
have_bem             = M.have_bem;
have_fem             = M.have_fem;
have_bslaw           = M.have_bslaw;
have_sphere          = M.have_sphere;
have_bem_cond        = M.have_bem_cond;

% BEM raw output is T/nAm for MSG, V/nAm for ESG — scale to the reporting unit.
% (r² and heatmaps are scale-invariant; this only affects absolute-amplitude
% outputs.) Set to 1 if your ESG pipeline already saved microvolts.
if sensor_is_meg
    bem_unit_scale = 1e15;   % T/nAm -> fT/nAm
else
    bem_unit_scale = 1e6;    % V/nAm -> uV/nAm
end


% =========================================================================
% SOURCE SHIFT PARAMETERS
% =========================================================================
% 24 configurations: 3 bundles × 8 random shifts
% Bundle 1 — small  (~2 mm):  U(1,3)  mm per axis (each axis independently)
% Bundle 2 — medium (~5 mm):  U(3,7)  mm per axis
% Bundle 3 — large  (~10 mm): U(7,13) mm per axis
% Shifts mesh_wm, mesh_bone, sources_cent.pos; torso/heart/lungs unchanged.

n_source_bundles = 3;
n_source_shifts  = 8;
n_src_total      = n_source_bundles * n_source_shifts;   % = 24

% PASTE SHIFT VECTORS HERE after running pt_generate_source_shifts
% Each cell is one bundle; each row is one [dx dy dz] shift in mm.
source_shift_vectors = {
    % Bundle 1 — ~2 mm (small): 8×3 matrix — PASTE from pt_generate_source_shifts output
     [+1.7, -2.9, -2.5;
     +1.1, -2.7, +2.2;
     -2.7, -1.4, +1.4;
     -1.9, -1.6, -2.2;
     +1.9, +2.6, -1.4;
     +2.2, +1.3, +1.1;
     -1.6, -1.2, -2.4;
     +1.1, -2.8, +1.5],...
    % Bundle 2 — ~5 mm (medium): 8×3 matrix
     [+5.2, +3.7, +6.;
     -5.4, -6.7, -3.4;
     -4.6, -4.1, +6.3;
     +3.6, +6.2, -3.3;
     +3.0, +6.3, -5.8;
     +4.4, -3.5, -6.5;
     +4.2, +4.3, -5.9;
     +3.5, +5.9, -6.0],...
    % Bundle 3 — ~10 mm (large): 8×3 matrix
      [-10.1, -9.6, +7.2;
      -8.9, -10.1, +12.4;
      -8.4, +7.5, +8.7;
      -10.8, +12.2, +11.8;
      -11.8, -12.4, -8.9;
      +11.9, -12.2, -7.0;
      -7.7, +9.0, +12.7;
      -9.2, -12.8, -12.8]
};

sensitivity_ref_key   = [base_geom_name '_source_original'];
sensitivity_keys      = cell(1, n_src_total);
sensitivity_labels    = cell(1, n_src_total);
source_bundle_idx     = zeros(1, n_src_total);
source_shift_idx      = zeros(1, n_src_total);

for b = 1:n_source_bundles
    for s = 1:n_source_shifts
        idx = (b-1)*n_source_shifts + s;
        sensitivity_keys{idx}   = sprintf('%s_source_bundle%d_shift%d', ...
            base_geom_name, b, s);
        sensitivity_labels{idx} = sprintf('Bundle %d  shift %d', b, s);
        source_bundle_idx(idx)  = b;
        source_shift_idx(idx)   = s;
    end
end

source_bundle_display = {'~2 mm (small)', '~5 mm (medium)', '~10 mm (large)'};

% Bundle colours (orange family — distinct from sensor blue family)
source_bundle_colors = [
    0.99, 0.75, 0.44;   % Bundle 1 — light orange
    0.93, 0.54, 0.13;   % Bundle 2 — orange
    0.70, 0.25, 0.05;   % Bundle 3 — dark red-orange
];


% =========================================================================
% SENSOR SHIFT PARAMETERS
% =========================================================================
% 24 configurations: 3 bundles × 8 random shifts
% Bundle 1 — small  (~2 mm):  U(1,3)  mm per axis
% Bundle 2 — medium (~5 mm):  U(3,7)  mm per axis
% Bundle 3 — large  (~10 mm): U(7,13) mm per axis

n_sensor_bundles = 3;
n_sensor_shifts  = 8;   % shifts per bundle
n_sen_total      = n_sensor_bundles * n_sensor_shifts;   % = 24

% PASTE SHIFT VECTORS HERE after running pt_generate_sensor_shifts
% Each cell is one bundle; each row is one [dx dy dz] shift in mm.
sensor_shift_vectors = {
    % Bundle 1 — ~2 mm (small): 8×3 matrix
    [ 1.7, -2.9, -2.5;
      1.1, -2.7,  2.2;
     -2.7, -1.4,  1.4;
     -1.9, -1.6, -2.2;
      1.9,  2.6, -1.4;
      2.2,  1.3,  1.1;
     -1.6, -1.2, -2.4;
      1.1, -2.8,  1.5], ...
    % Bundle 2 — ~5 mm (medium): 8×3 matrix
    [ 5.2,  3.7,  6.9;
     -5.4, -6.7, -3.4;
     -4.6, -4.1,  6.3;
      3.6,  6.2, -3.3;
      3.0,  6.3, -5.8;
      4.4, -3.5, -6.5;
      4.2,  4.3, -5.9;
      3.5,  5.9, -6.0], ...
    % Bundle 3 — ~10 mm (large): 8×3 matrix
    [-10.1,  -9.6,   7.2;
      -8.9, -10.1,  12.4;
      -8.4,   7.5,   8.7;
     -10.8,  12.2,  11.8;
     -11.8, -12.4,  -8.9;
      11.9, -12.2,  -7.0;
      -7.7,   9.0,  12.7;
      -9.2, -12.8, -12.8]
};
% Reference and shifted model keys
sensor_sensitivity_ref_key   = [base_geom_name '_sensor_original'];
sensor_sensitivity_keys      = cell(1, n_sen_total);
sensor_sensitivity_labels    = cell(1, n_sen_total);
sensor_sensitivity_bundle_idx = zeros(1, n_sen_total);
sensor_sensitivity_shift_idx  = zeros(1, n_sen_total);

for b = 1:n_sensor_bundles
    for s = 1:n_sensor_shifts
        idx = (b-1)*n_sensor_shifts + s;
        sensor_sensitivity_keys{idx}       = sprintf('%s_sensor_bundle%d_shift%d', ...
            base_geom_name, b, s);
        sensor_sensitivity_labels{idx}     = sprintf('Bundle %d  shift %d', b, s);
        sensor_sensitivity_bundle_idx(idx) = b;
        sensor_sensitivity_shift_idx(idx)  = s;
    end
end

sensor_bundle_display = {'~2 mm (small)', '~5 mm (medium)', '~10 mm (large)'};

% Bundle colours (light-to-dark blue family)
sensor_bundle_colors = [
    0.20, 0.60, 0.86;   % Bundle 1 — light blue
    0.05, 0.36, 0.65;   % Bundle 2 — mid blue
    0.00, 0.18, 0.40;   % Bundle 3 — dark blue
];


% =========================================================================
% CONDUCTIVITY PERTURBATION PARAMETERS (BEM only)
% =========================================================================
% 24 configurations: 3 bundles × 8 random perturbations of tissue σ values.
% Uses the same geometry as the unshifted source model; only the BEM
% conductivity of each compartment is scaled by a random factor.
%
% Bundle 1 — small  (up to +5%):  σ × (1 + U(0, 0.05))  per compartment
% Bundle 2 — medium (up to +10%): σ × (1 + U(0, 0.10))
% Bundle 3 — large  (up to +50%): σ × (1 + U(0, 0.50))
%
% Run run_conductivity_perturbation.m (in msg_fwd) to generate the files,
% then set have_bem_cond = true in pt_load_leadfields.

n_cond_bundles        = 3;
n_cond_shifts         = 8;
n_cond_total          = n_cond_bundles * n_cond_shifts;   % = 24
cond_bundle_pct       = [0.05, 0.10, 0.50];               % fractional range per bundle
cond_seed             = 99;    % must match seed in run_conductivity_perturbation.m
n_cond_compartments   = 5;     % SET THIS: number of BEM compartments in your geometry

cond_sensitivity_ref_key   = [base_geom_name '_source_original'];
cond_sensitivity_keys      = cell(1, n_cond_total);
cond_sensitivity_labels    = cell(1, n_cond_total);
cond_bundle_idx            = zeros(1, n_cond_total);
cond_shift_idx             = zeros(1, n_cond_total);

for b = 1:n_cond_bundles
    for s = 1:n_cond_shifts
        idx = (b-1)*n_cond_shifts + s;
        cond_sensitivity_keys{idx}   = sprintf('%s_cond_bundle%d_shift%d', ...
            base_geom_name, b, s);
        cond_sensitivity_labels{idx} = sprintf('Bundle %d  shift %d', b, s);
        cond_bundle_idx(idx)         = b;
        cond_shift_idx(idx)          = s;
    end
end

cond_bundle_display = {'+5% (small)', '+10% (medium)', '+50% (large)'};

% Bundle colours (green family — distinct from source orange and sensor blue)
cond_bundle_colors = [
    0.72, 0.92, 0.72;   % Bundle 1 — light green
    0.27, 0.68, 0.27;   % Bundle 2 — mid green
    0.10, 0.40, 0.10;   % Bundle 3 — dark green
];


% =========================================================================
% FORWARD MODEL METHODS FOR COMPARISON
% =========================================================================
% List the methods you ran in msg_fwd and want to compare.
% Labels and colours must be in the same order as fwd_methods.
% Method names must match the prefixes used in pt_load_leadfields:
%   'bslaw'  — Biot-Savart law (infinite homogeneous medium)
%   'sphere' — Single sphere (Sarvas analytical)
%   'bem'    — Boundary Element Method
%   'fem'    — Finite Element Method

fwd_methods       = {'bslaw', 'bem'};   % SET THIS: methods to compare
fwd_method_labels = {'Biot-Savart','BEM'};   % SET THIS

fwd_method_colors = [
    0.80, 0.15, 0.10;   % bslaw  — red
    0.10, 0.30, 0.80;   % sphere — blue
    0.10, 0.60, 0.20;   % bem    — green  
    0.55, 0.10, 0.75;   % fem    — purple 
];

fwd_method_styles = {'-', '--', ':', '-.'};   % line style per method


% =========================================================================
% LEADFIELD ORIENTATION LABELS
% =========================================================================

orientation_labels  = {'VD', 'RC', 'LR'};
orientation_display = {'Ventral-Dorsal', 'Rostral-Caudal', 'Left-Right'};


% =========================================================================
% SOURCE SPACING AND PLOT STYLING
% =========================================================================

src_spacing_mm   = 5;    % mm between adjacent source positions along cord

pub_line_width   = 2.0;
pub_marker_size  = 7;
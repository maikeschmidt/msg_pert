% config_pert - Shared configuration for all msg_pert perturbation scripts
%
% Defines paths, perturbation parameters, naming conventions, and plot
% styling used across the perturbation generation and analysis pipeline.
% Run as a script at the top of each msg_pert script.
%
% USAGE:
%   config_pert;
%
% MODALITIES:
%   MSG (triaxial magnetometers) and ESG (surface electrodes) are configured
%   side by side in mods_cfg. The active one is chosen by pt_modality, which
%   run_perturbation_analysis sets per iteration, so nothing has to be edited
%   between modalities. The selected modality's fields are then unpacked into
%   the loose variables listed below.
%
% VARIABLES DEFINED:
%   Paths (per modality):
%     geoms_path           - Path to original geometry .mat from msg_coreg
%     perturbed_geoms_path - Output path for perturbed geometry .mat files
%     forward_fields_base  - Path to leadfield .mat files (from msg_fwd runs)
%     save_base_dir        - Base path for saving figures and tables
%     combined_results_dir - Output for the combined MSG-vs-ESG comparison
%
%   Sensor description and method availability (per modality):
%     sensor_n_axes        - 3 for triaxial MSG, 2 for ESG. Declared, not
%                            inferred: an electrode count can also divide by 3
%     sensor_is_meg        - true for MSG, false for ESG
%     have_bem / have_fem / have_bslaw / have_sphere / have_bem_cond
%                          - which forward models were computed in msg_fwd
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
%   Loading (both modalities):
%     analysis_array       - 'back' or 'front': the only array loaded
%     unit_scale_mode      - 'auto' (per-file, from magnitude) or 'fixed'
%     mods_cfg.<m>.axis_names / axis_slot / radial_axis / sensor_shift_dims
%                          - sensor axes and how they correspond across modalities
%
%   Staged analysis (stages 1-6, see RUN_ORDER.md):
%     stage_results_dir    - output root; one numbered folder per stage
%     pert_types, pert_type_display, pert_type_colors, bundle_names
%     metric_orientations, headline_orientation, ori_titles
%     stats_n_perm, stats_n_boot, stats_alpha, stats_seed
%     source_sensor_paired - source and sensor shift k share a shift vector
%     msg_esg_paired       - per type: is shift k the same change in MSG and ESG
%     modality_compare_method - forward model used for MSG vs ESG
%     baseline_models, baseline_pairs, baseline_topo_src_mm   (stage 1)
%     noise_systems, noise_* waveform and sweep settings        (stage 5)
%     fig_resolution       - PNG resolution (dpi)
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
mods_cfg.msg.display              = 'MSG';
mods_cfg.msg.axis_names           = {'X', 'Y', 'Z'};
mods_cfg.msg.axis_slot            = [1 2 3];   % comparison slot per axis (see below)
mods_cfg.msg.radial_axis          = 3;
mods_cfg.msg.sensor_shift_dims    = [1 2 3];   % axes the sensor shift was applied in
mods_cfg.msg.unit                 = 'fT/nAm';

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
mods_cfg.esg.display              = 'ESG';
mods_cfg.esg.axis_names           = {'Tangential', 'Radial'};
mods_cfg.esg.axis_slot            = [1 3];     % tangential <-> MSG X, radial <-> MSG Z
mods_cfg.esg.radial_axis          = 2;
mods_cfg.esg.sensor_shift_dims    = [1 2];     % ESG sensor shifts hold Z at 0
mods_cfg.esg.unit                 = 'uV/nAm';

% AXIS SLOTS. MSG and ESG do not have the same sensor axes, so an MSG-vs-ESG
% comparison needs to know which axes correspond. Slots follow the MSG
% triaxial convention (1 = X, 2 = Y, 3 = Z); an axis with the same slot in
% both modalities is compared directly. ESG tangential sits with MSG X and ESG
% radial with MSG Z. Every comparison is ALSO made on the whole array (all
% axes stacked), which needs no matching and is the headline cross-modality
% number.

% ---- Combined (MSG vs ESG comparison output) ----------------------------
combined_results_dir = 'D:\Simulations\Pertubations\results\combined';   % SET THIS

% Which modalities the master loop should run
pert_modalities = {'msg', 'esg'};   % SET THIS: e.g. {'msg'} to run MSG only

% Which sensor array the analysis uses. msg_fwd writes a front and a back
% file per geometry; only files for this array are loaded. Loading both
% under one key would let whichever is read second silently replace the
% first.
analysis_array = 'back';            % SET THIS: 'back' or 'front'

% How leadfields are brought to fT/nAm (MSG) and uV/nAm (ESG) on load.
%   'auto'  - per file, from its magnitude (pt_unit_scale). Needed whenever
%             files come from different writers: the conductivity files carry
%             an extra x1e15, and raw BEM output may be per nA*m or per A*m.
%   'fixed' - bem_unit_scale below for BEM, 1 for everything else (the old
%             behaviour; safe for r2 only, which ignores scale).
unit_scale_mode = 'auto';


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
axis_names           = M.axis_names;
radial_axis          = M.radial_axis;

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
% then set have_bem_cond = true for that modality in the block above.

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

% Methods to compare, PER MODALITY. ESG has no Biot-Savart (that is a magnetic
% model), so the ESG list is BEM only — otherwise the analysis would look for
% bslaw_ ESG leadfields that do not exist.
if sensor_is_meg
    fwd_methods       = {'bslaw', 'bem'};   % SET THIS: MSG methods
    fwd_method_labels = {'Biot-Savart', 'BEM'};
else
    fwd_methods       = {'bem'};            % ESG: BEM only
    fwd_method_labels = {'BEM'};
end

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

% =========================================================================
% STAGED ANALYSIS (stages 1-6, see RUN_ORDER.md)
% =========================================================================
% Everything below is read by the staged scripts: pt_compute_metrics,
% pt_baseline, pt_within_modality, pt_compare_modalities, pt_noise_simulate,
% pt_noise_analyse and pt_summary_table. Each stage writes to its own
% numbered subfolder of stage_results_dir.

stage_results_dir = 'D:\Simulations\Pertubations\results\staged';   % SET THIS

% Perturbation types, in reporting order, and how they are labelled
pert_types        = {'source', 'sensor', 'cond'};
pert_type_display = struct('source', 'Source space', ...
                           'sensor', 'Sensor array', ...
                           'cond',   'Conductivity');
pert_type_colors  = struct('source', [0.90 0.55 0.10], ...   % orange
                           'sensor', [0.20 0.45 0.80], ...   % blue
                           'cond',   [0.25 0.62 0.35]);      % green
bundle_names      = {'small', 'medium', 'large'};

% Dipole orientations reported by every stage. 'ALL' stacks the three
% orientations into one vector (msg_fwd's concatenated convention) and is the
% single number to quote when one is wanted per comparison.
metric_orientations  = {'VD', 'RC', 'LR', 'ALL'};
headline_orientation = 'ALL';
ori_titles = struct('VD', 'Ventral-Dorsal', 'RC', 'Rostral-Caudal', ...
                    'LR', 'Left-Right', 'ALL', 'All orientations');

% Modality display colours
modality_colors = struct('msg', [0.49 0.18 0.56], ...   % purple
                         'esg', [0.20 0.63 0.35]);      % green

% ---- Statistics ---------------------------------------------------------
% The unit of observation is one perturbation REALISATION (one random shift,
% or one random conductivity draw), summarised by its median over cord
% positions. With 8 realisations per bundle, bundle-level tests are small-n
% by design; the pooled rows (24 per type) carry more power.
stats_n_perm = 10000;   % permutations (tests between 8 vs 8 are enumerated exactly)
stats_n_boot = 10000;   % bootstrap draws for median CIs
stats_alpha  = 0.05;    % FDR level
stats_seed   = 2026;    % RNG seed for reproducible p-values and CIs

% Source and sensor shifts are drawn from the SAME shift vectors (see the two
% blocks above), so shift k of one is paired with shift k of the other. Set
% false if you regenerate them independently.
source_sensor_paired = true;

% Is shift k in MSG the same geometric change as shift k in ESG? True when
% both modalities were generated from the same shift vectors / conductivity
% seed. ESG sensor shifts hold Z at 0, so they are NOT the same shift.
msg_esg_paired = struct('source', true, 'sensor', false, 'cond', true);

% Forward model used for MSG-vs-ESG comparisons (must exist in both)
modality_compare_method = 'bem';

% ---- Stage 1: noise-free baseline ---------------------------------------
% {modality, method, display name}. The unperturbed geometry of each is shown.
baseline_models = {
    'msg', 'bslaw', 'MSG — Biot-Savart'
    'msg', 'bem',   'MSG — BEM'
    'esg', 'bem',   'ESG — BEM'
};
baseline_model_colors = [0.80 0.15 0.10; 0.49 0.18 0.56; 0.20 0.63 0.35];

% Forward-model-type comparisons: {modality, reference method, comparison
% method, label}. The reference is the RE denominator.
baseline_pairs = {
    'msg', 'bem', 'bslaw', 'Biot-Savart vs BEM'
};

% Cord positions (mm from the top of the source space) for the topoplots
baseline_topo_src_mm = [75, 250, 450];   % SET THIS

% ---- Stage 5: sensor noise ----------------------------------------------
% The measured field pattern is estimated from trial-averaged data by
% projecting onto the known source waveform w:
%     g_hat = Y w / (w'w) = g + noise,   noise s.d. = sigma / ||w||
% and scored with the same four metrics against the noise-free ORIGINAL
% field. At zero noise this reproduces stages 2-4 exactly.
%
% Noise floors and bandwidths match simulations/config_sim.m; keep the two
% in step if you change either.
noise_systems = struct( ...
    'label',        {'SQUID MSG',  'OP-MSG',     'ESG'}, ...
    'short',        {'squid_msg',  'op_msg',     'esg'}, ...
    'modality',     {'msg',        'msg',        'esg'}, ...
    'method',       {'bem',        'bem',        'bem'}, ...
    'density',      {5,            20,           1}, ...        % per sqrt(Hz)
    'unit_txt',     {'fT/sqrt(Hz)', 'fT/sqrt(Hz)', 'uV/sqrt(Hz)'}, ...
    'bandwidth_hz', {1000,         150,          Inf}, ...
    'color',        {[0.10 0.30 0.80], [0.10 0.60 0.20], [0.80 0.15 0.10]});

% Evoked source waveform: Gaussian-windowed sinusoid (as config_sim)
noise_fs        = 1500;     % Hz
noise_duration  = 0.100;    % s
noise_freq      = 70;       % Hz
noise_peak_nAm  = 3;        % nA*m
noise_latency   = 0.025;    % s
noise_env_sd    = 0.005;    % s
noise_n_trials  = 8000;     % trials averaged

% Noise levels as multiples of each system's own floor. 0 is always added,
% so the noise-free result is the first point of every curve.
noise_factors   = [0.125, 0.25, 0.5, 1, 2, 4, 8];
noise_n_real    = 20;       % noise realisations per level
noise_seed      = 2026;
noise_reference_factor = 1; % level used for the "realistic" comparisons

% A perturbation counts as detected in one realisation when its error
% exceeds this percentile of the unperturbed (noise-only) error at the
% same noise level. The critical noise level is where the detection rate
% falls below noise_detect_rate.
noise_detect_pct  = 95;
noise_detect_rate = 0.8;

% Systems compared against each other at every noise level
noise_compare_pairs = {'squid_msg', 'esg'; 'op_msg', 'esg'; 'squid_msg', 'op_msg'};

% Figure resolution (dpi)
fig_resolution = 600;

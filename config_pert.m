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
%   Source shift parameters (18 configs: ±2,±4,±6 mm × X,Y,Z):
%     source_shift_mm      - Shift magnitudes, e.g. [2 4 6] mm
%     source_shift_axes    - Axis labels: {'X','Y','Z'}
%     sensitivity_ref_key  - Key for the original (unshifted) source model
%     sensitivity_keys     - [1 x 18] cell array of shifted model keys
%     sensitivity_labels   - [1 x 18] display labels matching sensitivity_keys
%     sensitivity_markers  - [1 x 18] marker styles
%     sensitivity_styles   - [1 x 18] line styles ('+' = solid, '-' = dashed)
%     sensitivity_shift_axis - [1 x 18] axis index (1=X, 2=Y, 3=Z)
%     sensitivity_axis_colors - [3 x 3] RGB colours per shift axis
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
%   Plot styling:
%     pub_line_width       - Line width for publication figures (default: 2.0)
%     pub_marker_size      - Marker size for publication figures (default: 7)
%
% NOTES:
%   - Set the four path variables and base_geom_name before running any script
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
% USER CONFIGURATION — set these paths before running any script
% =========================================================================

geoms_path           = 'D:\Simulations\Pertubations\geometries';   % SET THIS: path to original geometry .mat file
perturbed_geoms_path = 'D:\Simulations\Pertubations\geometries';   % SET THIS: output path for perturbed geometry files
forward_fields_base  = 'D:\Simulations\Pertubations\fields';   % SET THIS: path to leadfield .mat files (from msg_fwd)
save_base_dir        = 'D:\Simulations\Pertubations\results';   % SET THIS: base path for figures and tables

base_geom_name       = 'original';   % SET THIS: short stem used in file names,
                             %   WITHOUT the leading 'geometries_' prefix.
                             %   e.g. if files are named
                             %     leadfield_geometries_sub001_source_original_bslaw_back.mat
                             %   set base_geom_name = 'sub001'
                             %   For the default dataset: base_geom_name = 'original'


% =========================================================================
% SOURCE SHIFT PARAMETERS
% =========================================================================
% 18 configurations: ±2, ±4, ±6 mm × X, Y, Z axes

source_shift_mm   = [2 4 6];          % magnitudes in mm
source_shift_axes = {'X', 'Y', 'Z'};  % axis labels

% Expanded signed shifts: +2,+4,+6,-2,-4,-6 per axis, ordered X,Y,Z
signs  = [+1 +1 +1 -1 -1 -1];
mags   = [2   4   6  2   4   6];
n_src_shifts   = 18;
n_per_axis     = 6;

% Auto-build key, label, marker, style, axis-index arrays
sensitivity_ref_key   = [base_geom_name '_source_original'];
sensitivity_keys      = cell(1, n_src_shifts);
sensitivity_labels    = cell(1, n_src_shifts);
sensitivity_markers   = {'o','s','^','o','s','^', ...
                          'o','s','^','o','s','^', ...
                          'o','s','^','o','s','^'};
sensitivity_styles    = {'-', '-', '-', '--', '--', '--', ...
                          '-', '-', '-', '--', '--', '--', ...
                          '-', '-', '-', '--', '--', '--'};
sensitivity_shift_axis = [ones(1,6) 2*ones(1,6) 3*ones(1,6)];

for ax_i = 1:3
    for sh_i = 1:6
        idx = (ax_i-1)*6 + sh_i;
        sign_char = 'p'; if signs(sh_i) < 0, sign_char = 'n'; end
        sensitivity_keys{idx}   = sprintf('%s_source_%s_%s%dmm', ...
            base_geom_name, source_shift_axes{ax_i}, sign_char, mags(sh_i));
        sensitivity_labels{idx} = sprintf('%s%+dmm', ...
            source_shift_axes{ax_i}, signs(sh_i)*mags(sh_i));
    end
end

% Axis colours (X=blue, Y=orange, Z=green — colour-blind safe)
sensitivity_axis_colors = [
    0.00, 0.45, 0.70;   % X — blue
    0.90, 0.62, 0.00;   % Y — orange
    0.00, 0.62, 0.45;   % Z — green
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
sensor_shift_vectors = {};   % SET THIS (paste from pt_generate_sensor_shifts output)

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
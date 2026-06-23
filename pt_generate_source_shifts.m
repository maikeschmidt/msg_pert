% pt_generate_source_shifts - Generate geometry files for source-space perturbations
%
% Applies systematic shifts to the spinal cord mesh (mesh_wm), bone mesh
% (mesh_bone), and source positions (sources_cent.pos) in the original
% geometry struct. Produces 18 geometry files: ±2, ±4, ±6 mm along each
% of the X, Y, and Z anatomical axes.
%
% The torso, heart, and lung meshes are NOT shifted — only the structures
% that move with the spinal cord are perturbed. This models the effect of
% uncertainty or variability in anatomical understanding of the cord/bone
% relative to the fixed sensor array and outer body surface.
%
% Saves one .mat file per shift and prints the complete list of filenames
% for pasting into msg_fwd.
%
% USAGE:
%   pt_generate_source_shifts
%   pt_generate_source_shifts(S)
%
% INPUT (set in script or pass via S):
%   S.geoms_path           - Path to original geometry .mat file
%   S.perturbed_geoms_path - Output directory for shifted geometry files
%   S.base_geom_name       - Stem of original geometry file
%   S.source_shift_mm      - Shift magnitudes in mm (default: [2 4 6])
%
% OUTPUTS (saved to S.perturbed_geoms_path):
%   geometries_<base>_source_original.mat           — unshifted reference copy
%   geometries_<base>_source_X_p2mm.mat             — +2 mm along X
%   geometries_<base>_source_X_n2mm.mat             — -2 mm along X
%   ... (18 total shifted files + 1 original copy)
%
% ALSO PRINTS:
%   List of all geometry file stems to paste into msg_fwd's config_models.m
%   (or run_bem_leadfields.m) to trigger forward model computation.
%
% DEPENDENCIES:
%   config_pert
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
% Date:   June 2026

function pt_generate_source_shifts(S)

if nargin < 1, S = struct(); end

config_pert;

if ~isfield(S, 'geoms_path'),           S.geoms_path           = geoms_path;           end
if ~isfield(S, 'perturbed_geoms_path'), S.perturbed_geoms_path = perturbed_geoms_path; end
if ~isfield(S, 'base_geom_name'),       S.base_geom_name       = base_geom_name;       end
if ~isfield(S, 'source_shift_mm'),      S.source_shift_mm      = source_shift_mm;      end

if isempty(S.geoms_path) || isempty(S.base_geom_name)
    error('Set S.geoms_path and S.base_geom_name (or configure config_pert.m).');
end
if isempty(S.perturbed_geoms_path)
    error('Set S.perturbed_geoms_path for output geometry files.');
end

if ~isfolder(S.perturbed_geoms_path), mkdir(S.perturbed_geoms_path); end

% Load original geometry
geom_file = fullfile(S.geoms_path, [S.base_geom_name '.mat']);
if ~isfile(geom_file)
    error('Original geometry file not found: %s', geom_file);
end
geom = load(geom_file);
fprintf('\nLoaded: %s\n', geom_file);
fprintf('  mesh_wm:   %d vertices\n', size(geom.mesh_wm.vertices, 1));
fprintf('  mesh_bone: %d vertices\n', size(geom.mesh_bone.vertices, 1));
fprintf('  sources:   %d positions\n', size(geom.sources_cent.pos, 1));

% Save unshifted reference copy
ref_name    = [S.base_geom_name '_source_original'];
ref_outfile = fullfile(S.perturbed_geoms_path, ['geometries_' ref_name '.mat']);
save(ref_outfile, '-struct', 'geom', '-v7.3');
fprintf('\nSaved reference (unshifted): %s\n', ref_outfile);

% Build shift list: ±2,±4,±6 mm × X,Y,Z
axes     = {'X', 'Y', 'Z'};
signs    = [+1  +1  +1  -1  -1  -1];
mags     = [2    4   6   2   4   6 ];
saved_names = {ref_name};   % collect for filename list

fprintf('\n=== Generating source shift geometries ===\n');
fprintf('  Shifting: mesh_wm, mesh_bone, sources_cent.pos\n');
fprintf('  Fixed:    mesh_torso, mesh_heart, mesh_lungs\n\n');

for ax_i = 1:3
    for sh_i = 1:6
        delta      = zeros(1, 3);
        delta(ax_i) = signs(sh_i) * mags(sh_i);

        sign_char   = 'p'; if signs(sh_i) < 0, sign_char = 'n'; end
        model_name  = sprintf('%s_source_%s_%s%dmm', ...
            S.base_geom_name, axes{ax_i}, sign_char, mags(sh_i));
        outfile     = fullfile(S.perturbed_geoms_path, ...
            ['geometries_' model_name '.mat']);

        if isfile(outfile)
            fprintf('  Already exists: %s — skipping.\n', model_name);
            saved_names{end+1} = model_name; %#ok<AGROW>
            continue;
        end

        % Apply shift to cord, bone, and sources only
        geom_shifted = geom;
        geom_shifted.mesh_wm.vertices      = geom.mesh_wm.vertices      + delta;
        geom_shifted.mesh_bone.vertices    = geom.mesh_bone.vertices    + delta;
        geom_shifted.sources_cent.pos      = geom.sources_cent.pos      + delta;

        save(outfile, '-struct', 'geom_shifted', '-v7.3');
        saved_names{end+1} = model_name; %#ok<AGROW>

        fprintf('  [%s %s%dmm]  delta=[%+.0f, %+.0f, %+.0f] mm  → %s\n', ...
            axes{ax_i}, sign_char, mags(sh_i), delta(1), delta(2), delta(3), ...
            model_name);
    end
end

% Print filename list for msg_fwd
fprintf('\n\n');
fprintf('=================================================================\n');
fprintf('  PASTE THESE FILENAMES INTO msg_fwd config_models.m\n');
fprintf('  (or into the filenames cell in run_bem_leadfields.m)\n');
fprintf('=================================================================\n');
fprintf('\nfilenames = {\n');
for k = 1:numel(saved_names)
    fprintf("    '%s', ...\n", saved_names{k});
end
fprintf('};\n\n');
fprintf('Total geometries: %d (1 original + %d shifted)\n\n', ...
    numel(saved_names), numel(saved_names)-1);

end

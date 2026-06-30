% pt_generate_source_shifts - Generate geometry files for source-space perturbations
%
% Applies random 3-D shifts to the spinal cord mesh (mesh_wm), bone mesh
% (mesh_bone), and source positions (sources_cent.pos) in the original
% geometry struct. Produces 24 geometry files across 3 error bundles
% (8 random shifts each). Each axis is shifted independently so there is
% no preferred direction.
%
% Shift bundle definitions:
%   Bundle 1 — small  (~2 mm):  U(1,3)  mm per axis + random sign
%   Bundle 2 — medium (~5 mm):  U(3,7)  mm per axis + random sign
%   Bundle 3 — large  (~10 mm): U(7,13) mm per axis + random sign
%
% The torso, heart, and lung meshes are NOT shifted — only the structures
% that move with the spinal cord are perturbed. This models uncertainty or
% variability in the anatomical location of the cord/bone relative to the
% fixed outer body surface and sensor array.
%
% Shift vectors are randomly generated and printed so they can be pasted
% back into config_pert.m as source_shift_vectors.
%
% USAGE:
%   pt_generate_source_shifts
%   pt_generate_source_shifts(S)
%
% INPUT (set in script or pass via S):
%   S.geoms_path           - Path to original geometry .mat file
%   S.perturbed_geoms_path - Output directory for shifted geometry files
%   S.base_geom_name       - Stem of original geometry file
%   S.seed                 - Random seed for reproducible shifts (default: 42)
%
% OUTPUTS (saved to S.perturbed_geoms_path):
%   geometries_<base>_source_original.mat           — unshifted reference copy
%   geometries_<base>_source_bundle1_shift1.mat     — Bundle 1, shift 1
%   ...
%   geometries_<base>_source_bundle3_shift8.mat     — Bundle 3, shift 8
%
% ALSO PRINTS:
%   Generated shift vectors (for pasting into config_pert.m)
%   List of all geometry file stems for pasting into msg_fwd
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
if ~isfield(S, 'seed'),                 S.seed                 = 42;                   end

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
fprintf('\nSaved reference (unshifted): %s\n\n', ref_outfile);

% Bundle shift range definitions [lower, upper] mm per axis (absolute)
bundle_ranges = [1, 3; 3, 7; 7, 13];
n_bundles = 3;
n_shifts  = 8;

% Generate shift vectors (each axis independently randomised)
rng(S.seed);
shift_vectors = cell(n_bundles, 1);

for b = 1:n_bundles
    lo = bundle_ranges(b, 1);
    hi = bundle_ranges(b, 2);
    vecs = zeros(n_shifts, 3);
    for s = 1:n_shifts
        mag   = lo + (hi - lo) * rand(1, 3);       % magnitude per axis
        signs = (rand(1, 3) > 0.5) * 2 - 1;        % random sign per axis
        vecs(s, :) = mag .* signs;
    end
    shift_vectors{b} = vecs;
end

% Print shift vectors for pasting into config_pert
bundle_names = {'Bundle 1 — small  (~2mm)', ...
                'Bundle 2 — medium (~5mm)', ...
                'Bundle 3 — large  (~10mm)'};
fprintf('=================================================================\n');
fprintf('  PASTE THIS INTO config_pert.m  →  source_shift_vectors\n');
fprintf('=================================================================\n');
fprintf('source_shift_vectors = {\n');
for b = 1:n_bundles
    fprintf('    %% %s\n', bundle_names{b});
    fprintf('    [');
    for s = 1:n_shifts
        v = shift_vectors{b}(s, :);
        if s < n_shifts
            fprintf('%+.2f, %+.2f, %+.2f;\n     ', v(1), v(2), v(3));
        else
            fprintf('%+.2f, %+.2f, %+.2f], ...\n', v(1), v(2), v(3));
        end
    end
end
fprintf('};\n\n');

% Generate and save shifted geometry files
fprintf('=== Generating source shift geometries ===\n');
fprintf('  Shifting: mesh_wm, mesh_bone, sources_cent.pos\n');
fprintf('  Fixed:    mesh_torso, mesh_heart, mesh_lungs\n\n');

saved_names = {ref_name};

for b = 1:n_bundles
    for s = 1:n_shifts
        dxyz = shift_vectors{b}(s, :);

        model_name = sprintf('%s_source_bundle%d_shift%d', S.base_geom_name, b, s);
        outfile    = fullfile(S.perturbed_geoms_path, ['geometries_' model_name '.mat']);

        if isfile(outfile)
            fprintf('  Already exists: %s — skipping.\n', model_name);
            saved_names{end+1} = model_name; %#ok<AGROW>
            continue;
        end

        geom_shifted = geom;
        geom_shifted.mesh_wm.vertices   = geom.mesh_wm.vertices   + dxyz;
        geom_shifted.mesh_bone.vertices = geom.mesh_bone.vertices + dxyz;
        geom_shifted.sources_cent.pos   = geom.sources_cent.pos   + dxyz;

        save(outfile, '-struct', 'geom_shifted', '-v7.3');
        saved_names{end+1} = model_name; %#ok<AGROW>

        fprintf('  [Bundle %d  Shift %d]  [%+.1f, %+.1f, %+.1f] mm  → %s\n', ...
            b, s, dxyz(1), dxyz(2), dxyz(3), model_name);
    end
end

% Print filename list for msg_fwd
fprintf('\n');
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

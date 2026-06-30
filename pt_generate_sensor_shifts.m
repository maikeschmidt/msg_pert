% pt_generate_sensor_shifts - Generate geometry files for sensor-array perturbations
%
% Applies random rigid-body shifts to the sensor array in the original
% geometry struct. Produces 24 geometry files across 3 error bundles
% (8 random shifts each). All anatomy meshes remain fixed — only sensor
% positions (coilpos / chanpos) are displaced.
%
% Shift bundle definitions:
%   Bundle 1 — small  (~2 mm):  U(1,3)  mm per axis + random sign
%   Bundle 2 — medium (~5 mm):  U(3,7)  mm per axis + random sign
%   Bundle 3 — large  (~10 mm): U(7,13) mm per axis + random sign
%
% Shift vectors are randomly generated and then printed so they can be
% pasted back into config_pert.m for use in the displacement-vs-r² plots.
%
% USAGE:
%   pt_generate_sensor_shifts
%   pt_generate_sensor_shifts(S)
%
% INPUT (set in script or pass via S):
%   S.geoms_path           - Path to original geometry .mat file
%   S.perturbed_geoms_path - Output directory for shifted geometry files
%   S.base_geom_name       - Stem of original geometry file
%   S.seed                 - Random seed for reproducible shifts (default: 42)
%
% OUTPUTS (saved to S.perturbed_geoms_path):
%   geometries_<base>_sensor_original.mat              — unshifted reference
%   geometries_<base>_sensor_bundle1_shift1.mat        — Bundle 1, shift 1
%   ...
%   geometries_<base>_sensor_bundle3_shift8.mat        — Bundle 3, shift 8
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

function pt_generate_sensor_shifts(S)

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

% Detect sensor arrays in geometry
arrays = detect_sensor_arrays(geom, S.base_geom_name);
if isempty(arrays)
    error('No sensor arrays found in geometry: %s', S.base_geom_name);
end
fprintf('  Detected sensor arrays: ');
for a = 1:numel(arrays); fprintf('%s  ', arrays{a}.label); end
fprintf('\n');

% Save unshifted reference copy
ref_name    = [S.base_geom_name '_sensor_original'];
ref_outfile = fullfile(S.perturbed_geoms_path, ['geometries_' ref_name '.mat']);
save(ref_outfile, '-struct', 'geom', '-v7.3');
fprintf('Saved reference (unshifted): %s\n\n', ref_outfile);

% Bundle shift range definitions [lower, upper] mm per axis (absolute)
bundle_ranges = [1, 3; 3, 7; 7, 13];

% Generate shift vectors
rng(S.seed);
n_bundles = 3;
n_shifts  = 8;
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
fprintf('=================================================================\n');
fprintf('  PASTE THIS INTO config_pert.m  →  sensor_shift_vectors\n');
fprintf('=================================================================\n');
fprintf('sensor_shift_vectors = {\n');
bundle_names = {'Bundle 1 — small  (~2mm)',  ...
                'Bundle 2 — medium (~5mm)', ...
                'Bundle 3 — large  (~10mm)'};
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
fprintf('=== Generating sensor shift geometries ===\n');
fprintf('  Shifting: sensor array coilpos / chanpos\n');
fprintf('  Fixed:    all anatomy meshes\n\n');

saved_names = {ref_name};

for b = 1:n_bundles
    for s = 1:n_shifts
        dxyz = shift_vectors{b}(s, :);

        model_name = sprintf('%s_sensor_bundle%d_shift%d', ...
            S.base_geom_name, b, s);
        outfile = fullfile(S.perturbed_geoms_path, ...
            ['geometries_' model_name '.mat']);

        if isfile(outfile)
            fprintf('  Already exists: %s — skipping.\n', model_name);
            saved_names{end+1} = model_name; %#ok<AGROW>
            continue;
        end

        geom_shifted = geom;
        for a = 1:numel(arrays)
            field  = arrays{a}.field;
            grad_s = geom_shifted.(field);
            % ESG (surface electrodes): constrain to XY plane only — Z held fixed.
            % Detect by content: ESG has elecpos, MSG/OPM has coilpos.
            shift = dxyz;
            if isfield(grad_s, 'elecpos') && ~isfield(grad_s, 'coilpos')
                shift(3) = 0;
            end
            if isfield(grad_s, 'coilpos')
                grad_s.coilpos = grad_s.coilpos + shift;
            elseif isfield(grad_s, 'elecpos')
                grad_s.elecpos = grad_s.elecpos + shift;
            end
            if isfield(grad_s, 'chanpos')
                grad_s.chanpos = grad_s.chanpos + shift;
            end
            geom_shifted.(field) = grad_s;
        end

        save(outfile, '-struct', 'geom_shifted', '-v7.3');
        saved_names{end+1} = model_name; %#ok<AGROW>

        fprintf('  [Bundle %d  Shift %d]  [%.1f, %.1f, %.1f] mm  → %s\n', ...
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


% ---- Local function ----

function arrays = detect_sensor_arrays(geom, model_name)
% Detect sensor arrays in geometry. Returns struct array with fields:
%   .field   — fieldname in geom struct
%   .label   — 'front', 'back', or 'experimental'
%   .xy_only — true for ESG (sensors_2axis); shift constrained to XY plane
%
% ESG surface electrodes use sensors_2axis fields — shifts in X/Y only.
% MSG/OPM coils use coils_3axis / coils_2axis — full 3D shifts.
    arrays = {};

    % Experimental / custom single array
    if isfield(geom, 'experimental_sensors')
        arrays{end+1} = struct('field', 'experimental_sensors', ...
            'label', 'experimental', 'xy_only', false);
    end

    % Front array — check MSG (3axis/2axis coils) then ESG (sensors_2axis)
    if isfield(geom, 'front_coils_3axis')
        arrays{end+1} = struct('field', 'front_coils_3axis', ...
            'label', 'front', 'xy_only', false);
    elseif isfield(geom, 'front_coils_2axis')
        arrays{end+1} = struct('field', 'front_coils_2axis', ...
            'label', 'front', 'xy_only', false);
    elseif isfield(geom, 'front_sensors_2axis')
        arrays{end+1} = struct('field', 'front_sensors_2axis', ...
            'label', 'front', 'xy_only', true);
    elseif isfield(geom, 'front_sensors')
        arrays{end+1} = struct('field', 'front_sensors', ...
            'label', 'front', 'xy_only', false);
    end

    % Back array
    if isfield(geom, 'back_coils_3axis')
        arrays{end+1} = struct('field', 'back_coils_3axis', ...
            'label', 'back', 'xy_only', false);
    elseif isfield(geom, 'back_coils_2axis')
        arrays{end+1} = struct('field', 'back_coils_2axis', ...
            'label', 'back', 'xy_only', false);
    elseif isfield(geom, 'back_sensors_2axis')
        arrays{end+1} = struct('field', 'back_sensors_2axis', ...
            'label', 'back', 'xy_only', true);
    elseif isfield(geom, 'back_sensors')
        arrays{end+1} = struct('field', 'back_sensors', ...
            'label', 'back', 'xy_only', false);
    end

    % Also scan for any field containing 'sensors_2axis' not caught above
    all_fields = fieldnames(geom);
    for fi = 1:numel(all_fields)
        fn = all_fields{fi};
        if contains(fn, 'sensors_2axis') && ...
                ~any(strcmp(fn, cellfun(@(a) a.field, arrays, 'UniformOutput', false)))
            arrays{end+1} = struct('field', fn, 'label', fn, 'xy_only', true); %#ok<AGROW>
            fprintf('  Detected additional ESG array (XY-only): %s\n', fn);
        end
    end

    if isempty(arrays)
        warning('No sensor array fields found in geometry: %s', model_name);
    end

    % Report XY-only constraint for any ESG arrays
    for a = 1:numel(arrays)
        if arrays{a}.xy_only
            fprintf('  ESG array detected (%s) — shifts constrained to XY plane (Z=0)\n', ...
                arrays{a}.field);
        end
    end
end

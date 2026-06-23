% pt_load_leadfields - Load and organise perturbed leadfields for msg_pert
%
% Loads all available leadfield files for the original and perturbed geometry
% configurations defined in config_pert. Supports any combination of BEM,
% FEM, Biot-Savart, and single-sphere forward models. Organises everything
% into a single leadfields_organised.mat file ready for pt_compute_rsq and
% the rest of the msg_pert analysis pipeline.
%
% This replaces running load_and_organise_leadfields from msg_fwd. Run it
% here after all forward models have been computed in msg_fwd on the shifted
% geometry files.
%
% USAGE:
%   pt_load_leadfields
%
% WORKFLOW:
%   1. Configure which methods are available (BEM, FEM, Biot-Savart, sphere)
%   2. Configure path for each method (defaults to forward_fields_base)
%   3. Script scans for all files matching each geometry × method combination
%   4. Organises into leadfields struct and saves leadfields_organised.mat
%
% GEOMETRY NAMES LOADED:
%   Source shifts:  sensitivity_ref_key   (original) + sensitivity_keys (18)
%   Sensor shifts:  sensor_sensitivity_ref_key (original) + sensor_sensitivity_keys (24)
%   Keys with no matching files are skipped with a warning.
%
% FILE NAMING CONVENTIONS (matching msg_fwd output):
%   BEM:   <forward_fields_base>/geometries_<geom>/leadfield_<geom>_bem_<array>.mat
%          Variable: leadfield_cord   | Scale: 1e15 (T/nAm → fT/nAm)
%   FEM:   <forward_fields_base>/geometries_<geom>/cord_leadfield_<geom>_fem_<array>.mat
%          Variable: leadfield_ft     | Scale: 1 (already fT/nAm)
%   BS:    <bslaw_path>/leadfield_<geom>_bslaw_<array>.mat
%          Variable: leadfield_bs     | Scale: 1 (already fT/nAm)
%   Sphere:<sphere_path>/leadfield_<geom>_sphere_<array>.mat
%          Variable: leadfield_sphere | Scale: 1 (already fT/nAm)
%
%   BEM and FEM files are looked up in per-geometry subfolders under the
%   respective base path. Biot-Savart and sphere files are looked up in a
%   flat folder (no subfolders) — consistent with msg_fwd/simpler_models.
%
% OUTPUT FILE:
%   <forward_fields_base>/leadfields_organised.mat  containing:
%     leadfields          — struct with one field per loaded key
%                           (e.g. leadfields.bem_sub001_source_X_p2mm_front)
%                           Each field contains .VD, .RC, .LR cell arrays
%                           [n_sensor_axes x n_sources], plus metadata fields
%     abs_max_per_source  — struct of peak absolute amplitudes per source
%     loaded_models       — cell array of all successfully loaded keys
%
% DEPENDENCIES:
%   config_pert           — paths, geometry key lists, orientation labels
%   pt_add_functions      — adds msg_fwd/functions/ (organise_leadfield) to path
%   organise_leadfield()  — from msg_fwd/functions/
%
% NOTES:
%   - Only methods with have_<method> = true are searched
%   - Missing files produce a warning and are skipped (not an error)
%   - The original (unshifted) geometry is always included under both
%     sensitivity_ref_key and sensor_sensitivity_ref_key
%   - If both source-shift and sensor-shift originals are the same geometry,
%     only one copy is loaded (keys differ so both appear in leadfields)
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
%
% This file is part of the MSG Perturbation Toolbox (msg_pert).
% Used in conjunction with msg_coreg and msg_fwd:
%   https://github.com/maikeschmidt/msg_coreg
%   https://github.com/maikeschmidt/msg_fwd

clearvars
close all
clc

config_pert;
pt_add_functions;


% =========================================================================
% USER CONFIGURATION — which forward models are available?
% =========================================================================
% Set have_<method> = true for each method you ran in msg_fwd.
% All four can be true simultaneously.

have_bem    = true;    % BEM via Helsinki BEM Framework
have_fem    = false;   % FEM via DUNEuro
have_bslaw  = false;   % Biot-Savart (infinite homogeneous space)
have_sphere = false;   % Single sphere (Sarvas analytical solution)

% Output paths for each method.
% BEM and FEM: files are in per-geometry subfolders under the base path.
% Biot-Savart and sphere: files are in a flat folder (no subfolders).
% Defaults to forward_fields_base — override if you stored outputs elsewhere.

bem_path    = forward_fields_base;   % override if BEM output is elsewhere
fem_path    = forward_fields_base;   % override if FEM output is elsewhere
bslaw_path  = forward_fields_base;   % override if BS  output is elsewhere
sphere_path = forward_fields_base;   % override if sphere output is elsewhere


% =========================================================================
% BUILD LIST OF ALL GEOMETRY NAMES TO LOAD
% =========================================================================
% Collect all unique geometry names: original + source shifts + sensor shifts.
% Each key is a full geometry stem (e.g. 'geometries_sub001_source_X_p2mm').
% Duplicates are removed (ref keys may overlap if you ran both shift types
% on the same original geometry file).

source_geom_names = [sensitivity_ref_key, sensitivity_keys];
sensor_geom_names = [sensor_sensitivity_ref_key, sensor_sensitivity_keys];
all_geom_names    = unique([source_geom_names, sensor_geom_names], 'stable');

fprintf('pt_load_leadfields\n');
fprintf('  %d unique geometry names to search\n', numel(all_geom_names));
fprintf('  Methods: BEM=%d  FEM=%d  BS=%d  Sphere=%d\n\n', ...
    have_bem, have_fem, have_bslaw, have_sphere);


% =========================================================================
% LOAD ALL LEADFIELD FILES
% =========================================================================

lf_raw    = struct();   % raw FieldTrip leadfield structs, keyed by model_key
n_loaded  = 0;

for g = 1:numel(all_geom_names)
    geom_full = all_geom_names{g};   % full geometry stem, e.g. 'geometries_sub001_...'
    geom_short = regexprep(geom_full, '^geometries[_-]?', '');

    fprintf('  [%d/%d] %s\n', g, numel(all_geom_names), geom_short);

    % ------------------------------------------------------------------
    % BEM
    % Files: <bem_path>/geometries_<geom>/leadfield_<geom_short>_bem_<array>.mat
    % Variable: leadfield_cord
    % ------------------------------------------------------------------
    if have_bem
        bem_subdir = fullfile(bem_path, ['geometries_' geom_short]);
        bem_files  = dir(fullfile(bem_subdir, ...
            ['leadfield_' geom_short '_bem_*.mat']));

        for bf = 1:numel(bem_files)
            fname = bem_files(bf).name;
            tok   = regexp(fname, ...
                ['leadfield_' geom_short '_bem_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr      = tok{1}{1};
            key      = ['bem_' geom_short '_' arr];

            tmp = load(fullfile(bem_subdir, fname), 'leadfield_cord');
            if ~isfield(tmp, 'leadfield_cord')
                warning('Variable leadfield_cord not found in: %s', fname);
                continue
            end
            lf_raw.(key)   = tmp.leadfield_cord;
            lf_raw.(key).unit_scale = 1e15;   % T/nAm → fT/nAm
            n_loaded = n_loaded + 1;
            fprintf('    BEM: %s\n', key);
        end

        if numel(bem_files) == 0
            % no warning — geometry may only have sensor or source shifts
        end
    end

    % ------------------------------------------------------------------
    % FEM
    % Files: <fem_path>/geometries_<geom>/cord_leadfield_<geom_short>_fem_<array>.mat
    % Variable: leadfield_ft
    % ------------------------------------------------------------------
    if have_fem
        fem_subdir = fullfile(fem_path, ['geometries_' geom_short]);
        fem_files  = dir(fullfile(fem_subdir, ...
            ['cord_leadfield_' geom_short '_fem_*.mat']));

        for ff = 1:numel(fem_files)
            fname = fem_files(ff).name;
            tok   = regexp(fname, ...
                ['cord_leadfield_' geom_short '_fem_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr      = tok{1}{1};
            key      = ['fem_' geom_short '_' arr];

            tmp = load(fullfile(fem_subdir, fname), 'leadfield_ft');
            if ~isfield(tmp, 'leadfield_ft')
                warning('Variable leadfield_ft not found in: %s', fname);
                continue
            end
            lf_raw.(key)   = tmp.leadfield_ft;
            lf_raw.(key).unit_scale = 1;   % already fT/nAm
            n_loaded = n_loaded + 1;
            fprintf('    FEM: %s\n', key);
        end
    end

    % ------------------------------------------------------------------
    % BIOT-SAVART
    % Files: <bslaw_path>/leadfield_<geom_short>_bslaw_<array>.mat  (flat)
    % Variable: leadfield_bs
    % ------------------------------------------------------------------
    if have_bslaw
        bs_files = dir(fullfile(bslaw_path, ...
            ['leadfield_' geom_short '_bslaw_*.mat']));

        for bf = 1:numel(bs_files)
            fname = bs_files(bf).name;
            tok   = regexp(fname, ...
                ['leadfield_' geom_short '_bslaw_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr      = tok{1}{1};
            key      = ['bslaw_' geom_short '_' arr];

            tmp = load(fullfile(bslaw_path, fname), 'leadfield_bs');
            if ~isfield(tmp, 'leadfield_bs')
                warning('Variable leadfield_bs not found in: %s', fname);
                continue
            end
            lf_raw.(key)   = tmp.leadfield_bs;
            lf_raw.(key).unit_scale = 1;   % already fT/nAm
            n_loaded = n_loaded + 1;
            fprintf('    BS:  %s\n', key);
        end
    end

    % ------------------------------------------------------------------
    % SINGLE SPHERE
    % Files: <sphere_path>/leadfield_<geom_short>_sphere_<array>.mat  (flat)
    % Variable: leadfield_sphere
    % ------------------------------------------------------------------
    if have_sphere
        sp_files = dir(fullfile(sphere_path, ...
            ['leadfield_' geom_short '_sphere_*.mat']));

        for sf = 1:numel(sp_files)
            fname = sp_files(sf).name;
            tok   = regexp(fname, ...
                ['leadfield_' geom_short '_sphere_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr      = tok{1}{1};
            key      = ['sphere_' geom_short '_' arr];

            tmp = load(fullfile(sphere_path, fname), 'leadfield_sphere');
            if ~isfield(tmp, 'leadfield_sphere')
                warning('Variable leadfield_sphere not found in: %s', fname);
                continue
            end
            lf_raw.(key)   = tmp.leadfield_sphere;
            lf_raw.(key).unit_scale = 1;   % already fT/nAm
            n_loaded = n_loaded + 1;
            fprintf('    Sp:  %s\n', key);
        end
    end
end

fprintf('\nLoaded %d leadfield files.\n\n', n_loaded);

if n_loaded == 0
    error(['No leadfield files found. Check that:\n' ...
           '  1. forward_fields_base (and bem/fem/bslaw/sphere paths) point to the\n' ...
           '     correct directories\n' ...
           '  2. The correct have_<method> flags are set to true\n' ...
           '  3. msg_fwd has been run on the perturbed geometry files']);
end


% =========================================================================
% ORGANISE INTO LEADFIELDS STRUCT
% =========================================================================
% Uses organise_leadfield() from msg_fwd/functions/ (added to path by
% pt_add_functions). Splits each leadfield matrix into VD/RC/LR cell arrays
% and computes peak absolute amplitude per source.

fprintf('Organising leadfields by source and orientation...\n');
leadfields = struct();
abs_max_per_source = struct();

raw_keys = fieldnames(lf_raw);

for k = 1:numel(raw_keys)
    key        = raw_keys{k};
    lf_struct  = lf_raw.(key);
    unit_scale = lf_struct.unit_scale;
    lf_struct  = rmfield(lf_struct, 'unit_scale');   % remove temp field

    [leadfields, abs_max_per_source] = organise_leadfield( ...
        leadfields, abs_max_per_source, lf_struct, key, unit_scale, orientation_labels);

    fprintf('  Organised: %s\n', key);
end

loaded_models = fieldnames(leadfields);
fprintf('\nOrganised %d model configurations.\n', numel(loaded_models));


% =========================================================================
% SAVE
% =========================================================================

outfile = fullfile(forward_fields_base, 'leadfields_organised.mat');
save(outfile, 'leadfields', 'abs_max_per_source', 'loaded_models', '-v7.3');
fprintf('\nSaved: %s\n', outfile);
fprintf('\nNext: run pt_compute_rsq (or run_perturbation_analysis)\n');

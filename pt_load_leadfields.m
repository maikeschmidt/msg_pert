% pt_load_leadfields - Load and organise perturbed leadfields for msg_pert
%
% Loads all available leadfield files for the original and perturbed geometry
% configurations defined in config_pert. Supports any combination of BEM,
% FEM, Biot-Savart, and single-sphere forward models. Organises everything
% into a single leadfields_organised.mat file ready for pt_compute_rsq and
% the rest of the msg_pert analysis pipeline.
%
% Each forward model method is stored under a distinct key prefix so that
% results from different methods can be compared directly:
%   bslaw_<geom>   — Biot-Savart
%   sphere_<geom>  — Single sphere
%   bem_<geom>     — BEM
%   fem_<geom>     — FEM
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
%   BEM:   <bem_path>/geometries_<geom_short>/leadfield_<geom_short>_bem_<array>.mat
%          Variable: leadfield_cord   | Scale: 1e15 (T/nAm → fT/nAm)
%   FEM:   <fem_path>/geometries_<geom_short>/cord_leadfield_<geom_short>_fem_<array>.mat
%          Variable: leadfield_ft     | Scale: 1 (already fT/nAm)
%   BS:    <bslaw_path>/leadfield_geometries_<geom_full>_bslaw_<array>.mat
%          Variable: leadfield_bs     | Scale: 1 (already fT/nAm)
%   Sphere:<sphere_path>/leadfield_geometries_<geom_full>_sphere_<array>.mat
%          Variable: leadfield_sphere | Scale: 1 (already fT/nAm)
%
%   BEM and FEM files live in per-geometry subfolders; BS and sphere files
%   are in a flat folder (no subfolders) — consistent with msg_fwd output.
%
% OUTPUT FILE:
%   <forward_fields_base>/leadfields_organised.mat  containing:
%     leadfields          — struct with one field per loaded key, e.g.
%                           leadfields.bslaw_original_source_original
%                           Each field: .VD/.RC/.LR cell arrays plus metadata
%     abs_max_per_source  — struct of peak absolute amplitudes per source
%     loaded_models       — cell array of all successfully loaded keys
%
% DEPENDENCIES:
%   config_pert           — paths, geometry key lists, orientation labels
%   pt_add_functions      — adds msg_fwd/functions/ to path
%   organise_leadfield()  — from msg_fwd/functions/
%
% NOTES:
%   - Only methods with have_<method> = true are searched
%   - Missing files produce a warning and are skipped (not an error)
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

clearvars
close all
clc

config_pert;
pt_add_functions;


% =========================================================================
% USER CONFIGURATION — which forward models are available?
% =========================================================================

have_bem      = true;    % BEM via Helsinki BEM Framework
have_fem      = false;   % FEM via DUNEuro
have_bslaw    = true;    % Biot-Savart (infinite homogeneous space)
have_sphere   = false;    % Single sphere (Sarvas analytical solution)
have_bem_cond = false;   % BEM with perturbed tissue conductivities (run_conductivity_perturbation)

% Output paths for each method.
% BEM and FEM: files are in per-geometry subfolders under the base path.
% Biot-Savart and sphere: files are in a flat folder (no subfolders).
% Defaults to forward_fields_base — override if you stored outputs elsewhere.

bem_path      = 'D:\Simulations\Pertubations\fields\bem';
fem_path      = forward_fields_base;
bslaw_path    = 'D:\Simulations\Pertubations\fields\bs_law';
sphere_path   = '';
bem_cond_path = 'D:\Simulations\Pertubations\fields\bem_cond_msg';   % SET THIS: from run_conductivity_perturbation


% =========================================================================
% BUILD LIST OF ALL GEOMETRY NAMES TO LOAD
% =========================================================================

source_geom_names = [sensitivity_ref_key, sensitivity_keys];
sensor_geom_names = [sensor_sensitivity_ref_key, sensor_sensitivity_keys];
cond_geom_names   = [cond_sensitivity_ref_key, cond_sensitivity_keys];
all_geom_names    = unique([source_geom_names, sensor_geom_names, cond_geom_names], 'stable');

fprintf('pt_load_leadfields\n');
fprintf('  %d unique geometry names to search\n', numel(all_geom_names));
fprintf('  Methods: BEM=%d  FEM=%d  BS=%d  Sphere=%d  BEM-cond=%d\n\n', ...
    have_bem, have_fem, have_bslaw, have_sphere, have_bem_cond);


% =========================================================================
% LOAD ALL LEADFIELD FILES
% =========================================================================

n_loaded  = 0;
fprintf('Loading and organising leadfields...\n');
leadfields = struct();
abs_max_per_source = struct();

for g = 1:numel(all_geom_names)
    geom_full  = all_geom_names{g};   % e.g. 'original_source_original'
    geom_short = regexprep(geom_full, '^geometries[_-]?', '');

    fprintf('  [%d/%d] %s\n', g, numel(all_geom_names), geom_full);

    % ------------------------------------------------------------------
    % BEM
    % Key: bem_<geom_full>
    % Files: <bem_path>/geometries_<geom_short>/leadfield_<geom_short>_bem_<array>.mat
    % Variable: leadfield_cord | Scale: 1e15
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
            arr = tok{1}{1};
            key = ['bem_' geom_full];

            tmp = load(fullfile(bem_subdir, fname), 'leadfield_cord');
            if ~isfield(tmp, 'leadfield_cord')
                warning('Variable leadfield_cord not found in: %s', fname);
                continue
            end
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_cord, ...
                key, 1e15, orientation_labels);
            n_loaded = n_loaded + 1;
            fprintf('    BEM: %s (%s)\n', key, arr);
        end
    end

    % ------------------------------------------------------------------
    % FEM
    % Key: fem_<geom_full>
    % Files: <fem_path>/geometries_<geom_short>/cord_leadfield_<geom_short>_fem_<array>.mat
    % Variable: leadfield_ft | Scale: 1
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
            arr = tok{1}{1};
            key = ['fem_' geom_full];

            tmp = load(fullfile(fem_subdir, fname), 'leadfield_ft');
            if ~isfield(tmp, 'leadfield_ft')
                warning('Variable leadfield_ft not found in: %s', fname);
                continue
            end
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_ft, ...
                key, 1, orientation_labels);
            n_loaded = n_loaded + 1;
            fprintf('    FEM: %s (%s)\n', key, arr);
        end
    end

    % ------------------------------------------------------------------
    % BIOT-SAVART
    % Key: bslaw_<geom_full>
    % Files: <bslaw_path>/leadfield_geometries_<geom_full>_bslaw_<array>.mat
    % Variable: leadfield_bs | Scale: 1
    % ------------------------------------------------------------------
    if have_bslaw
        bs_files = dir(fullfile(bslaw_path, ...
            ['leadfield_geometries_' geom_full '_bslaw_*.mat']));

        for bf = 1:numel(bs_files)
            fname = bs_files(bf).name;
            tok   = regexp(fname, ...
                ['leadfield_geometries_' geom_full '_bslaw_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr = tok{1}{1};
            key = ['bslaw_' geom_full];

            tmp = load(fullfile(bslaw_path, fname), 'leadfield_bs');
            if ~isfield(tmp, 'leadfield_bs')
                warning('Variable leadfield_bs not found in: %s', fname);
                continue
            end
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_bs, ...
                key, 1, orientation_labels);
            n_loaded = n_loaded + 1;
            fprintf('    BS:  %s (%s)\n', key, arr);
        end
    end

    % ------------------------------------------------------------------
    % BEM — CONDUCTIVITY PERTURBATION
    % Key: bem_cond_<geom_full>
    % Files: <bem_cond_path>/leadfield_<geom_short>_bem_cond_<array>.mat
    % Variable: leadfield_cord | Scale: 1e15
    % geom_full for cond files is always the base source-original geometry;
    % the perturbation variant is encoded in the filename suffix.
    % ------------------------------------------------------------------
    if have_bem_cond
        cond_subdir = fullfile(bem_cond_path, ['geometries_' geom_short]);
        cond_files  = dir(fullfile(cond_subdir, ...
            ['leadfield_' geom_short '_bem_cond_*.mat']));

        for cf = 1:numel(cond_files)
            fname = cond_files(cf).name;
            tok   = regexp(fname, ...
                ['leadfield_' geom_short '_bem_cond_(bundle\d+_shift\d+)_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            bund_shift = tok{1}{1};   % e.g. 'bundle1_shift1'
            arr        = tok{1}{2};   % e.g. 'front' or 'back'
            key        = ['bem_cond_' geom_full '_' bund_shift];

            tmp = load(fullfile(cond_subdir, fname), 'leadfield_cord');
            if ~isfield(tmp, 'leadfield_cord')
                warning('Variable leadfield_cord not found in: %s', fname);
                continue
            end
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_cord, ...
                key, 1e15, orientation_labels);
            n_loaded = n_loaded + 1;
            fprintf('    BEM-cond: %s (%s)\n', key, arr);
        end
    end

    % ------------------------------------------------------------------
    % SINGLE SPHERE
    % Key: sphere_<geom_full>
    % Files: <sphere_path>/leadfield_geometries_<geom_full>_sphere_<array>.mat
    % Variable: leadfield_sphere | Scale: 1
    % ------------------------------------------------------------------
    if have_sphere
        sp_files = dir(fullfile(sphere_path, ...
            ['leadfield_geometries_' geom_full '_sphere_*.mat']));

        for sf = 1:numel(sp_files)
            fname = sp_files(sf).name;
            tok   = regexp(fname, ...
                ['leadfield_geometries_' geom_full '_sphere_(.+)\.mat'], 'tokens');
            if isempty(tok); continue; end
            arr = tok{1}{1};
            key = ['sphere_' geom_full];

            tmp = load(fullfile(sphere_path, fname), 'leadfield_sphere');
            if ~isfield(tmp, 'leadfield_sphere')
                warning('Variable leadfield_sphere not found in: %s', fname);
                continue
            end
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_sphere, ...
                key, 1, orientation_labels);
            n_loaded = n_loaded + 1;
            fprintf('    Sp:  %s (%s)\n', key, arr);
        end
    end
end

fprintf('\nLoaded and organised %d leadfield files.\n\n', n_loaded);

if n_loaded == 0
    error(['No leadfield files found. Check that:\n' ...
           '  1. forward_fields_base (and bem/fem/bslaw/sphere paths) point to the\n' ...
           '     correct directories\n' ...
           '  2. The correct have_<method> flags are set to true\n' ...
           '  3. msg_fwd has been run on the perturbed geometry files']);
end

loaded_models = fieldnames(leadfields);
fprintf('Organised %d model configurations.\n', numel(loaded_models));


% =========================================================================
% SAVE
% =========================================================================

outfile = fullfile(forward_fields_base, 'leadfields_organised.mat');
save(outfile, 'leadfields', 'abs_max_per_source', 'loaded_models', '-v7.3');
fprintf('\nSaved: %s\n', outfile);
fprintf('\nNext: run pt_compute_rsq\n');
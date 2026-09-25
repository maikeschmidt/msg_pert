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
%   Source shifts:  sensitivity_ref_key   (original) + sensitivity_keys (24)
%   Sensor shifts:  sensor_sensitivity_ref_key (original) + sensor_sensitivity_keys (24)
%   Keys with no matching files are skipped with a warning.
%
% FILE NAMING CONVENTIONS (matching msg_fwd output):
%   BEM:   <bem_path>/geometries_<geom_short>/leadfield_<geom_short>_bem_<array>.mat
%          Variable: leadfield_cord
%   FEM:   <fem_path>/geometries_<geom_short>/cord_leadfield_<geom_short>_fem_<array>.mat
%          Variable: leadfield_ft
%   BS:    <bslaw_path>/leadfield_geometries_<geom_full>_bslaw_<array>.mat
%          Variable: leadfield_bs
%   Sphere:<sphere_path>/leadfield_geometries_<geom_full>_sphere_<array>.mat
%          Variable: leadfield_sphere
%
%   BEM and FEM files live in per-geometry subfolders; BS and sphere files
%   are in a flat folder (no subfolders) — consistent with msg_fwd output.
%
% OUTPUT FILE:
%   <forward_fields_base>/leadfield_scale_log.csv   factor applied to each file
%   <forward_fields_base>/leadfields_organised.mat  containing:
%     leadfields          — struct with one field per loaded key, e.g.
%                           leadfields.bslaw_original_source_original
%                           Each field: .VD/.RC/.LR cell arrays plus metadata
%     abs_max_per_source  — struct of peak absolute amplitudes per source
%     loaded_models       — cell array of all successfully loaded keys
%     scale_log           — table: key, file, scale, reason
%     analysis_array_loaded — the array these lead fields belong to
%
% DEPENDENCIES:
%   config_pert           — paths, geometry key lists, orientation labels
%   pt_add_functions      — adds msg_fwd/functions/ to path
%   organise_leadfield()  — from msg_fwd/functions/
%
% ARRAY AND UNITS:
%   Only files for analysis_array (config_pert) are loaded. Each file is
%   scaled on its own (unit_scale_mode = 'auto', see pt_unit_scale), because
%   raw BEM output, Biot-Savart output and conductivity-perturbation files do
%   not share one unit convention. Every choice is written to
%   leadfield_scale_log.csv, and any perturbed lead field more than 10x
%   larger or smaller than its reference is flagged as a units problem.
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
% CONFIGURATION — all from config_pert (per active modality)
% =========================================================================
% config_pert defines both MSG and ESG side by side and activates the one
% selected by pt_modality (the master loop sets it; standalone defaults to MSG).
% The variables used below all come from there:
%   have_bem / have_fem / have_bslaw / have_sphere / have_bem_cond
%   bem_path / fem_path / bslaw_path / sphere_path / bem_cond_path
%   sensor_n_axes / sensor_is_meg / bem_unit_scale / forward_fields_base
%
% To run a single modality standalone, call pt_modality('set','esg') before
% pt_load_leadfields, or edit active_modality's fallback in config_pert.
fprintf('  Modality: %s  (%d axes, is_meg=%d)\n', ...
    active_modality, sensor_n_axes, sensor_is_meg);


% =========================================================================
% BUILD LIST OF ALL GEOMETRY NAMES TO LOAD
% =========================================================================

% BEM-cond files are all stored under the original source geometry —
% they are NOT separate geometry variants, so they are handled in their
% own section below (not in the main geometry loop).
source_geom_names = [sensitivity_ref_key, sensitivity_keys];
sensor_geom_names = [sensor_sensitivity_ref_key, sensor_sensitivity_keys];
all_geom_names    = unique([source_geom_names, sensor_geom_names], 'stable');

fprintf('pt_load_leadfields\n');
fprintf('  %d unique geometry names to search\n', numel(all_geom_names));
fprintf('  Methods: BEM=%d  FEM=%d  BS=%d  Sphere=%d  BEM-cond=%d\n\n', ...
    have_bem, have_fem, have_bslaw, have_sphere, have_bem_cond);


% =========================================================================
% LOAD ALL LEADFIELD FILES
% =========================================================================

n_loaded  = 0;

% Unit scale per file. In 'auto' mode each file is scaled from its own
% magnitude (pt_unit_scale), because the writers feeding this loader do not
% share one convention — see config_pert. Every decision is logged and saved.
if strcmp(unit_scale_mode, 'auto')
    scale_for = @(lf, method) pt_unit_scale(lf, sensor_is_meg, method);
else
    scale_for = @(lf, method) fixed_scale(method, bem_unit_scale);
end
scale_log = cell(0, 4);   % {key, file, scale, reason}
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
            ['leadfield_' geom_short '_bem_' analysis_array '.mat']));

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
            [sc, last_why] = scale_for(tmp.leadfield_cord, 'bem');
            scale_log(end+1, :) = {key, fname, sc, last_why}; %#ok<SAGROW>
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_cord, ...
                key, sc, orientation_labels, ...
                sensor_n_axes, sensor_is_meg);
            n_loaded = n_loaded + 1;
            fprintf('    BEM: %s (%s) %s\n', key, arr, last_why);
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
            ['cord_leadfield_' geom_short '_fem_' analysis_array '.mat']));

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
            [sc, last_why] = scale_for(tmp.leadfield_ft, 'fem');
            scale_log(end+1, :) = {key, fname, sc, last_why}; %#ok<SAGROW>
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_ft, ...
                key, sc, orientation_labels, ...
                sensor_n_axes, sensor_is_meg);
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
            ['leadfield_geometries_' geom_full '_bslaw_' analysis_array '.mat']));

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
            [sc, last_why] = scale_for(tmp.leadfield_bs, 'bslaw');
            scale_log(end+1, :) = {key, fname, sc, last_why}; %#ok<SAGROW>
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_bs, ...
                key, sc, orientation_labels, ...
                sensor_n_axes, sensor_is_meg);
            n_loaded = n_loaded + 1;
            fprintf('    BS:  %s (%s)\n', key, arr);
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
            ['leadfield_geometries_' geom_full '_sphere_' analysis_array '.mat']));

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
            [sc, last_why] = scale_for(tmp.leadfield_sphere, 'sphere');
            scale_log(end+1, :) = {key, fname, sc, last_why}; %#ok<SAGROW>
            [leadfields, abs_max_per_source] = organise_leadfield( ...
                leadfields, abs_max_per_source, tmp.leadfield_sphere, ...
                key, sc, orientation_labels, ...
                sensor_n_axes, sensor_is_meg);
            n_loaded = n_loaded + 1;
            fprintf('    Sp:  %s (%s)\n', key, arr);
        end
    end
end

% =========================================================================
% BEM — CONDUCTIVITY PERTURBATION (standalone section)
% =========================================================================
% All cond leadfields live under a single geometry folder (the unshifted
% original). Keys: bem_cond_<ref>_bundle<B>_shift<S>
% Files: <bem_cond_path>/geometries_<ref_short>/leadfield_<ref_short>_bem_cond_bundle<B>_shift<S>_<array>.mat

if have_bem_cond
    fprintf('\n  --- BEM Conductivity Perturbation ---\n');
    ref_short    = regexprep(cond_sensitivity_ref_key, '^geometries[_-]?', '');
    cond_subdir  = fullfile(bem_cond_path, ['geometries_' ref_short]);
    cond_files   = dir(fullfile(cond_subdir, ...
        ['leadfield_' ref_short '_bem_cond_*_' analysis_array '.mat']));

    if isempty(cond_files)
        warning('No BEM-cond files found in: %s', cond_subdir);
    end

    for cf = 1:numel(cond_files)
        fname = cond_files(cf).name;
        tok   = regexp(fname, ...
            ['leadfield_' ref_short '_bem_cond_(bundle\d+_shift\d+)_(.+)\.mat'], 'tokens');
        if isempty(tok); continue; end
        bund_shift = tok{1}{1};   % e.g. 'bundle1_shift1'
        arr        = tok{1}{2};   % e.g. 'front' or 'back'
        key        = ['bem_cond_' cond_sensitivity_ref_key '_' bund_shift];

        tmp = load(fullfile(cond_subdir, fname), 'leadfield_cord');
        if ~isfield(tmp, 'leadfield_cord')
            warning('Variable leadfield_cord not found in: %s', fname);
            continue
        end
        [sc, last_why] = scale_for(tmp.leadfield_cord, 'bem_cond');
        scale_log(end+1, :) = {key, fname, sc, last_why}; %#ok<SAGROW>
        [leadfields, abs_max_per_source] = organise_leadfield( ...
            leadfields, abs_max_per_source, tmp.leadfield_cord, ...
            key, sc, orientation_labels, ...
            sensor_n_axes, sensor_is_meg);
        n_loaded = n_loaded + 1;
        fprintf('    BEM-cond: %s (%s) %s\n', key, arr, last_why);
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
% MAGNITUDE CHECK AGAINST EACH REFERENCE
% =========================================================================
% A geometric shift or a conductivity change of the size simulated here
% changes the overall field magnitude by tens of percent at most. A factor
% of 10 between a perturbed leadfield and its reference is a units problem,
% and it would make every RE and lnMAG value meaningless.

n_flag = 0;
for k = 1:numel(loaded_models)
    key = loaded_models{k};
    rk  = regexprep(key, '^bem_cond_(.*)_bundle\d+_shift\d+$', 'bem_$1');
    rk  = regexprep(rk, '_(source|sensor)_bundle\d+_shift\d+$', '_$1_original');
    if strcmp(rk, key) || ~isfield(leadfields, rk), continue; end
    r = median_norm(leadfields.(key)) / median_norm(leadfields.(rk));
    if r > 10 || r < 0.1
        warning('pt_load_leadfields:magnitude', ...
            ['%s is %.3g x its reference %s. That is a units mismatch, not ' ...
             'a perturbation effect. Check leadfield_scale_log.csv.'], key, r, rk);
        n_flag = n_flag + 1;
    end
end
if n_flag == 0
    fprintf('Magnitude check: every perturbed leadfield is within 10x of its reference.\n');
end


% =========================================================================
% SAVE
% =========================================================================

outfile = fullfile(forward_fields_base, 'leadfields_organised.mat');
scale_log = cell2table(scale_log, 'VariableNames', {'key', 'file', 'scale', 'reason'});
analysis_array_loaded = analysis_array; %#ok<NASGU>
save(outfile, 'leadfields', 'abs_max_per_source', 'loaded_models', ...
    'scale_log', 'analysis_array_loaded', '-v7.3');
writetable(scale_log, fullfile(forward_fields_base, 'leadfield_scale_log.csv'));
fprintf('\nSaved: %s\n', outfile);
fprintf('\nNext: run pt_compute_rsq\n');


% ---- Local functions ----

function [s, why] = fixed_scale(method, bem_unit_scale)
    if startsWith(method, 'bem')
        s = bem_unit_scale;
    else
        s = 1;
    end
    why = sprintf('fixed x%g', s);
end

function n = median_norm(E)
    ax = E.n_sensor_axes;
    v  = zeros(1, E.n_sources);
    for s = 1:E.n_sources
        v(s) = norm(E.VD{ax, s});
    end
    n = median(v);
end

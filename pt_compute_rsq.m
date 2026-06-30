% pt_compute_rsq - Compute per-source r² for perturbation analysis
%
% Loads all leadfields for shifted models (from msg_fwd output), computes
% the squared Pearson correlation (r²) between each shifted model and the
% original reference at every source position, for each forward model method
% defined in config_pert (fwd_methods).
%
% Supports two analysis modes, run independently or together:
%   SOURCE mode — cord/bone shifted ±2, ±4, ±6 mm per axis (18 configs)
%   SENSOR mode — sensor array shifted in 3 bundles × 8 realisations (24 configs)
%
% r² is computed independently per method so results can be compared across
% forward models in pt_plot_curves.
%
% USAGE:
%   pt_compute_rsq
%
% OUTPUTS (saved to <forward_fields_base>):
%   pert_source_rsq.mat  — source mode results (per method)
%   pert_sensor_rsq.mat  — sensor mode results (per method)
%
% DEPENDENCIES:
%   config_pert, leadfields_organised.mat
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_pert
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
% Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

% SETTINGS

run_source = true;    % SET THIS: compute source shift r²
run_sensor = true;    % SET THIS: compute sensor shift r²

% INITIALISE

config_pert;
pt_add_functions;

load(fullfile(forward_fields_base, 'leadfields_organised.mat'), ...
    'leadfields', 'abs_max_per_source', 'loaded_models');

n_methods = numel(fwd_methods);


%% SOURCE MODE

if run_source

    fprintf('\n SOURCE SHIFT r² \n');

    % Determine geometry dimensions from the first available method
    ref_key_found = '';
    for m = 1:n_methods
        candidate = [fwd_methods{m} '_' sensitivity_ref_key];
        if isfield(leadfields, candidate)
            ref_key_found = candidate;
            break
        end
    end

    if isempty(ref_key_found)
        warning('Source reference not found for any method — skipping source mode.');
    else
        % Use first found method to get shared geometry info
        [valid_keys_geom, valid_labels, valid_idx] = validate_keys( ...
            sensitivity_keys, sensitivity_labels, leadfields, fwd_methods{1});

        valid_source_bundle_idx = source_bundle_idx(valid_idx);
        valid_source_shift_idx  = source_shift_idx(valid_idx);

        [n_sources, n_axes, src_range, n_src_plot, distances, marker_idx] = ...
            get_dimensions(leadfields, ref_key_found, orientation_labels, src_spacing_mm);

        % Compute r² per method
        rsq_by_method  = struct();
        loaded_methods = {};
        for m = 1:n_methods
            method  = fwd_methods{m};
            ref_key = [method '_' sensitivity_ref_key];
            if ~isfield(leadfields, ref_key)
                warning('Source reference not found for method %s (%s) — skipping.', ...
                    method, ref_key);
                continue
            end
            method_shift_keys = cellfun(@(k) [method '_' k], ...
                sensitivity_keys(valid_idx), 'UniformOutput', false);
            [method_keys, ~, ~] = validate_keys_direct(method_shift_keys, leadfields);

            min_sensors = get_min_sensors(leadfields, ref_key, method_keys, orientation_labels);

            fprintf('  [%s] Reference: %s  |  %d valid shifts\n', ...
                method, ref_key, numel(method_keys));
            rsq_by_method.(method) = compute_rsq(leadfields, ref_key, ...
                method_keys, orientation_labels, src_range, n_axes, min_sensors);
            loaded_methods{end+1} = method; %#ok<AGROW>
        end

        outfile = fullfile(forward_fields_base, 'pert_source_rsq.mat');
        save(outfile, ...
            'rsq_by_method', 'loaded_methods', ...
            'valid_keys_geom', 'valid_labels', 'valid_idx', ...
            'valid_source_bundle_idx', 'valid_source_shift_idx', ...
            'n_sources', 'n_axes', 'src_range', 'n_src_plot', ...
            'distances', 'marker_idx', '-v7.3');
        fprintf('  Saved: %s\n', outfile);
    end
end


%% SENSOR MODE

if run_sensor

    fprintf('\n SENSOR SHIFT r² \n');

    ref_key_found = '';
    for m = 1:n_methods
        candidate = [fwd_methods{m} '_' sensor_sensitivity_ref_key];
        if isfield(leadfields, candidate)
            ref_key_found = candidate;
            break
        end
    end

    if isempty(ref_key_found)
        warning('Sensor reference not found for any method — skipping sensor mode.');
    else
        [valid_keys_geom, valid_labels, valid_idx] = validate_keys( ...
            sensor_sensitivity_keys, sensor_sensitivity_labels, leadfields, fwd_methods{1});

        valid_bundle_idx = sensor_sensitivity_bundle_idx(valid_idx);
        valid_shift_idx  = sensor_sensitivity_shift_idx(valid_idx);

        [n_sources, n_axes, src_range, n_src_plot, distances, marker_idx] = ...
            get_dimensions(leadfields, ref_key_found, orientation_labels, src_spacing_mm);

        rsq_by_method  = struct();
        loaded_methods = {};
        for m = 1:n_methods
            method  = fwd_methods{m};
            ref_key = [method '_' sensor_sensitivity_ref_key];
            if ~isfield(leadfields, ref_key)
                warning('Sensor reference not found for method %s (%s) — skipping.', ...
                    method, ref_key);
                continue
            end
            method_shift_keys = cellfun(@(k) [method '_' k], ...
                sensor_sensitivity_keys(valid_idx), 'UniformOutput', false);
            [method_keys, ~, ~] = validate_keys_direct(method_shift_keys, leadfields);

            min_sensors = get_min_sensors(leadfields, ref_key, method_keys, orientation_labels);

            fprintf('  [%s] Reference: %s  |  %d valid shifts\n', ...
                method, ref_key, numel(method_keys));
            rsq_by_method.(method) = compute_rsq(leadfields, ref_key, ...
                method_keys, orientation_labels, src_range, n_axes, min_sensors);
            loaded_methods{end+1} = method; %#ok<AGROW>
        end

        outfile = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
        save(outfile, ...
            'rsq_by_method', 'loaded_methods', ...
            'valid_keys_geom', 'valid_labels', 'valid_idx', ...
            'valid_bundle_idx', 'valid_shift_idx', ...
            'n_sources', 'n_axes', 'src_range', 'n_src_plot', ...
            'distances', 'marker_idx', '-v7.3');
        fprintf('  Saved: %s\n', outfile);
    end
end

fprintf('\n pt_compute_rsq complete \n');


% ---- Local functions ----

function rsq_store = compute_rsq(leadfields, ref_key, valid_keys, ...
    orientation_labels, src_range, n_axes, min_sensors)
    n_valid    = numel(valid_keys);
    n_src_plot = numel(src_range);
    rsq_store  = struct();
    for ori_idx = 1:numel(orientation_labels)
        ori_label = orientation_labels{ori_idx};
        rsq_mat   = nan(n_valid, n_src_plot, n_axes);
        for ax = 1:n_axes
            for i = 1:n_valid
                for si = 1:n_src_plot
                    src_idx = src_range(si);
                    vecA    = leadfields.(ref_key).(ori_label){ax, src_idx}(1:min_sensors);
                    vecB    = leadfields.(valid_keys{i}).(ori_label){ax, src_idx}(1:min_sensors);
                    tmp     = corrcoef(vecA, vecB);
                    rsq_mat(i, si, ax) = tmp(1, 2)^2;
                end
            end
        end
        rsq_store.(ori_label) = rsq_mat;
        fprintf('    %s: done\n', ori_label);
    end
end

function [valid_keys, valid_labels, valid_idx] = validate_keys(all_geom_keys, all_labels, leadfields, method)
    n_total   = numel(all_geom_keys);
    valid_idx = false(1, n_total);
    for i = 1:n_total
        key = [method '_' all_geom_keys{i}];
        valid_idx(i) = isfield(leadfields, key);
    end
    valid_keys   = all_geom_keys(valid_idx);
    valid_labels = all_labels(valid_idx);
    if numel(valid_keys) == 0
        warning('No valid shift models found for method %s.', method);
    end
    fprintf('  [%s] Valid geometry keys: %d of %d\n', method, numel(valid_keys), n_total);
end

function [valid_keys, dummy1, dummy2] = validate_keys_direct(full_keys, leadfields)
    dummy1 = {}; dummy2 = {};
    n = numel(full_keys);
    mask = false(1, n);
    for i = 1:n
        mask(i) = isfield(leadfields, full_keys{i});
        if ~mask(i)
            warning('Model not found, skipping: %s', full_keys{i});
        end
    end
    valid_keys = full_keys(mask);
end

function [n_sources, n_axes, src_range, n_src_plot, distances, marker_idx] = ...
        get_dimensions(leadfields, ref_key, orientation_labels, src_spacing_mm)
    n_sources  = leadfields.(ref_key).n_sources;
    n_axes     = leadfields.(ref_key).n_sensor_axes;
    src_range  = 2:(n_sources - 1);
    n_src_plot = numel(src_range);
    distances  = src_range * src_spacing_mm;
    marker_idx = 1:5:n_src_plot;
end

function min_sensors = get_min_sensors(leadfields, ref_key, valid_keys, orientation_labels)
    min_sensors = numel(leadfields.(ref_key).(orientation_labels{1}){1, 1});
    for i = 1:numel(valid_keys)
        min_sensors = min(min_sensors, ...
            numel(leadfields.(valid_keys{i}).(orientation_labels{1}){1, 1}));
    end
end
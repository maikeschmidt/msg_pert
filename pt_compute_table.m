% pt_compute_table - Write summary tables for perturbation analysis
%
% Loads pre-computed r² data and writes formatted summary tables as both
% .txt (human-readable) and .csv (for Excel/R import).
%
% Computes for each shifted model:
%   Median r²          — median across source positions
%   Min r²             — worst-case across source positions
%   First drop < 0.99  — cord position where r² first falls below 0.99
%   First drop < 0.95  — cord position where r² first falls below 0.95
%
% USAGE:
%   pt_compute_table
%
% DEPENDENCIES:
%   config_pert, pert_source_rsq.mat, pert_sensor_rsq.mat
%
% OUTPUTS:
%   <save_base_dir>/perturbation_analysis/source/source_rsq_table.txt/.csv
%   <save_base_dir>/perturbation_analysis/sensor/sensor_rsq_table.txt/.csv
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

run_source = true;
run_sensor = true;

config_pert;
pt_add_functions;

%% SOURCE TABLE

if run_source
    fprintf('SOURCE PERTURBATION TABLE\n');
    src_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
    if ~isfile(src_file)
        error('Source r² file not found: %s\nRun pt_compute_rsq first.', src_file);
    end
    load(src_file);

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'source');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    groups       = source_bundle_display;
    group_header = 'Error bundle';
    group_fn     = @(i) source_bundle_display{valid_source_bundle_idx(i)};

    header_lines = {
        '=== SOURCE POSITION PERTURBATION — SUMMARY TABLE ===', ...
        'Cord mesh + bone mesh + source positions shifted ±2, ±4, ±6 mm per axis', ...
    };

    for m = 1:numel(loaded_methods)
        method    = loaded_methods{m};
        rsq_store = rsq_by_method.(method);
        prefix    = sprintf('source_%s', method);
        write_table(rsq_store, valid_labels, orientation_labels, orientation_display, ...
            src_range, src_spacing_mm, n_axes, groups, group_header, group_fn, ...
            sensitivity_ref_key, save_dir, prefix, header_lines);
    end
    fprintf('Source table complete.\n\n');
end

%% SENSOR TABLE

if run_sensor
    fprintf('SENSOR PERTURBATION TABLE\n');
    sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
    if ~isfile(sen_file)
        error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
    end
    load(sen_file);

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    groups       = sensor_bundle_display;
    group_header = 'Error bundle';
    group_fn     = @(i) sensor_bundle_display{valid_bundle_idx(i)};

    header_lines = {
        '=== SENSOR ARRAY PERTURBATION — SUMMARY TABLE ===', ...
        'Sensor array positions randomly displaced in three error bundles', ...
        '  Bundle 1 — small  (~2mm):  U(1,3)  mm per axis', ...
        '  Bundle 2 — medium (~5mm):  U(3,7)  mm per axis', ...
        '  Bundle 3 — large  (~10mm): U(7,13) mm per axis', ...
    };

    for m = 1:numel(loaded_methods)
        method    = loaded_methods{m};
        rsq_store = rsq_by_method.(method);
        prefix    = sprintf('sensor_%s', method);
        write_table(rsq_store, valid_labels, orientation_labels, orientation_display, ...
            src_range, src_spacing_mm, n_axes, groups, group_header, group_fn, ...
            sensor_sensitivity_ref_key, save_dir, prefix, header_lines);
    end
    fprintf('Sensor table complete.\n\n');
end

fprintf('pt_compute_table complete.\n');


% ---- Local function ----

function write_table(rsq_store, valid_labels, orientation_labels, orientation_display, ...
    src_range, src_spacing_mm, n_axes, groups, group_header, group_fn, ...
    ref_key, save_dir, file_prefix, header_lines)

    n_valid       = numel(valid_labels);
    row_model     = {};  row_group   = {};  row_sens_ax  = {};  row_ori      = {};
    row_med_rsq   = [];  row_min_rsq = [];  row_min_mm   = [];
    row_drop99_mm = [];  row_drop95_mm = [];

    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};
            for i = 1:n_valid
                rsq_vec = squeeze(rsq_store.(ori_label)(i, :, sens_ax));
                med_rsq = median(rsq_vec, 'omitnan');
                [min_rsq, min_idx] = min(rsq_vec);
                min_mm  = src_range(min_idx) * src_spacing_mm;

                b99 = find(rsq_vec < 0.99, 1, 'first');
                b95 = find(rsq_vec < 0.95, 1, 'first');
                drop99_mm = NaN; if ~isempty(b99), drop99_mm = src_range(b99) * src_spacing_mm; end
                drop95_mm = NaN; if ~isempty(b95), drop95_mm = src_range(b95) * src_spacing_mm; end

                row_model{end+1}     = valid_labels{i};
                row_group{end+1}     = group_fn(i);
                row_sens_ax{end+1}   = sprintf('Axis %d', sens_ax);
                row_ori{end+1}       = orientation_display{ori_idx};
                row_med_rsq(end+1)   = med_rsq;
                row_min_rsq(end+1)   = min_rsq;
                row_min_mm(end+1)    = min_mm;
                row_drop99_mm(end+1) = drop99_mm;
                row_drop95_mm(end+1) = drop95_mm;
            end
        end
    end

    T = table(row_model', row_group', row_sens_ax', row_ori', ...
        round(row_med_rsq', 4), round(row_min_rsq', 4), ...
        row_min_mm', row_drop99_mm', row_drop95_mm', ...
        'VariableNames', {'ShiftModel','Group','SensorAxis','LeadfieldOrientation', ...
            'Median_Rsq','Min_Rsq','Min_Rsq_Position_mm', ...
            'First_Drop_Below_0p99_mm','First_Drop_Below_0p95_mm'});

    csv_path = fullfile(save_dir, [file_prefix '_rsq_table.csv']);
    writetable(T, csv_path);
    fprintf('  Saved CSV: %s\n', csv_path);

    txt_path = fullfile(save_dir, [file_prefix '_rsq_table.txt']);
    fid = fopen(txt_path, 'w');
    for k = 1:numel(header_lines); fprintf(fid, '%s\n', header_lines{k}); end
    fprintf(fid, 'Generated : %s\n', datestr(now));
    fprintf(fid, 'Reference : %s\n\n', ref_key);
    fprintf(fid, 'r²        = (Pearson r)^2 per source position\n');
    fprintf(fid, 'Drop<0.99 : cord position where r² first falls below 0.99\n');
    fprintf(fid, 'Drop<0.95 : cord position where r² first falls below 0.95\n');
    fprintf(fid, 'never     : r² remained above threshold throughout\n\n');
    fprintf(fid, 'SOURCE SPACING: %d mm\n', src_spacing_mm);
    fprintf(fid, 'EDGE SOURCES  : first and last excluded\n\n');

    divider = repmat('=', 1, 110);
    subdiv  = repmat('-', 1, 80);

    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            fprintf(fid, '%s\n', divider);
            fprintf(fid, 'SENSOR AXIS %d  |  ORIENTATION: %s\n', ...
                sens_ax, orientation_display{ori_idx});
            fprintf(fid, '%s\n', divider);

            for g = 1:numel(groups)
                fprintf(fid, '\n  %s: %s\n', group_header, groups{g});
                fprintf(fid, '  %s\n', subdiv);
                fprintf(fid, '  %-18s  %10s  %10s  %14s  %15s  %15s\n', ...
                    'Model', 'Median r²', 'Min r²', 'Min pos (mm)', ...
                    'Drop<0.99 (mm)', 'Drop<0.95 (mm)');
                fprintf(fid, '  %s\n', subdiv);

                mask  = strcmp(T.SensorAxis, sprintf('Axis %d', sens_ax)) & ...
                        strcmp(T.LeadfieldOrientation, orientation_display{ori_idx}) & ...
                        strcmp(T.Group, groups{g});
                T_sub = T(mask, :);

                if height(T_sub) == 0
                    fprintf(fid, '  [no data]\n'); continue;
                end
                for r = 1:height(T_sub)
                    s99 = 'never'; if ~isnan(T_sub.First_Drop_Below_0p99_mm(r)), s99 = sprintf('%d mm', T_sub.First_Drop_Below_0p99_mm(r)); end
                    s95 = 'never'; if ~isnan(T_sub.First_Drop_Below_0p95_mm(r)), s95 = sprintf('%d mm', T_sub.First_Drop_Below_0p95_mm(r)); end
                    fprintf(fid, '  %-18s  %10.4f  %10.4f  %14d  %15s  %15s\n', ...
                        T_sub.ShiftModel{r}, T_sub.Median_Rsq(r), T_sub.Min_Rsq(r), ...
                        T_sub.Min_Rsq_Position_mm(r), s99, s95);
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, '\n');
        end
    end

    fprintf(fid, '%s\nEND OF TABLE\n', divider);
    fclose(fid);
    fprintf('  Saved TXT: %s\n', txt_path);
end

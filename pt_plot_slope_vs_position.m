% pt_plot_slope_vs_position - Slope of r² vs perturbation along the spinal cord
%
% Reads trend tables produced by pt_plot_displacement and plots the linear-fit
% slope (dr²/d_displacement or dr²/d_%_change) as a function of source position
% over the full spinal cord length.
%
% One figure per sensor axis; three sub-panels (one per dipole orientation).
% Each forward model drawn as a separate line.
%
% Handles sensor shifts, source shifts, and conductivity perturbation.
%
% USAGE:
%   pt_plot_slope_vs_position
%
% DEPENDENCIES:
%   config_pert
%   sensor_disp_trend_table.tsv (from pt_plot_displacement)
%   source_disp_trend_table.tsv (from pt_plot_displacement)
%   cond_disp_trend_table.tsv   (from pt_plot_displacement)
%
% OUTPUTS (saved to <save_base_dir>/perturbation_analysis/<mode>/):
%   <mode>_slope_vs_position_sensorax<N>.png/.fig
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

run_sensor = true;
run_source = true;
run_cond   = true;

% INITIALISE

config_pert;
pt_add_functions;

n_ori = numel(orientation_labels);


%% =========================================================================
%% SENSOR
%% =========================================================================

if run_sensor
    tbl_file = fullfile(save_base_dir, 'perturbation_analysis', 'sensor', ...
        'sensor_disp_trend_table.tsv');
    rsq_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
    y_label  = 'Slope  (r² mm⁻¹)';

    plot_slope_mode(tbl_file, rsq_file, save_dir, 'sensor', ...
        'Rate of r² change per unit sensor displacement', y_label, ...
        fwd_method_labels, fwd_method_colors, orientation_labels, orientation_display, n_ori);
end


%% =========================================================================
%% SOURCE
%% =========================================================================

if run_source
    tbl_file = fullfile(save_base_dir, 'perturbation_analysis', 'source', ...
        'source_disp_trend_table.tsv');
    rsq_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'source');
    y_label  = 'Slope  (r² mm⁻¹)';

    plot_slope_mode(tbl_file, rsq_file, save_dir, 'source', ...
        'Rate of r² change per unit source displacement', y_label, ...
        fwd_method_labels, fwd_method_colors, orientation_labels, orientation_display, n_ori);
end


%% =========================================================================
%% CONDUCTIVITY
%% =========================================================================

if run_cond
    tbl_file = fullfile(save_base_dir, 'perturbation_analysis', 'cond', ...
        'cond_disp_trend_table.tsv');
    rsq_file = fullfile(forward_fields_base, 'pert_cond_rsq.mat');
    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'cond');
    y_label  = 'Slope  (r² per % change)';

    plot_slope_mode(tbl_file, rsq_file, save_dir, 'cond', ...
        'Rate of r² change per unit conductivity increase', y_label, ...
        fwd_method_labels, fwd_method_colors, orientation_labels, orientation_display, n_ori);
end

fprintf('\npt_plot_slope_vs_position complete.\n');


%% =========================================================================
%% LOCAL FUNCTION
%% =========================================================================

function plot_slope_mode(tbl_file, rsq_file, save_dir, mode_name, ...
    fig_title, y_label, fwd_method_labels, fwd_method_colors, ...
    orientation_labels, orientation_display, n_ori)

    fprintf('\n=== SLOPE VS POSITION: %s ===\n', upper(mode_name));

    if ~isfile(tbl_file)
        warning('Trend table not found: %s\nRun pt_plot_displacement first.', tbl_file);
        return
    end
    if ~isfile(rsq_file)
        warning('r² file not found: %s\nRun pt_compute_rsq first.', rsq_file);
        return
    end

    T      = readtable(tbl_file, 'Delimiter', '\t', 'FileType', 'text');
    tmp    = load(rsq_file, 'n_axes');
    n_axes = tmp.n_axes;
    clear tmp

    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    fprintf('  Loaded %d rows from %s\n', height(T), tbl_file);

    % Methods present in the table — preserve config order
    methods_in_table = unique(T.Method);
    ordered_methods  = {};
    ordered_colors   = [];
    for k = 1:numel(fwd_method_labels)
        lbl = fwd_method_labels{k};
        if any(strcmp(methods_in_table, lbl))
            ordered_methods{end+1} = lbl; %#ok<AGROW>
            ordered_colors(end+1,:) = fwd_method_colors(k,:); %#ok<AGROW>
        end
    end
    for k = 1:numel(methods_in_table)
        if ~any(strcmp(ordered_methods, methods_in_table{k}))
            ordered_methods{end+1} = methods_in_table{k}; %#ok<AGROW>
            ordered_colors(end+1,:) = [0.5 0.5 0.5]; %#ok<AGROW>
        end
    end
    n_methods = numel(ordered_methods);

    for sens_ax = 1:n_axes
        fig = figure('Color', 'w', 'Position', [100, 100, 1400, 420]);
        tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
        title(tl, sprintf('%s — Sensor axis %d of %d', fig_title, sens_ax, n_axes), ...
            'FontSize', 13, 'FontWeight', 'bold');
        xlabel(tl, 'Distance along the spinal cord (mm)', 'FontSize', 12);
        ylabel(tl, y_label, 'FontSize', 12);

        leg_h   = gobjects(n_methods, 1);
        leg_lbl = ordered_methods;

        for ori_idx = 1:n_ori
            ori_name = orientation_display{ori_idx};
            ax_panel = nexttile(tl, ori_idx);
            hold(ax_panel, 'on');

            for m_idx = 1:n_methods
                mlabel = ordered_methods{m_idx};
                col    = ordered_colors(m_idx,:);

                row_mask = strcmp(T.Method, mlabel) & ...
                           strcmp(T.Orientation, ori_name) & ...
                           T.SensorAxis == sens_ax;
                sub = T(row_mask, :);

                if isempty(sub)
                    leg_h(m_idx) = plot(ax_panel, NaN, NaN, '-', ...
                        'Color', col, 'LineWidth', 2.0);
                    continue
                end

                [src_sorted, si] = sort(sub.SourcePosition_mm);
                slope_sorted     = sub.Slope(si);

                h = plot(ax_panel, src_sorted, slope_sorted, '-', ...
                    'Color', col, 'LineWidth', 2.0);
                leg_h(m_idx) = h;
            end

            yline(ax_panel, 0, '--k', 'LineWidth', 1.0, 'Alpha', 0.4);

            title(ax_panel, ori_name, 'FontSize', 12, 'FontWeight', 'bold');
            if ori_idx == 1
                ylabel(ax_panel, y_label, 'FontSize', 11);
            end
            xlabel(ax_panel, 'Distance along the spinal cord (mm)', 'FontSize', 11);
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
            hold(ax_panel, 'off');
        end

        ax_last = nexttile(tl, n_ori);
        lgd = legend(ax_last, leg_h, leg_lbl, ...
            'Location', 'eastoutside', 'FontSize', 11);
        lgd.Box = 'off';
        title(lgd, 'Forward model');

        fname = sprintf('%s_slope_vs_position_sensorax%d', mode_name, sens_ax);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end

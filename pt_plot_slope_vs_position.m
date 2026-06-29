% pt_plot_slope_vs_position - Slope of r² vs displacement along spinal cord
%
% Reads sensor_disp_trend_table.tsv (produced by pt_plot_displacement) and
% plots the linear-fit slope (dr²/d_displacement) as a function of source
% position along the spinal cord.
%
% One figure per sensor axis; three sub-panels per figure (one per dipole
% orientation). Each forward model is drawn as a separate line.
%
% USAGE:
%   pt_plot_slope_vs_position
%
% DEPENDENCIES:
%   config_pert, sensor_disp_trend_table.tsv (from pt_plot_displacement)
%
% OUTPUTS (saved to <save_base_dir>/perturbation_analysis/sensor/):
%   slope_vs_position_sensorax<N>.png/.fig
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

% INITIALISE

config_pert;
pt_add_functions;

n_ori = numel(orientation_labels);

save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
tbl_file = fullfile(save_dir, 'sensor_disp_trend_table.tsv');

if ~isfile(tbl_file)
    error('Trend table not found: %s\nRun pt_plot_displacement first.', tbl_file);
end

T = readtable(tbl_file, 'Delimiter', '\t', 'FileType', 'text');

% Load n_axes from sensor r² file (not defined by config_pert)
sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
if ~isfile(sen_file)
    error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
end
tmp = load(sen_file, 'n_axes');
n_axes = tmp.n_axes;
clear tmp

fprintf('=== SLOPE VS POSITION PLOTS ===\n');
fprintf('  Loaded %d rows from %s\n', height(T), tbl_file);

% All source positions present in the table (sorted)
all_src_mm = unique(T.SourcePosition_mm);
all_src_mm = sort(all_src_mm);

% Forward models present in the table (preserving config order where possible)
methods_in_table = unique(T.Method);
% Build ordered list: config order first, then any extras
ordered_methods = {};
ordered_colors  = [];
for k = 1:numel(fwd_method_labels)
    lbl = fwd_method_labels{k};
    if any(strcmp(methods_in_table, lbl))
        ordered_methods{end+1} = lbl; %#ok<AGROW>
        ordered_colors(end+1,:) = fwd_method_colors(k,:); %#ok<AGROW>
    end
end
% Any method in table not in config gets a grey fallback
for k = 1:numel(methods_in_table)
    if ~any(strcmp(ordered_methods, methods_in_table{k}))
        ordered_methods{end+1} = methods_in_table{k}; %#ok<AGROW>
        ordered_colors(end+1,:) = [0.5 0.5 0.5]; %#ok<AGROW>
    end
end
n_methods = numel(ordered_methods);

% Line styles to distinguish methods if colours are similar
line_styles = {'-', '--', ':', '-.'};

% One figure per sensor axis
for sens_ax = 1:n_axes
    fig = figure('Color', 'w', 'Position', [100, 100, 1400, 420]);
    tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
    title(tl, sprintf('Rate of r² change per unit displacement — Sensor axis %d of %d', ...
        sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');
    xlabel(tl, 'Distance along the spinal cord (mm)', 'FontSize', 12);
    ylabel(tl, 'Slope  (r² mm⁻¹)', 'FontSize', 12);

    leg_h   = gobjects(n_methods, 1);
    leg_lbl = ordered_methods;

    for ori_idx = 1:n_ori
        ori_name = orientation_display{ori_idx};
        ax_panel = nexttile(tl, ori_idx);
        hold(ax_panel, 'on');

        for m_idx = 1:n_methods
            mlabel = ordered_methods{m_idx};
            col    = ordered_colors(m_idx,:);
            lstyle = line_styles{mod(m_idx-1, numel(line_styles)) + 1};

            % Filter table rows for this method, orientation, sensor axis
            row_mask = strcmp(T.Method, mlabel) & ...
                       strcmp(T.Orientation, ori_name) & ...
                       T.SensorAxis == sens_ax;
            sub = T(row_mask, :);

            if isempty(sub)
                leg_h(m_idx) = plot(ax_panel, NaN, NaN, lstyle, ...
                    'Color', col, 'LineWidth', 2.0);
                continue
            end

            % Sort by source position
            [src_sorted, si] = sort(sub.SourcePosition_mm);
            slope_sorted     = sub.Slope(si);

            h = plot(ax_panel, src_sorted, slope_sorted, lstyle, ...
                'Color', col, 'LineWidth', 2.0, ...
                'Marker', 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'w');
            leg_h(m_idx) = h;
        end

        yline(ax_panel, 0, '--k', 'LineWidth', 1.0, 'Alpha', 0.4);

        title(ax_panel, ori_name, 'FontSize', 12, 'FontWeight', 'bold');
        if ori_idx == 1
            ylabel(ax_panel, 'Slope  (r² mm⁻¹)', 'FontSize', 11);
        end
        xlabel(ax_panel, 'Distance along the spinal cord (mm)', 'FontSize', 11);
        grid(ax_panel, 'on');
        set(ax_panel, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
        hold(ax_panel, 'off');
    end

    % Shared legend on rightmost panel
    ax_last = nexttile(tl, n_ori);
    lgd = legend(ax_last, leg_h, leg_lbl, ...
        'Location', 'eastoutside', 'FontSize', 11);
    lgd.Box = 'off';
    title(lgd, 'Forward model');

    fname = sprintf('slope_vs_position_sensorax%d', sens_ax);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

fprintf('\npt_plot_slope_vs_position complete.\n');
fprintf('Figures saved to: %s\n', save_dir);
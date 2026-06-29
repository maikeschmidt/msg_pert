% pt_plot_displacement - Median sensor displacement vs r² plots (sensor mode)
%
% Sensor mode only. For selected source points along the cord, plots r²
% against the median absolute sensor displacement per shift realisation.
% Requires sensor_shift_vectors to be set in config_pert.m.
%
% USAGE:
%   pt_plot_displacement
%
% DEPENDENCIES:
%   config_pert, pert_sensor_rsq.mat
%
% OUTPUTS (saved to <save_base_dir>/perturbation_analysis/sensor/):
%   sensor_disp_vs_rsq_source<N>mm.png/.fig   — per source position
%   sensor_disp_vs_rsq_combined_sensorax<N>.png/.fig — all sources per sensor axis
%   sensor_disp_trend_table.csv               — line-fit coefficients per source
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

target_min_mm = 5;   % SET THIS: lower bound of source selection (mm)
target_max_mm = 600;   % SET THIS: upper bound of source selection (mm)

% INITIALISE

config_pert;
pt_add_functions;

n_ori = numel(orientation_labels);

sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
if ~isfile(sen_file)
    error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
end
load(sen_file);   %#ok<LOAD>

save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

fprintf('=== DISPLACEMENT vs R² PLOTS ===\n');

% Compute median displacement per shift realisation
if ~isempty(sensor_shift_vectors) && numel(sensor_shift_vectors) == n_sensor_bundles
    median_displacements = nan(1, numel(valid_keys_geom));
    for b = 1:n_sensor_bundles
        bundle_mask = valid_bundle_idx == b;
        bund_rows   = find(bundle_mask);
        for i = 1:numel(bund_rows)
            s   = valid_shift_idx(bund_rows(i));
            vec = sensor_shift_vectors{b}(s, :);
            median_displacements(bund_rows(i)) = median(abs(vec));
        end
    end
    x_label_disp    = 'Median sensor displacement (mm)';
    use_actual_disp = true;
    fprintf('  Using actual shift vectors for x-axis.\n');
else
    warning(['sensor_shift_vectors not set in config_pert.m — ' ...
             'using shift index as x-axis. ' ...
             'Paste shift vectors from pt_generate_sensor_shifts output.']);
    median_displacements = 1:numel(valid_keys_geom);
    x_label_disp         = 'Shift realisation index';
    use_actual_disp      = false;
end

% Source selection
cord_positions_mm = src_range * src_spacing_mm;
source_sel_mask   = cord_positions_mm >= target_min_mm & ...
                    cord_positions_mm <= target_max_mm;
source_sel_idx    = find(source_sel_mask);
source_sel_mm     = cord_positions_mm(source_sel_idx);
n_sources_sel     = numel(source_sel_idx);

fprintf('  Source selection: %d–%d mm — %d source points\n', ...
    target_min_mm, target_max_mm, n_sources_sel);
if n_sources_sel == 0
    error('No sources found between %d and %d mm.', target_min_mm, target_max_mm);
end

% Bundle x-axis shading regions and overall x range
if use_actual_disp
    bundle_x_ranges = zeros(n_sensor_bundles, 2);
    for b = 1:n_sensor_bundles
        vals = median_displacements(valid_bundle_idx == b);
        bundle_x_ranges(b,:) = [min(vals)*0.85, max(vals)*1.15];
    end
    x_max = max(median_displacements) * 1.2;
else
    x_max = numel(valid_keys_geom) + 1;
end

bundle_shade_colors = [0.85, 0.93, 0.97; 0.75, 0.87, 0.95; 0.65, 0.78, 0.90];
src_cmap = parula(n_sources_sel);
[sorted_disp, sort_idx] = sort(median_displacements);

method_label_map = containers.Map(fwd_methods, fwd_method_labels);
n_loaded_methods = numel(loaded_methods);


%% =========================================================================
%% INDIVIDUAL SOURCE FIGURES
%% =========================================================================
% Layout: n_axes rows x n_ori cols.
% Colour = bundle. Y-axis: same limits per row across all orientation cols.
% Legend: horizontal rows of 4 at the bottom.

% fprintf('Generating individual source figures...\n');
% 
% bundle_markers_ind = {'o', 's', '^'};
% 
% for src_sel = 1:n_sources_sel
%     src_plot_idx = source_sel_idx(src_sel);
%     src_mm       = source_sel_mm(src_sel);
% 
%     % Pre-compute per-row y-limits (same across all ori cols for that row)
%     row_ylim = zeros(n_axes, 2);
%     for sens_ax = 1:n_axes
%         all_y = [];
%         for m_idx = 1:n_loaded_methods
%             method    = loaded_methods{m_idx};
%             rsq_store = rsq_by_method.(method);
%             for ori_idx = 1:n_ori
%                 ori_label = orientation_labels{ori_idx};
%                 vals = squeeze(rsq_store.(ori_label)(:, src_plot_idx, sens_ax));
%                 all_y = [all_y; vals(:)]; %#ok<AGROW>
%             end
%         end
%         all_y(isnan(all_y)) = [];
%         if isempty(all_y)
%             row_ylim(sens_ax,:) = [0, 1];
%         else
%             ylo = max(0, min(all_y) - 0.02);
%             yhi = min(1, max(all_y) + 0.02);
%             if yhi - ylo < 0.01; yhi = ylo + 0.01; end
%             row_ylim(sens_ax,:) = [ylo, yhi];
%         end
%     end
% 
%     fig = figure('Color', 'w', 'Position', [100, 100, 1400, 1050]);
%     tl  = tiledlayout(n_axes, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
%     title(tl, sprintf('Sensor Shift — Source at %d mm along spinal cord', src_mm), ...
%         'FontSize', 14, 'FontWeight', 'bold');
% 
%     % Collect legend handles: bundle x method
%     n_leg      = n_sensor_bundles * n_loaded_methods;
%     leg_h      = gobjects(n_leg, 1);
%     leg_lbl    = cell(n_leg, 1);
%     leg_filled = false(n_leg, 1);
% 
%     for sens_ax = 1:n_axes
%         for ori_idx = 1:n_ori
%             ori_label = orientation_labels{ori_idx};
%             ax_panel  = nexttile(tl, (sens_ax-1)*n_ori + ori_idx);
%             hold(ax_panel, 'on');
% 
%             if use_actual_disp
%                 for b = 1:n_sensor_bundles
%                     xr = bundle_x_ranges(b,:);
%                     patch(ax_panel, [xr(1) xr(2) xr(2) xr(1)], ...
%                         [row_ylim(sens_ax,1) row_ylim(sens_ax,1) ...
%                          row_ylim(sens_ax,2) row_ylim(sens_ax,2)], ...
%                         bundle_shade_colors(b,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%                 end
%             end
% 
%             for m_idx = 1:n_loaded_methods
%                 method    = loaded_methods{m_idx};
%                 rsq_store = rsq_by_method.(method);
%                 mcol      = fwd_method_colors(m_idx, :);
% 
%                 for b = 1:n_sensor_bundles
%                     bund_rows = find(valid_bundle_idx == b);
%                     bmarker   = bundle_markers_ind{min(b, numel(bundle_markers_ind))};
%                     x_vals    = median_displacements(bund_rows);
%                     y_vals    = squeeze(rsq_store.(ori_label)(bund_rows, src_plot_idx, sens_ax))';
%                     h = scatter(ax_panel, x_vals, y_vals, 55, ...
%                         'Marker', bmarker, ...
%                         'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', ...
%                         'LineWidth', 0.8);
% 
%                     li = (b-1)*n_loaded_methods + m_idx;
%                     if ~leg_filled(li)
%                         leg_h(li)      = h;
%                         leg_lbl{li}    = sprintf('%s — %s', ...
%                             method_label_map(method), sensor_bundle_display{b});
%                         leg_filled(li) = true;
%                     end
%                 end
%             end
% 
%             yline(ax_panel, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4);
%             yline(ax_panel, 0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.4 0.4 0.4]);
%             yline(ax_panel, 0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.6 0.6 0.6]);
% 
%             if sens_ax == 1
%                 title(ax_panel, orientation_display{ori_idx}, 'FontSize', 12, 'FontWeight', 'bold');
%             end
%             if ori_idx == 1
%                 ylabel(ax_panel, sprintf('Sensor axis %d\nr²', sens_ax), 'FontSize', 11);
%             end
%             if sens_ax == n_axes
%                 xlabel(ax_panel, x_label_disp, 'FontSize', 11);
%             end
% 
%             xlim(ax_panel, [0, x_max]);
%             ylim(ax_panel, row_ylim(sens_ax,:));
%             grid(ax_panel, 'on');
%             set(ax_panel, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
%             hold(ax_panel, 'off');
%         end
%     end
% 
%     % Legend: 4 entries per row at the bottom
%     lgd = legend(leg_h, leg_lbl, ...
%         'Orientation', 'horizontal', 'NumColumns', 4, ...
%         'FontSize', 10, 'Box', 'off');
%     lgd.Layout.Tile = 'south';
% 
%         fname = sprintf('sensor_disp_vs_rsq_source%dmm', src_mm);
%     % exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
%     % saveas(fig, fullfile(save_dir, [fname '.fig']));
%     % close(fig);
%     fprintf('  Saved: %s\n', fname);
% end


%% =========================================================================
%% COMBINED FIGURES (one per sensor axis)
%% =========================================================================
% Layout: n_loaded_methods rows x n_ori cols.
% Each row = one forward model; each col = one dipole orientation.
% Y-axis auto-scaled per tile. Coloured by source position (parula).
% No in-figure r annotations — saved to CSV table instead.

fprintf('\nGenerating combined figures...\n');

trend_rows = {};

for sens_ax = 1:n_axes
    fig_h = max(500, n_loaded_methods * 320 + 150);
    fig = figure('Color', 'w', 'Position', [100, 100, 1800, fig_h]);
    tl  = tiledlayout(n_loaded_methods, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
    title(tl, sprintf('Sensor Shift — Sensor axis %d of %d  |  Sources %d–%d mm', ...
        sens_ax, n_axes, target_min_mm, target_max_mm), ...
        'FontSize', 13, 'FontWeight', 'bold');
    xlabel(tl, x_label_disp, 'FontSize', 12);

    for m_idx = 1:n_loaded_methods
        method    = loaded_methods{m_idx};
        rsq_store = rsq_by_method.(method);
        mlabel    = method_label_map(method);

        for ori_idx = 1:n_ori
            ori_label = orientation_labels{ori_idx};
            ax_panel  = nexttile(tl, (m_idx-1)*n_ori + ori_idx);
            hold(ax_panel, 'on');

            % Auto y-limits for this tile
            all_y = [];
            for src_sel = 1:n_sources_sel
                src_plot_idx = source_sel_idx(src_sel);
                vals = squeeze(rsq_store.(ori_label)(sort_idx, src_plot_idx, sens_ax));
                all_y = [all_y; vals(:)]; %#ok<AGROW>
            end
            all_y(isnan(all_y)) = [];
            if isempty(all_y)
                panel_ylim = [0, 1];
            else
                ylo = max(0, min(all_y) - 0.02);
                yhi = min(1, max(all_y) + 0.02);
                if yhi - ylo < 0.01; yhi = ylo + 0.01; end
                panel_ylim = [ylo, yhi];
            end

            if use_actual_disp
                for b = 1:n_sensor_bundles
                    xr = bundle_x_ranges(b,:);
                    patch(ax_panel, [xr(1) xr(2) xr(2) xr(1)], ...
                        [panel_ylim(1) panel_ylim(1) panel_ylim(2) panel_ylim(2)], ...
                        bundle_shade_colors(b,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                    if m_idx == 1
                        text(ax_panel, mean(xr), panel_ylim(2) - 0.002, sensor_bundle_display{b}, ...
                            'HorizontalAlignment', 'center', 'FontSize', 9, ...
                            'Color', sensor_bundle_colors(b,:), 'FontWeight', 'bold', ...
                            'VerticalAlignment', 'top');
                    end
                end
            end

            leg_h = gobjects(n_sources_sel, 1);
            for src_sel = 1:n_sources_sel
                src_plot_idx = source_sel_idx(src_sel);
                src_mm_val   = source_sel_mm(src_sel);
                col          = src_cmap(src_sel,:);

                rsq_vals_sorted = squeeze(rsq_store.(ori_label)(sort_idx, src_plot_idx, sens_ax));

                scatter(ax_panel, sorted_disp, rsq_vals_sorted, 35, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
                    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);

                valid_pts = ~isnan(sorted_disp) & ~isnan(rsq_vals_sorted');
                if sum(valid_pts) >= 2
                    p     = polyfit(sorted_disp(valid_pts), rsq_vals_sorted(valid_pts), 1);
                    x_fit = linspace(min(sorted_disp), max(sorted_disp), 200);
                    y_fit = polyval(p, x_fit);
                    leg_h(src_sel) = plot(ax_panel, x_fit, y_fit, '-', ...
                        'Color', col, 'LineWidth', 2.0);
                    r_val = corr(sorted_disp(valid_pts)', rsq_vals_sorted(valid_pts), 'Type', 'Pearson');
                    trend_rows{end+1} = {src_mm_val, orientation_display{ori_idx}, sens_ax, ...
                        mlabel, p(1), p(2), r_val}; %#ok<AGROW>
                else
                    leg_h(src_sel) = plot(ax_panel, NaN, NaN, '-', 'Color', col, 'LineWidth', 2.0);
                    trend_rows{end+1} = {src_mm_val, orientation_display{ori_idx}, sens_ax, ...
                        mlabel, NaN, NaN, NaN}; %#ok<AGROW>
                end
            end

            yline(ax_panel, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
            yline(ax_panel, 0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                'Color', [0.4 0.4 0.4], 'Label', 'r²=0.99', ...
                'LabelHorizontalAlignment', 'left', 'FontSize', 9);
            yline(ax_panel, 0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.4, ...
                'Color', [0.6 0.6 0.6], 'Label', 'r²=0.95', ...
                'LabelHorizontalAlignment', 'left', 'FontSize', 9);

            if m_idx == 1
                title(ax_panel, orientation_display{ori_idx}, 'FontSize', 13, 'FontWeight', 'bold');
            end
            if ori_idx == 1
                ylabel(ax_panel, {mlabel; 'r²'}, 'FontSize', 12);
            end

            xlim(ax_panel, [0, x_max]);
            ylim(ax_panel, panel_ylim);
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');

            if ori_idx == n_ori
                lgd = legend(ax_panel, leg_h, ...
                    arrayfun(@(mm) sprintf('%d mm', mm), source_sel_mm, 'UniformOutput', false), ...
                    'Location', 'eastoutside', 'FontSize', 11);
                lgd.Box = 'off';
                title(lgd, 'Source (mm)');
            end
            hold(ax_panel, 'off');
        end
    end

    fname = sprintf('sensor_disp_vs_rsq_combined_sensorax%d', sens_ax);
    % exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
    % saveas(fig, fullfile(save_dir, [fname '.fig']));
    % close(fig);
    fprintf('  Saved: %s\n', fname);
end


%% =========================================================================
%% TREND TABLE
%% =========================================================================
% One row per source x orientation x sensor axis.
% Columns: source position, orientation, sensor axis, method, slope, intercept, r.

fprintf('\nSaving trend table...\n');

if ~isempty(trend_rows)
    src_pos_col   = cellfun(@(r) r{1}, trend_rows)';
    ori_col       = cellfun(@(r) r{2}, trend_rows, 'UniformOutput', false)';
    saxes_col     = cellfun(@(r) r{3}, trend_rows)';
    method_col    = cellfun(@(r) r{4}, trend_rows, 'UniformOutput', false)';
    slope_col     = cellfun(@(r) r{5}, trend_rows)';
    intercept_col = cellfun(@(r) r{6}, trend_rows)';
    r_col         = cellfun(@(r) r{7}, trend_rows)';

    T = table(src_pos_col, ori_col, saxes_col, method_col, ...
              slope_col, intercept_col, r_col, ...
        'VariableNames', {'SourcePosition_mm', 'Orientation', 'SensorAxis', 'Method', ...
                          'Slope', 'Intercept', 'PearsonR'});

    tbl_file = fullfile(save_dir, 'sensor_disp_trend_table.tsv');
    writetable(T, tbl_file, 'Delimiter', '\t', 'FileType', 'text');
    fprintf('  Saved: %s\n', tbl_file);
    disp(T);
else
    fprintf('  No trend data collected — table skipped.\n');
end

fprintf('\npt_plot_displacement complete.\n');
fprintf('Figures saved to: %s\n', save_dir);
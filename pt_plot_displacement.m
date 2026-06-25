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
%   sensor_disp_vs_rsq_source<N>mm.png/.fig
%   sensor_disp_vs_rsq_combined_sensorax<N>.png/.fig
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

target_min_mm = 40;   % SET THIS: lower bound of source selection (mm)
target_max_mm = 85;   % SET THIS: upper bound of source selection (mm)

% INITIALISE

config_pert;
pt_add_functions;

sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
if ~isfile(sen_file)
    error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
end
load(sen_file);

save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

fprintf('=== DISPLACEMENT vs R² PLOTS ===\n');

n_valid_shifts = numel(valid_keys_geom);

% Compute median displacement per shift realisation
if ~isempty(sensor_shift_vectors) && numel(sensor_shift_vectors) == n_sensor_bundles
    median_displacements = nan(1, n_valid_shifts);
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
             'using shift index as x-axis.\n' ...
             'Paste shift vectors from pt_generate_sensor_shifts output.']);
    median_displacements = 1:n_valid_shifts;
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
    x_max = numel(valid_keys) + 1;
end

bundle_shade_colors = [0.85, 0.93, 0.97; 0.75, 0.87, 0.95; 0.65, 0.78, 0.90];
src_cmap = parula(n_sources_sel);
[sorted_disp, sort_idx] = sort(median_displacements);

method_label_map = containers.Map(fwd_methods, fwd_method_labels);
n_loaded_methods = numel(loaded_methods);

for m_idx = 1:n_loaded_methods
    method    = loaded_methods{m_idx};
    rsq_store = rsq_by_method.(method);
    mlabel    = method_label_map(method);
    fprintf('\n[%s]\n', mlabel);

    %% Individual source figures

    fprintf('  Generating individual source figures...\n');

    for src_sel = 1:n_sources_sel
        src_plot_idx = source_sel_idx(src_sel);
        src_mm       = source_sel_mm(src_sel);

        fig = figure('Color', 'w', 'Position', [100, 100, 1400, 1050]);
        tl  = tiledlayout(n_axes, numel(orientation_labels), ...
            'TileSpacing', 'compact', 'Padding', 'loose');
        title(tl, sprintf('[%s]  Sensor Shift — Source at %d mm along spinal cord', ...
            mlabel, src_mm), 'FontSize', 14, 'FontWeight', 'bold');
        leg_h_bundles = gobjects(n_sensor_bundles, 1);

        for sens_ax = 1:n_axes
            for ori_idx = 1:numel(orientation_labels)
                ori_label = orientation_labels{ori_idx};
                ax_panel  = nexttile(tl);
                hold(ax_panel, 'on');

                if use_actual_disp
                    for b = 1:n_sensor_bundles
                        xr = bundle_x_ranges(b,:);
                        patch(ax_panel, [xr(1) xr(2) xr(2) xr(1)], [0 0 1.05 1.05], ...
                            bundle_shade_colors(b,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                    end
                end

                for b = 1:n_sensor_bundles
                    bund_rows = find(valid_bundle_idx == b);
                    col = sensor_bundle_colors(b,:);
                    x_vals = median_displacements(bund_rows);
                    y_vals = squeeze(rsq_store.(ori_label)(bund_rows, src_plot_idx, sens_ax))';
                    leg_h_bundles(b) = scatter(ax_panel, x_vals, y_vals, 60, ...
                        'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
                        'LineWidth', 0.8, 'DisplayName', sensor_bundle_display{b});
                end

                yline(ax_panel, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4);
                yline(ax_panel, 0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.4 0.4 0.4]);
                yline(ax_panel, 0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.6 0.6 0.6]);

                if sens_ax == 1; title(ax_panel, orientation_display{ori_idx}, 'FontSize', 12, 'FontWeight', 'bold'); end
                if ori_idx == 1; ylabel(ax_panel, sprintf('Sensor axis %d\nr²', sens_ax), 'FontSize', 11); end
                if sens_ax == n_axes; xlabel(ax_panel, x_label_disp, 'FontSize', 11); end

                xlim(ax_panel, [0, x_max]); ylim(ax_panel, [0, 1.05]);
                grid(ax_panel, 'on');
                set(ax_panel, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
                hold(ax_panel, 'off');
            end
        end

        lgd = legend(leg_h_bundles, sensor_bundle_display, 'Orientation', 'horizontal', ...
            'FontSize', 11, 'Box', 'off');
        lgd.Layout.Tile = 'south';

        fname = sprintf('sensor_disp_vs_rsq_%s_source%dmm', method, src_mm);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end

    %% Combined figures (one per sensor axis)

    fprintf('  Generating combined figures...\n');

    for sens_ax = 1:n_axes
        fig = figure('Color', 'w', 'Position', [100, 100, 1800, 650]);
        tl  = tiledlayout(1, numel(orientation_labels), 'TileSpacing', 'compact', 'Padding', 'loose');
        title(tl, sprintf('[%s]  Sensor Shift — Sensor axis %d of %d  |  Sources %d–%d mm', ...
            mlabel, sens_ax, n_axes, target_min_mm, target_max_mm), ...
            'FontSize', 13, 'FontWeight', 'bold');

        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};
            ax_panel  = nexttile(tl);
            hold(ax_panel, 'on');

            if use_actual_disp
                for b = 1:n_sensor_bundles
                    xr = bundle_x_ranges(b,:);
                    patch(ax_panel, [xr(1) xr(2) xr(2) xr(1)], [0 0 1.05 1.05], ...
                        bundle_shade_colors(b,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                    text(ax_panel, mean(xr), 1.03, sensor_bundle_display{b}, ...
                        'HorizontalAlignment', 'center', 'FontSize', 9, ...
                        'Color', sensor_bundle_colors(b,:), 'FontWeight', 'bold');
                end
            end

            leg_h = gobjects(n_sources_sel, 1);
            for src_sel = 1:n_sources_sel
                src_plot_idx = source_sel_idx(src_sel);
                src_mm       = source_sel_mm(src_sel);
                col          = src_cmap(src_sel,:);

                rsq_vals_sorted = squeeze(rsq_store.(ori_label)(sort_idx, src_plot_idx, sens_ax));
                scatter(ax_panel, sorted_disp, rsq_vals_sorted, 35, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
                    'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);

                valid_pts = ~isnan(sorted_disp) & ~isnan(rsq_vals_sorted');
                if sum(valid_pts) >= 2
                    p     = polyfit(sorted_disp(valid_pts), rsq_vals_sorted(valid_pts), 1);
                    x_fit = linspace(min(sorted_disp), max(sorted_disp), 200);
                    y_fit = polyval(p, x_fit);
                    leg_h(src_sel) = plot(ax_panel, x_fit, y_fit, '-', ...
                        'Color', col, 'LineWidth', 2.0, ...
                        'DisplayName', sprintf('%d mm', src_mm));
                    r_val = corr(sorted_disp(valid_pts)', rsq_vals_sorted(valid_pts), 'Type', 'Pearson');
                    y_ann = 0.06 + (src_sel - 1) * 0.065;
                    text(ax_panel, x_max * 0.02, y_ann, sprintf('%d mm: r=%.2f', src_mm, r_val), ...
                        'FontSize', 8, 'Color', col, 'FontWeight', 'bold');
                else
                    leg_h(src_sel) = plot(ax_panel, NaN, NaN, '-', 'Color', col, 'LineWidth', 2.0, ...
                        'DisplayName', sprintf('%d mm', src_mm));
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

            title(ax_panel, orientation_display{ori_idx}, 'FontSize', 13, 'FontWeight', 'bold');
            xlabel(ax_panel, x_label_disp, 'FontSize', 13);
            if ori_idx == 1
                ylabel(ax_panel, {'r²  vs original leadfield'; '(1.0 = no effect)'}, 'FontSize', 12);
            end

            xlim(ax_panel, [0, x_max]); ylim(ax_panel, [0, 1.05]);
            grid(ax_panel, 'on');
            set(ax_panel, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');

            lgd = legend(ax_panel, leg_h, ...
                arrayfun(@(mm) sprintf('%d mm along cord', mm), source_sel_mm, 'UniformOutput', false), ...
                'Location', 'eastoutside', 'FontSize', 11);
            lgd.Box = 'off';
            title(lgd, 'Source position');
            hold(ax_panel, 'off');
        end

        fname = sprintf('sensor_disp_vs_rsq_%s_combined_sensorax%d', method, sens_ax);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end

    fprintf('  [%s] Done\n', mlabel);
end

fprintf('\npt_plot_displacement complete.\n');
fprintf('Figures saved to: %s\n', save_dir);

% pt_plot_displacement - Displacement / perturbation vs r² plots
%
% For sensor and source shifts: x-axis is median absolute displacement (mm).
% For conductivity perturbation: x-axis is mean % conductivity change.
%
% Produces:
%   Individual per-source figures  — cervical region only (cervical_min/max_mm)
%   Combined figure (all sources)  — full cord, coloured by source position
%   Trend table (TSV)              — full cord, used by pt_plot_slope_vs_position
%
% USAGE:
%   pt_plot_displacement
%
% DEPENDENCIES:
%   config_pert, pert_sensor_rsq.mat, pert_source_rsq.mat, pert_cond_rsq.mat
%
% OUTPUTS (saved to <save_base_dir>/perturbation_analysis/<mode>/):
%   <mode>_disp_vs_rsq_source<N>mm.png/.fig   — individual (cervical only)
%   <mode>_disp_vs_rsq_combined_sensorax<N>.png/.fig
%   <mode>_disp_trend_table.tsv
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

cervical_min_mm = 40;   % individual per-source figures produced only for this range
cervical_max_mm = 85;

% INITIALISE

config_pert;
pt_add_functions;

n_ori = numel(orientation_labels);

method_label_map = containers.Map(fwd_methods, fwd_method_labels);
bundle_shade_colors = [0.85, 0.93, 0.97; 0.75, 0.87, 0.95; 0.65, 0.78, 0.90];


%% =========================================================================
%% LOCAL FUNCTION: run one mode
%% =========================================================================
% Defined inline below — called for each of the three modes.


%% =========================================================================
%% SENSOR MODE
%% =========================================================================

if run_sensor
    fprintf('=== SENSOR DISPLACEMENT vs R² ===\n');

    sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
    if ~isfile(sen_file)
        error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
    end
    load(sen_file);   %#ok<LOAD>

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    % Displacement per realisation (median of abs shift vector)
    if ~isempty(sensor_shift_vectors) && numel(sensor_shift_vectors) == n_sensor_bundles
        median_displacements = nan(1, numel(valid_keys_geom));
        for b = 1:n_sensor_bundles
            bund_rows = find(valid_bundle_idx == b);
            for i = 1:numel(bund_rows)
                s   = valid_shift_idx(bund_rows(i));
                vec = sensor_shift_vectors{b}(s, :);
                median_displacements(bund_rows(i)) = median(abs(vec));
            end
        end
        x_label_str  = 'Median sensor displacement (mm)';
        use_actual   = true;
    else
        warning('sensor_shift_vectors not set — using shift index as x-axis.');
        median_displacements = 1:numel(valid_keys_geom);
        x_label_str  = 'Shift realisation index';
        use_actual   = false;
    end

    bundle_display_use = sensor_bundle_display;
    bundle_colors_use  = sensor_bundle_colors;
    bundle_idx_use     = valid_bundle_idx;
    n_bundles_use      = n_sensor_bundles;

    run_disp_mode(rsq_by_method, loaded_methods, valid_keys_geom, ...
        median_displacements, x_label_str, use_actual, ...
        bundle_idx_use, n_bundles_use, bundle_display_use, bundle_colors_use, ...
        bundle_shade_colors, orientation_labels, orientation_display, ...
        src_range, src_spacing_mm, n_axes, distances, marker_idx, ...
        cervical_min_mm, cervical_max_mm, ...
        method_label_map, pub_line_width, pub_marker_size, ...
        save_dir, 'sensor', 'Sensor shift');
end


%% =========================================================================
%% SOURCE MODE
%% =========================================================================

if run_source
    fprintf('\n=== SOURCE DISPLACEMENT vs R² ===\n');

    src_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
    if ~isfile(src_file)
        error('Source r² file not found: %s\nRun pt_compute_rsq first.', src_file);
    end
    load(src_file);   %#ok<LOAD>

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'source');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    if ~isempty(source_shift_vectors) && numel(source_shift_vectors) == n_source_bundles
        median_displacements = nan(1, numel(valid_keys_geom));
        for b = 1:n_source_bundles
            bund_rows = find(valid_source_bundle_idx == b);
            for i = 1:numel(bund_rows)
                s   = valid_source_shift_idx(bund_rows(i));
                vec = source_shift_vectors{b}(s, :);
                median_displacements(bund_rows(i)) = median(abs(vec));
            end
        end
        x_label_str = 'Median source displacement (mm)';
        use_actual  = true;
    else
        warning('source_shift_vectors not set — using shift index as x-axis.');
        median_displacements = 1:numel(valid_keys_geom);
        x_label_str = 'Shift realisation index';
        use_actual  = false;
    end

    bundle_display_use = source_bundle_display;
    bundle_colors_use  = source_bundle_colors;
    bundle_idx_use     = valid_source_bundle_idx;
    n_bundles_use      = n_source_bundles;

    run_disp_mode(rsq_by_method, loaded_methods, valid_keys_geom, ...
        median_displacements, x_label_str, use_actual, ...
        bundle_idx_use, n_bundles_use, bundle_display_use, bundle_colors_use, ...
        bundle_shade_colors, orientation_labels, orientation_display, ...
        src_range, src_spacing_mm, n_axes, distances, marker_idx, ...
        cervical_min_mm, cervical_max_mm, ...
        method_label_map, pub_line_width, pub_marker_size, ...
        save_dir, 'source', 'Source shift');
end


%% =========================================================================
%% CONDUCTIVITY MODE
%% =========================================================================

if run_cond
    fprintf('\n=== CONDUCTIVITY PERTURBATION vs R² ===\n');

    cond_file = fullfile(forward_fields_base, 'pert_cond_rsq.mat');
    if ~isfile(cond_file)
        error('Cond r² file not found: %s\nRun pt_compute_rsq first.', cond_file);
    end
    load(cond_file);   %#ok<LOAD>

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'cond');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    % x-axis: mean % conductivity change per realisation (saved by pt_compute_rsq)
    median_displacements = valid_cond_mean_pct;
    x_label_str          = 'Mean conductivity increase (%)';
    use_actual           = true;

    bundle_display_use = cond_bundle_display;
    bundle_colors_use  = cond_bundle_colors;
    bundle_idx_use     = valid_cond_bundle_idx;
    n_bundles_use      = n_cond_bundles;

    % Only one method for cond (BEM)
    run_disp_mode(rsq_by_method, loaded_methods, valid_cond_keys, ...
        median_displacements, x_label_str, use_actual, ...
        bundle_idx_use, n_bundles_use, bundle_display_use, bundle_colors_use, ...
        bundle_shade_colors, orientation_labels, orientation_display, ...
        src_range, src_spacing_mm, n_axes, distances, marker_idx, ...
        cervical_min_mm, cervical_max_mm, ...
        method_label_map, pub_line_width, pub_marker_size, ...
        save_dir, 'cond', 'Conductivity perturbation');
end

fprintf('\npt_plot_displacement complete.\n');


%% =========================================================================
%% LOCAL FUNCTION
%% =========================================================================

function run_disp_mode(rsq_by_method, loaded_methods, valid_keys, ...
    median_displacements, x_label_str, use_actual, ...
    bundle_idx_use, n_bundles_use, bundle_display_use, bundle_colors_use, ...
    bundle_shade_colors, orientation_labels, orientation_display, ...
    src_range, src_spacing_mm, n_axes, distances, marker_idx, ...
    cervical_min_mm, cervical_max_mm, ...
    method_label_map, pub_line_width, pub_marker_size, ...
    save_dir, mode_name, mode_title)

n_ori            = numel(orientation_labels);
n_loaded_methods = numel(loaded_methods);
cord_positions_mm = src_range * src_spacing_mm;

% All source positions (for trend table + combined figure)
all_src_idx = 1:numel(cord_positions_mm);
all_src_mm  = cord_positions_mm;

% Cervical source positions (for individual figures only)
cerv_mask   = cord_positions_mm >= cervical_min_mm & cord_positions_mm <= cervical_max_mm;
cerv_idx    = find(cerv_mask);
cerv_mm     = cord_positions_mm(cerv_idx);
n_cerv      = numel(cerv_idx);

fprintf('  Full cord: %d source positions\n', numel(all_src_mm));
fprintf('  Cervical (%d–%d mm): %d source positions (individual figures)\n', ...
    cervical_min_mm, cervical_max_mm, n_cerv);

% Bundle x-axis shading
if use_actual
    bundle_x_ranges = zeros(n_bundles_use, 2);
    for b = 1:n_bundles_use
        vals = median_displacements(bundle_idx_use == b);
        bundle_x_ranges(b,:) = [min(vals)*0.85, max(vals)*1.15];
    end
    x_max = max(median_displacements) * 1.2;
else
    x_max = numel(valid_keys) + 1;
    bundle_x_ranges = zeros(n_bundles_use, 2);
end

src_cmap = parula(max(numel(all_src_mm), 2));
[sorted_disp, sort_idx] = sort(median_displacements);

trend_rows = {};

% ----------------------------------------------------------------
% INDIVIDUAL PER-SOURCE FIGURES — cervical range only
% ----------------------------------------------------------------
if n_cerv > 0
    fprintf('  Generating individual source figures (cervical)...\n');
    bundle_markers_ind = {'o', 's', '^'};

    for src_sel = 1:n_cerv
        src_plot_idx = cerv_idx(src_sel);
        src_mm       = cerv_mm(src_sel);

        % Pre-compute y-limits per row (sensor axis)
        row_ylim = zeros(n_axes, 2);
        for sens_ax = 1:n_axes
            all_y = [];
            for m_idx = 1:n_loaded_methods
                method    = loaded_methods{m_idx};
                rsq_store = rsq_by_method.(method);
                for ori_idx = 1:n_ori
                    ori_label = orientation_labels{ori_idx};
                    vals = squeeze(rsq_store.(ori_label)(:, src_plot_idx, sens_ax));
                    all_y = [all_y; vals(:)]; %#ok<AGROW>
                end
            end
            all_y(isnan(all_y)) = [];
            if isempty(all_y)
                row_ylim(sens_ax,:) = [0, 1];
            else
                ylo = max(0, min(all_y) - 0.02);
                yhi = min(1, max(all_y) + 0.02);
                if yhi - ylo < 0.01; yhi = ylo + 0.01; end
                row_ylim(sens_ax,:) = [ylo, yhi];
            end
        end

        fig = figure('Color', 'w', 'Position', [100, 100, 1400, 1050]);
        tl  = tiledlayout(n_axes, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
        title(tl, sprintf('%s — Source at %d mm along spinal cord', mode_title, src_mm), ...
            'FontSize', 14, 'FontWeight', 'bold');

        n_leg      = n_bundles_use * n_loaded_methods;
        leg_h      = gobjects(n_leg, 1);
        leg_lbl    = cell(n_leg, 1);
        leg_filled = false(n_leg, 1);

        for sens_ax = 1:n_axes
            for ori_idx = 1:n_ori
                ori_label = orientation_labels{ori_idx};
                ax_panel  = nexttile(tl, (sens_ax-1)*n_ori + ori_idx);
                hold(ax_panel, 'on');

                if use_actual
                    for b = 1:n_bundles_use
                        xr = bundle_x_ranges(b,:);
                        patch(ax_panel, [xr(1) xr(2) xr(2) xr(1)], ...
                            [row_ylim(sens_ax,1) row_ylim(sens_ax,1) ...
                             row_ylim(sens_ax,2) row_ylim(sens_ax,2)], ...
                            bundle_shade_colors(min(b,3),:), ...
                            'EdgeColor', 'none', 'FaceAlpha', 0.5);
                    end
                end

                for m_idx = 1:n_loaded_methods
                    method    = loaded_methods{m_idx};
                    rsq_store = rsq_by_method.(method);
                    mcol      = get_method_color(method, fwd_methods, fwd_method_colors, m_idx);

                    for b = 1:n_bundles_use
                        bund_rows = find(bundle_idx_use == b);
                        bmarker   = bundle_markers_ind{min(b, numel(bundle_markers_ind))};
                        x_vals    = median_displacements(bund_rows);
                        y_vals    = squeeze(rsq_store.(ori_label)(bund_rows, src_plot_idx, sens_ax))';
                        h = scatter(ax_panel, x_vals, y_vals, 55, ...
                            'Marker', bmarker, ...
                            'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', ...
                            'LineWidth', 0.8);

                        li = (b-1)*n_loaded_methods + m_idx;
                        if ~leg_filled(li)
                            try; mlbl = method_label_map(method); catch; mlbl = method; end
                            leg_h(li)      = h;
                            leg_lbl{li}    = sprintf('%s — %s', mlbl, bundle_display_use{b});
                            leg_filled(li) = true;
                        end
                    end
                end

                yline(ax_panel, 1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.4);
                yline(ax_panel, 0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.4 0.4 0.4]);
                yline(ax_panel, 0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.4, 'Color', [0.6 0.6 0.6]);

                if sens_ax == 1
                    title(ax_panel, orientation_display{ori_idx}, 'FontSize', 12, 'FontWeight', 'bold');
                end
                if ori_idx == 1
                    ylabel(ax_panel, sprintf('Sensor axis %d\nr²', sens_ax), 'FontSize', 11);
                end
                if sens_ax == n_axes
                    xlabel(ax_panel, x_label_str, 'FontSize', 11);
                end

                xlim(ax_panel, [0, x_max]);
                ylim(ax_panel, row_ylim(sens_ax,:));
                grid(ax_panel, 'on');
                set(ax_panel, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
                hold(ax_panel, 'off');
            end
        end

        lgd = legend(leg_h(leg_filled), leg_lbl(leg_filled), ...
            'Orientation', 'horizontal', 'NumColumns', 4, ...
            'FontSize', 10, 'Box', 'off');
        lgd.Layout.Tile = 'south';

        fname = sprintf('%s_disp_vs_rsq_source%dmm', mode_name, src_mm);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end
end

% ----------------------------------------------------------------
% COMBINED FIGURES (one per sensor axis) — full cord
% ----------------------------------------------------------------
fprintf('  Generating combined figures (full cord)...\n');

n_src_all = numel(all_src_mm);
src_cmap  = parula(max(n_src_all, 2));

for sens_ax = 1:n_axes
    fig_h = max(500, n_loaded_methods * 320 + 150);
    fig   = figure('Color', 'w', 'Position', [100, 100, 1800, fig_h]);
    tl    = tiledlayout(n_loaded_methods, n_ori, 'TileSpacing', 'compact', 'Padding', 'loose');
    title(tl, sprintf('%s — Sensor axis %d of %d  |  Full cord', ...
        mode_title, sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');
    xlabel(tl, x_label_str, 'FontSize', 12);

    for m_idx = 1:n_loaded_methods
        method    = loaded_methods{m_idx};
        rsq_store = rsq_by_method.(method);
        try; mlabel = method_label_map(method); catch; mlabel = method; end

        for ori_idx = 1:n_ori
            ori_label = orientation_labels{ori_idx};
            ax_panel  = nexttile(tl, (m_idx-1)*n_ori + ori_idx);
            hold(ax_panel, 'on');

            all_y = [];
            for si = 1:n_src_all
                vals = squeeze(rsq_store.(ori_label)(sort_idx, all_src_idx(si), sens_ax));
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

            if use_actual
                for b = 1:n_bundles_use
                    xr = bundle_x_ranges(b,:);
                    patch(ax_panel, [xr(1) xr(2) xr(2) xr(1)], ...
                        [panel_ylim(1) panel_ylim(1) panel_ylim(2) panel_ylim(2)], ...
                        bundle_shade_colors(min(b,3),:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                    if m_idx == 1
                        text(ax_panel, mean(xr), panel_ylim(2) - 0.002, ...
                            bundle_display_use{b}, ...
                            'HorizontalAlignment', 'center', 'FontSize', 9, ...
                            'Color', bundle_colors_use(b,:), 'FontWeight', 'bold', ...
                            'VerticalAlignment', 'top');
                    end
                end
            end

            leg_h = gobjects(n_src_all, 1);
            for si = 1:n_src_all
                src_plot_idx    = all_src_idx(si);
                src_mm_val      = all_src_mm(si);
                col             = src_cmap(si,:);
                rsq_vals_sorted = squeeze(rsq_store.(ori_label)(sort_idx, src_plot_idx, sens_ax));

                valid_pts = ~isnan(sorted_disp) & ~isnan(rsq_vals_sorted');
                if sum(valid_pts) >= 2
                    p     = polyfit(sorted_disp(valid_pts), rsq_vals_sorted(valid_pts), 1);
                    x_fit = linspace(min(sorted_disp), max(sorted_disp), 200);
                    y_fit = polyval(p, x_fit);
                    leg_h(si) = plot(ax_panel, x_fit, y_fit, '-', ...
                        'Color', col, 'LineWidth', 2.0);
                    r_val = corr(sorted_disp(valid_pts)', rsq_vals_sorted(valid_pts), 'Type', 'Pearson');
                    trend_rows{end+1} = {src_mm_val, orientation_display{ori_idx}, sens_ax, ...
                        mlabel, p(1), p(2), r_val}; %#ok<AGROW>
                else
                    leg_h(si) = plot(ax_panel, NaN, NaN, '-', 'Color', col, 'LineWidth', 2.0);
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
                    arrayfun(@(mm) sprintf('%d mm', mm), all_src_mm, 'UniformOutput', false), ...
                    'Location', 'eastoutside', 'FontSize', 9);
                lgd.Box = 'off';
                title(lgd, 'Source (mm)');
            end
            hold(ax_panel, 'off');
        end
    end

    fname = sprintf('%s_disp_vs_rsq_combined_sensorax%d', mode_name, sens_ax);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('    Saved: %s\n', fname);
end

% ----------------------------------------------------------------
% TREND TABLE — full cord
% ----------------------------------------------------------------
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

    tbl_file = fullfile(save_dir, [mode_name '_disp_trend_table.tsv']);
    writetable(T, tbl_file, 'Delimiter', '\t', 'FileType', 'text');
    fprintf('  Saved trend table: %s\n', tbl_file);
end

end % run_disp_mode


function col = get_method_color(method, fwd_methods, fwd_method_colors, fallback_idx)
    idx = find(strcmp(fwd_methods, method), 1);
    if ~isempty(idx)
        col = fwd_method_colors(idx, :);
    else
        col = fwd_method_colors(min(fallback_idx, size(fwd_method_colors,1)), :);
    end
end

% pt_plot_curves - Plot r² summary figures for perturbation analysis
%
% Loads pre-computed r² data (from pt_compute_rsq) and produces concise
% summary figures. Individual shift results are stored in tables (see
% pt_compute_table); only summary figures are produced here.
%
% FIGURES PRODUCED
%   Source shift (detail) — per forward model method, per sensor axis:
%     source_detail_<method>_sensorax<N>.png
%       3×3 tiled layout: rows = bundle (1/2/3), cols = orientation (VD/RC/LR)
%       Each tile: 8 individual shift lines in bundle colour
%
%   Source shift (summary) — per forward model method, per sensor axis:
%     source_summary_<method>_sensorax<N>.png
%       1×3 tiled layout: cols = orientation; 3 mean lines (one per bundle)
%
%   Source shift (cross-model) — per sensor axis:
%     source_crossmodel_sensorax<N>.png
%       1×3 tiled layout: colour = method, marker = bundle
%
%   Sensor shift (detail) — per forward model method, per sensor axis:
%     sensor_detail_<method>_sensorax<N>.png
%       3×3 tiled layout: rows = bundle (1/2/3), cols = orientation
%       Each tile: 8 individual shift lines in bundle colour
%
%   Sensor shift (summary) — per forward model method, per sensor axis:
%     sensor_summary_<method>_sensorax<N>.png
%       1×3 tiled layout: cols = orientation
%       Each tile: 3 mean lines (one per bundle) in bundle colour
%
%   Sensor shift (cross-model) — per sensor axis:
%     sensor_crossmodel_sensorax<N>.png
%       1×3 tiled layout: cols = orientation
%       Each tile: n_methods × 3 lines (method line style, bundle colour shade)
%
% USAGE:
%   pt_plot_curves
%
% DEPENDENCIES:
%   config_pert, pert_source_rsq.mat, pert_sensor_rsq.mat
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
run_cond   = true;

config_pert;
pt_add_functions;

n_ori = numel(orientation_labels);


%% =========================================================================
%% SOURCE MODE
%% =========================================================================

if run_source
    fprintf('SOURCE PERTURBATION FIGURES\n');

    src_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
    if ~isfile(src_file)
        error('Source r² file not found: %s\nRun pt_compute_rsq first.', src_file);
    end
    load(src_file);   %#ok<LOAD>

    n_loaded_methods = numel(loaded_methods);
    if n_loaded_methods == 0
        warning('No method results in pert_source_rsq.mat — skipping source figures.');
    else

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'source');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    method_label_map = containers.Map(fwd_methods, fwd_method_labels);
    n_src_bundles    = numel(source_bundle_display);

    % ----------------------------------------------------------------
    % SOURCE DETAIL — per method, per sensor axis
    % Layout: n_src_bundles rows × n_ori cols; 8 lines per tile (one per shift)
    % ----------------------------------------------------------------
    for m_idx = 1:n_loaded_methods
        method    = loaded_methods{m_idx};
        rsq_store = rsq_by_method.(method);
        mlabel    = method_label_map(method);

        for sens_ax = 1:n_axes
            fig = figure('Color', 'w', 'Position', [50, 50, 1400, 960]);
            tl  = tiledlayout(n_src_bundles, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('[%s]  Source shift  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
            xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);

            for bund = 1:n_src_bundles
                bund_mask = valid_source_bundle_idx == bund;
                bund_rows = find(bund_mask);
                bund_col  = source_bundle_colors(bund, :);

                for ori_idx = 1:n_ori
                    ori_label = orientation_labels{ori_idx};
                    ax = nexttile(tl, (bund-1)*n_ori + ori_idx);
                    hold(ax, 'on');

                    n_in_bundle = numel(bund_rows);
                    tile_rsq = zeros(n_in_bundle, numel(distances));
                    for i = 1:n_in_bundle
                        tile_rsq(i,:) = squeeze(rsq_store.(ori_label)(bund_rows(i), :, sens_ax));
                    end
                    leg_handles = gobjects(n_in_bundle, 1);
                    for i = 1:n_in_bundle
                        t   = (i-1) / max(n_in_bundle-1, 1);
                        col = bund_col * (0.5 + 0.5*(1-t)) + (1-(0.5+0.5*(1-t)))*[1 1 1];
                        col = min(1, max(0, col));
                        leg_handles(i) = plot(ax, distances, tile_rsq(i,:), ...
                            '-', 'Color', col, 'LineWidth', pub_line_width - 0.5, ...
                            'Marker', 'o', 'MarkerIndices', marker_idx, ...
                            'MarkerSize', pub_marker_size - 1, ...
                            'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                    end

                    add_ref_lines(ax);
                    xlim(ax, [distances(1), distances(end)]);
                    y_lo = max(0,    floor(min(tile_rsq(:))*100)/100 - 0.01);
                    y_hi = min(1.05, ceil(max(tile_rsq(:))*100)/100 + 0.01);
                    ylim(ax, [y_lo, y_hi]);
                    xticks(ax, 0:20:ceil(distances(end)));
                    grid(ax, 'on');
                    set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');

                    if ori_idx == 1
                        ylabel(ax, {source_bundle_display{bund}; 'r²'}, ...
                            'FontSize', 11, 'Color', bund_col);
                    end
                    if bund == 1
                        title(ax, orientation_display{ori_idx}, 'FontSize', 12);
                    end
                    if ori_idx == n_ori
                        shift_labels = sensitivity_labels(bund_mask);
                        lgd = legend(ax, leg_handles, shift_labels(1:n_in_bundle), ...
                            'Location', 'eastoutside', 'FontSize', 9);
                        lgd.Box = 'off';
                        title(lgd, 'Shift');
                    end
                    hold(ax, 'off');
                end
            end

            fname = sprintf('source_detail_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('  Saved: %s\n', fname);
        end

        % ----------------------------------------------------------------
        % SOURCE SUMMARY — per method, per sensor axis
        % Layout: 1 row × n_ori cols; one mean line per bundle
        % ----------------------------------------------------------------
        for sens_ax = 1:n_axes
            fig = figure('Color', 'w', 'Position', [50, 50, 1400, 420]);
            tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('[%s]  Source shift summary  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
            xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);

            for ori_idx = 1:n_ori
                ori_label = orientation_labels{ori_idx};
                ax = nexttile(tl, ori_idx);
                hold(ax, 'on');

                all_tile = [];
                leg_h = gobjects(n_src_bundles, 1);
                for bund = 1:n_src_bundles
                    bund_rows = find(valid_source_bundle_idx == bund);
                    bund_col  = source_bundle_colors(bund, :);
                    bund_rsq  = zeros(numel(bund_rows), numel(distances));
                    for i = 1:numel(bund_rows)
                        bund_rsq(i,:) = squeeze(rsq_store.(ori_label)(bund_rows(i), :, sens_ax));
                    end
                    mean_rsq = mean(bund_rsq, 1, 'omitnan');
                    all_tile = [all_tile; bund_rsq]; %#ok<AGROW>
                    leg_h(bund) = plot(ax, distances, mean_rsq, '-', ...
                        'Color', bund_col, 'LineWidth', pub_line_width + 0.5, ...
                        'Marker', 'o', 'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', bund_col, 'MarkerEdgeColor', bund_col);
                end

                add_ref_lines(ax);
                xlim(ax, [distances(1), distances(end)]);
                y_lo = max(0,    floor(min(all_tile(:))*100)/100 - 0.01);
                y_hi = min(1.05, ceil(max(all_tile(:))*100)/100 + 0.01);
                ylim(ax, [y_lo, y_hi]);
                xticks(ax, 0:20:ceil(distances(end)));
                grid(ax, 'on');
                set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');
                title(ax, orientation_display{ori_idx}, 'FontSize', 12);
                if ori_idx == 1
                    ylabel(ax, 'r²', 'FontSize', 11);
                end
                if ori_idx == n_ori
                    lgd = legend(ax, leg_h, source_bundle_display, ...
                        'Location', 'eastoutside', 'FontSize', 9);
                    lgd.Box = 'off';
                    title(lgd, 'Bundle');
                end
                hold(ax, 'off');
            end

            fname = sprintf('source_summary_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('  Saved: %s\n', fname);
        end
        fprintf('  [%s] Source figures done\n', method);
    end

    % ----------------------------------------------------------------
    % SOURCE CROSS-MODEL — per sensor axis
    % Layout: 1 row × n_ori cols; colour = method, shade = bundle
    % ----------------------------------------------------------------
    bundle_markers = {'o', 's', '^'};

    for sens_ax = 1:n_axes
        fig = figure('Color', 'w', 'Position', [50, 50, 1400, 420]);
        tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
        title(tl, sprintf('Source shift — cross-method  |  Sensor axis %d of %d', ...
            sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
        xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);

        for ori_idx = 1:n_ori
            ori_label = orientation_labels{ori_idx};
            ax = nexttile(tl, ori_idx);
            hold(ax, 'on');

            all_tile = [];
            leg_h   = gobjects(n_loaded_methods * n_src_bundles, 1);
            leg_lbl = cell(n_loaded_methods * n_src_bundles, 1);
            leg_filled = false(n_loaded_methods * n_src_bundles, 1);

            for m_idx = 1:n_loaded_methods
                method    = loaded_methods{m_idx};
                rsq_store = rsq_by_method.(method);
                mcol      = fwd_method_colors(m_idx, :);
                mlabel    = method_label_map(method);

                for bund = 1:n_src_bundles
                    bund_rows = find(valid_source_bundle_idx == bund);
                    bund_rsq  = zeros(numel(bund_rows), numel(distances));
                    for i = 1:numel(bund_rows)
                        bund_rsq(i,:) = squeeze(rsq_store.(ori_label)(bund_rows(i), :, sens_ax));
                    end
                    mean_rsq = mean(bund_rsq, 1, 'omitnan');
                    all_tile = [all_tile; bund_rsq]; %#ok<AGROW>

                    li = (m_idx-1)*n_src_bundles + bund;
                    h  = plot(ax, distances, mean_rsq, '-', ...
                        'Color', mcol, 'LineWidth', pub_line_width, ...
                        'Marker', bundle_markers{bund}, ...
                        'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w');
                    if ~leg_filled(li)
                        leg_h(li)      = h;
                        leg_lbl{li}    = sprintf('%s — %s', mlabel, source_bundle_display{bund});
                        leg_filled(li) = true;
                    end
                end
            end

            add_ref_lines(ax);
            xlim(ax, [distances(1), distances(end)]);
            if ~isempty(all_tile)
                y_lo = max(0,    floor(min(all_tile(:))*100)/100 - 0.01);
                y_hi = min(1.05, ceil(max(all_tile(:))*100)/100 + 0.01);
                ylim(ax, [y_lo, y_hi]);
            end
            xticks(ax, 0:20:ceil(distances(end)));
            grid(ax, 'on');
            set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');
            title(ax, orientation_display{ori_idx}, 'FontSize', 12);
            if ori_idx == 1; ylabel(ax, 'r²', 'FontSize', 11); end
            if ori_idx == n_ori
                lgd = legend(ax, leg_h(leg_filled), leg_lbl(leg_filled), ...
                    'Location', 'eastoutside', 'FontSize', 9);
                lgd.Box = 'off';
            end
            hold(ax, 'off');
        end

        fname = sprintf('source_crossmodel_sensorax%d', sens_ax);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end

    end
    fprintf('Source figures complete.\n\n');
end


%% =========================================================================
%% SENSOR MODE
%% =========================================================================

if run_sensor
    fprintf('SENSOR PERTURBATION FIGURES\n');

    sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
    if ~isfile(sen_file)
        error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
    end
    load(sen_file);   %#ok<LOAD>

    n_loaded_methods = numel(loaded_methods);
    if n_loaded_methods == 0
        warning('No method results in pert_sensor_rsq.mat — skipping sensor figures.');
    else

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    method_label_map = containers.Map(fwd_methods, fwd_method_labels);
    method_style_map = containers.Map(fwd_methods, fwd_method_styles(1:numel(fwd_methods)));

    % ----------------------------------------------------------------
    % SENSOR DETAIL — per method, per sensor axis
    % ----------------------------------------------------------------
    fprintf('  Sensor detail figures (per method)...\n');
    for m_idx = 1:n_loaded_methods
        method    = loaded_methods{m_idx};
        rsq_store = rsq_by_method.(method);
        mlabel    = method_label_map(method);

        for sens_ax = 1:n_axes
            fig = figure('Color', 'w', 'Position', [50, 50, 1400, 960]);
            tl  = tiledlayout(n_sensor_bundles, n_ori, ...
                'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('[%s]  Sensor shift — all realisations  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
            xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);
            ylabel(tl, 'r²  (shifted vs original)', 'FontSize', 12);

            for b = 1:n_sensor_bundles
                bcolor      = sensor_bundle_colors(b, :);
                bundle_mask = valid_bundle_idx == b;
                bund_rows   = find(bundle_mask);
                n_in_bundle = numel(bund_rows);

                for ori_idx = 1:n_ori
                    ori_label = orientation_labels{ori_idx};

                    ax = nexttile(tl, (b-1)*n_ori + ori_idx);
                    hold(ax, 'on');
                    leg_handles = gobjects(n_in_bundle, 1);

                    % Pre-collect tile data for y-axis scaling
                    tile_rsq = zeros(n_in_bundle, numel(distances));
                    for i = 1:n_in_bundle
                        tile_rsq(i,:) = squeeze(rsq_store.(ori_label)(bund_rows(i), :, sens_ax));
                    end

                    for i = 1:n_in_bundle
                        rsq_row = bund_rows(i);
                        t       = (i-1) / max(n_in_bundle-1, 1);
                        col     = bcolor * (0.5 + 0.5*(1-t)) + (1 - (0.5 + 0.5*(1-t))) * [1 1 1];
                        col     = min(1, max(0, col));
                        leg_handles(i) = plot(ax, distances, tile_rsq(i,:), ...
                            '-', 'Color', col, 'LineWidth', pub_line_width - 0.5, ...
                            'Marker', 'o', 'MarkerIndices', marker_idx, ...
                            'MarkerSize', pub_marker_size - 1, ...
                            'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                    end

                    add_ref_lines(ax);
                    xlim(ax, [distances(1), distances(end)]);
                    y_lo = max(0,   floor(min(tile_rsq(:)) * 100) / 100 - 0.01);
                    y_hi = min(1.05, ceil(max(tile_rsq(:)) * 100) / 100 + 0.01);
                    ylim(ax, [y_lo, y_hi]);
                    xticks(ax, 0:20:ceil(distances(end)));
                    grid(ax, 'on');
                    set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');

                    if ori_idx == 1
                        ylabel(ax, sensor_bundle_display{b}, 'FontSize', 11, ...
                            'Color', bcolor);
                    end
                    if b == 1
                        title(ax, orientation_display{ori_idx}, 'FontSize', 12);
                    end
                    if ori_idx == n_ori
                        shift_labels = valid_labels(bundle_mask);
                        lgd = legend(ax, leg_handles, shift_labels(1:n_in_bundle), ...
                            'Location', 'eastoutside', 'FontSize', 9);
                        lgd.Box = 'off';
                        title(lgd, 'Shift');
                    end
                    hold(ax, 'off');
                end
            end

            fname = sprintf('sensor_detail_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
        fprintf('  [%s] Sensor detail done\n', method);
    end

    % ----------------------------------------------------------------
    % SENSOR SUMMARY — per method, per sensor axis
    % ----------------------------------------------------------------
    fprintf('  Sensor summary figures (per method)...\n');
    for m_idx = 1:n_loaded_methods
        method    = loaded_methods{m_idx};
        rsq_store = rsq_by_method.(method);
        mlabel    = method_label_map(method);

        for sens_ax = 1:n_axes
            fig = figure('Color', 'w', 'Position', [50, 50, 1400, 420]);
            tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('[%s]  Sensor shift — bundle means  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
            xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);
            ylabel(tl, 'r²  (shifted vs original)', 'FontSize', 12);

            for ori_idx = 1:n_ori
                ori_label = orientation_labels{ori_idx};

                ax = nexttile(tl);
                hold(ax, 'on');
                leg_handles = gobjects(n_sensor_bundles, 1);
                tile_rsq    = zeros(n_sensor_bundles, numel(distances));

                for b = 1:n_sensor_bundles
                    bcolor      = sensor_bundle_colors(b, :);
                    bundle_mask = valid_bundle_idx == b;
                    bund_rows   = find(bundle_mask);

                    rsq_all  = squeeze(rsq_store.(ori_label)(bund_rows, :, sens_ax));
                    if numel(bund_rows) > 1
                        rsq_mean = mean(rsq_all, 1);
                    else
                        rsq_mean = rsq_all;
                    end
                    tile_rsq(b,:) = rsq_mean;

                    leg_handles(b) = plot(ax, distances, rsq_mean, ...
                        '-', 'Color', bcolor, 'LineWidth', pub_line_width + 0.5, ...
                        'Marker', 'o', 'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', bcolor, 'MarkerEdgeColor', bcolor);
                end

                add_ref_lines(ax);
                xlim(ax, [distances(1), distances(end)]);
                y_lo = max(0,   floor(min(tile_rsq(:)) * 100) / 100 - 0.01);
                y_hi = min(1.05, ceil(max(tile_rsq(:)) * 100) / 100 + 0.01);
                ylim(ax, [y_lo, y_hi]);
                xticks(ax, 0:20:ceil(distances(end)));
                grid(ax, 'on');
                set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');
                title(ax, orientation_display{ori_idx}, 'FontSize', 12);

                if ori_idx == n_ori
                    lgd = legend(ax, leg_handles, sensor_bundle_display, ...
                        'Location', 'eastoutside', 'FontSize', 10);
                    lgd.Box = 'off';
                    title(lgd, 'Bundle');
                end
                hold(ax, 'off');
            end

            fname = sprintf('sensor_summary_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
        fprintf('  [%s] Sensor summary done\n', method);
    end

    % ----------------------------------------------------------------
    % SENSOR CROSS-MODEL — per sensor axis (only when >1 method loaded)
    % ----------------------------------------------------------------
    if n_loaded_methods > 1
        fprintf('  Sensor cross-model figures...\n');
        for sens_ax = 1:n_axes
            fig = figure('Color', 'w', 'Position', [50, 50, 1400, 420]);
            tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('Sensor shift — cross-model bundle means  |  Sensor axis %d of %d', ...
                sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
            xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);
            ylabel(tl, 'r²  (shifted vs original)', 'FontSize', 12);

            % Colour = method, marker = bundle
            bundle_markers  = {'o', 's', '^'};   % circle / square / triangle
            n_lines     = n_loaded_methods * n_sensor_bundles;
            leg_handles = gobjects(n_lines, 1);
            leg_labels  = cell(n_lines, 1);

            for ori_idx = 1:n_ori
                ori_label = orientation_labels{ori_idx};

                ax = nexttile(tl);
                hold(ax, 'on');
                line_n   = 0;
                tile_rsq = zeros(n_lines, numel(distances));

                for m_idx = 1:n_loaded_methods
                    method    = loaded_methods{m_idx};
                    rsq_store = rsq_by_method.(method);
                    mlabel    = method_label_map(method);
                    mcolor    = fwd_method_colors(m_idx, :);

                    for b = 1:n_sensor_bundles
                        bundle_mask = valid_bundle_idx == b;
                        bund_rows   = find(bundle_mask);
                        bmarker     = bundle_markers{b};

                        rsq_all  = squeeze(rsq_store.(ori_label)(bund_rows, :, sens_ax));
                        if numel(bund_rows) > 1
                            rsq_mean = mean(rsq_all, 1);
                        else
                            rsq_mean = rsq_all;
                        end

                        line_n = line_n + 1;
                        tile_rsq(line_n,:) = rsq_mean;
                        h = plot(ax, distances, rsq_mean, ...
                            'LineStyle', '-', 'Color', mcolor, ...
                            'LineWidth', pub_line_width, ...
                            'Marker', bmarker, 'MarkerIndices', marker_idx, ...
                            'MarkerSize', pub_marker_size - 1, ...
                            'MarkerFaceColor', mcolor, 'MarkerEdgeColor', mcolor);
                        if ori_idx == 1
                            leg_handles(line_n) = h;
                            leg_labels{line_n}  = sprintf('%s — %s', ...
                                mlabel, sensor_bundle_display{b});
                        end
                    end
                end

                add_ref_lines(ax);
                xlim(ax, [distances(1), distances(end)]);
                y_lo = max(0,   floor(min(tile_rsq(:)) * 100) / 100 - 0.01);
                y_hi = min(1.05, ceil(max(tile_rsq(:)) * 100) / 100 + 0.01);
                ylim(ax, [y_lo, y_hi]);
                xticks(ax, 0:20:ceil(distances(end)));
                grid(ax, 'on');
                set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');
                title(ax, orientation_display{ori_idx}, 'FontSize', 12);

                if ori_idx == n_ori
                    lgd = legend(ax, leg_handles(1:n_lines), leg_labels(1:n_lines), ...
                        'Location', 'eastoutside', 'FontSize', 9);
                    lgd.Box = 'off';
                    title(lgd, 'Method — Bundle');
                end
                hold(ax, 'off');
            end

            fname = sprintf('sensor_crossmodel_sensorax%d', sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
        fprintf('  Sensor cross-model done\n');
    end

    end
    fprintf('Sensor figures complete.\n\n');
end

%% =========================================================================
%% CONDUCTIVITY MODE
%% =========================================================================

if run_cond
    fprintf('CONDUCTIVITY PERTURBATION FIGURES\n');

    cond_file = fullfile(forward_fields_base, 'pert_cond_rsq.mat');
    if ~isfile(cond_file)
        error('Cond r² file not found: %s\nRun pt_compute_rsq first.', cond_file);
    end
    load(cond_file);   %#ok<LOAD>

    if isempty(valid_cond_keys)
        warning('No conductivity perturbation results found — skipping cond figures.');
    else

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'cond');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    rsq_store     = rsq_by_method.bem;
    n_cond_bund   = numel(cond_bundle_display);

    % ----------------------------------------------------------------
    % COND DETAIL — per sensor axis
    % Layout: n_cond_bund rows × n_ori cols; 8 shaded individual lines per tile
    % ----------------------------------------------------------------
    fprintf('  Conductivity detail figures...\n');
    for sens_ax = 1:n_axes
        fig = figure('Color', 'w', 'Position', [50, 50, 1400, 960]);
        tl  = tiledlayout(n_cond_bund, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
        title(tl, sprintf('[BEM]  Conductivity perturbation — all realisations  |  Sensor axis %d of %d', ...
            sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
        xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);
        ylabel(tl, 'r²  (perturbed vs nominal BEM)', 'FontSize', 12);

        for bund = 1:n_cond_bund
            bund_mask   = valid_cond_bundle_idx == bund;
            bund_rows   = find(bund_mask);
            n_in_bundle = numel(bund_rows);
            bcolor      = cond_bundle_colors(bund, :);

            for ori_idx = 1:n_ori
                ori_label = orientation_labels{ori_idx};
                ax = nexttile(tl, (bund-1)*n_ori + ori_idx);
                hold(ax, 'on');
                leg_handles = gobjects(n_in_bundle, 1);

                tile_rsq = zeros(n_in_bundle, numel(distances));
                for i = 1:n_in_bundle
                    tile_rsq(i,:) = squeeze(rsq_store.(ori_label)(bund_rows(i), :, sens_ax));
                end

                for i = 1:n_in_bundle
                    t   = (i-1) / max(n_in_bundle-1, 1);
                    col = bcolor * (0.5 + 0.5*(1-t)) + (1 - (0.5 + 0.5*(1-t))) * [1 1 1];
                    col = min(1, max(0, col));
                    leg_handles(i) = plot(ax, distances, tile_rsq(i,:), ...
                        '-', 'Color', col, 'LineWidth', pub_line_width - 0.5, ...
                        'Marker', 'o', 'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size - 1, ...
                        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                end

                add_ref_lines(ax);
                xlim(ax, [distances(1), distances(end)]);
                y_lo = max(0,    floor(min(tile_rsq(:))*100)/100 - 0.01);
                y_hi = min(1.05, ceil(max(tile_rsq(:))*100)/100 + 0.01);
                ylim(ax, [y_lo, y_hi]);
                xticks(ax, 0:20:ceil(distances(end)));
                grid(ax, 'on');
                set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');

                if ori_idx == 1
                    ylabel(ax, {cond_bundle_display{bund}; 'r²'}, ...
                        'FontSize', 11, 'Color', bcolor);
                end
                if bund == 1
                    title(ax, orientation_display{ori_idx}, 'FontSize', 12);
                end
                if ori_idx == n_ori
                    shift_labels = valid_cond_labels(bund_mask);
                    lgd = legend(ax, leg_handles, shift_labels(1:n_in_bundle), ...
                        'Location', 'eastoutside', 'FontSize', 9);
                    lgd.Box = 'off';
                    title(lgd, 'Shift');
                end
                hold(ax, 'off');
            end
        end

        fname = sprintf('cond_detail_sensorax%d', sens_ax);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end

    % ----------------------------------------------------------------
    % COND SUMMARY — per sensor axis
    % Layout: 1 row × n_ori cols; one mean line per bundle with markers
    % ----------------------------------------------------------------
    fprintf('  Conductivity summary figures...\n');
    for sens_ax = 1:n_axes
        fig = figure('Color', 'w', 'Position', [50, 50, 1400, 420]);
        tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
        title(tl, sprintf('[BEM]  Conductivity perturbation — bundle means  |  Sensor axis %d of %d', ...
            sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
        xlabel(tl, 'Distance along spinal cord (mm)', 'FontSize', 12);
        ylabel(tl, 'r²  (perturbed vs nominal BEM)', 'FontSize', 12);

        for ori_idx = 1:n_ori
            ori_label = orientation_labels{ori_idx};
            ax = nexttile(tl);
            hold(ax, 'on');

            all_tile    = [];
            leg_handles = gobjects(n_cond_bund, 1);

            for bund = 1:n_cond_bund
                bund_rows = find(valid_cond_bundle_idx == bund);
                bcolor    = cond_bundle_colors(bund, :);

                rsq_all = squeeze(rsq_store.(ori_label)(bund_rows, :, sens_ax));
                if numel(bund_rows) > 1
                    rsq_mean = mean(rsq_all, 1);
                else
                    rsq_mean = rsq_all;
                end
                all_tile = [all_tile; rsq_all]; %#ok<AGROW>

                leg_handles(bund) = plot(ax, distances, rsq_mean, '-', ...
                    'Color', bcolor, 'LineWidth', pub_line_width + 0.5, ...
                    'Marker', 'o', 'MarkerIndices', marker_idx, ...
                    'MarkerSize', pub_marker_size, ...
                    'MarkerFaceColor', bcolor, 'MarkerEdgeColor', bcolor);
            end

            add_ref_lines(ax);
            xlim(ax, [distances(1), distances(end)]);
            y_lo = max(0,    floor(min(all_tile(:))*100)/100 - 0.01);
            y_hi = min(1.05, ceil(max(all_tile(:))*100)/100 + 0.01);
            ylim(ax, [y_lo, y_hi]);
            xticks(ax, 0:20:ceil(distances(end)));
            grid(ax, 'on');
            set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');
            title(ax, orientation_display{ori_idx}, 'FontSize', 12);
            if ori_idx == 1; ylabel(ax, 'r²', 'FontSize', 11); end
            if ori_idx == n_ori
                lgd = legend(ax, leg_handles, cond_bundle_display, ...
                    'Location', 'eastoutside', 'FontSize', 9);
                lgd.Box = 'off';
                title(lgd, 'Bundle');
            end
            hold(ax, 'off');
        end

        fname = sprintf('cond_summary_sensorax%d', sens_ax);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end

    end
    fprintf('Conductivity figures complete.\n\n');
end

fprintf('pt_plot_curves complete.\n');


% ---- Local functions ----

function add_ref_lines(ax)
    yline(ax, 1.00, '--k',  'LineWidth', 1.0, 'Alpha', 0.45, ...
        'Label', 'r²=1.00', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
    yline(ax, 0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.45, 'Color', [0.4 0.4 0.4], ...
        'Label', 'r²=0.99', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
    yline(ax, 0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.45, 'Color', [0.6 0.6 0.6], ...
        'Label', 'r²=0.95', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
end
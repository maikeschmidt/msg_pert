% pt_plot_curves - Plot r² vs distance along spinal cord for perturbation analysis
%
% Loads pre-computed r² data (from pt_compute_rsq) and produces curve figures
% showing how leadfield similarity degrades along the cord for each perturbation.
%
% SOURCE mode: grouped by shift axis (X/Y/Z), individual + overview figures
% SENSOR mode: grouped by bundle (~2mm/~5mm/~10mm), individual + overview figures
%
% Figures saved to <save_base_dir>/perturbation_analysis/source/ and .../sensor/
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

config_pert;
pt_add_functions;

%% SOURCE MODE

if run_source
    fprintf('SOURCE PERTURBATION CURVE FIGURES\n');

    src_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
    if ~isfile(src_file)
        error('Source r² file not found: %s\nRun pt_compute_rsq first.', src_file);
    end
    load(src_file);

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'source');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    shift_axis_short  = {'X (Left-Right)', 'Y (Rostral-Caudal)', 'Z (Ventral-Dorsal)'};
    shift_axis_labels = {'X', 'Y', 'Z'};

    % Individual figures — one per shift axis per orientation per sensor axis
    fprintf('  Generating individual figures...\n');
    for shift_ax = 1:3
        axis_mask  = valid_shift_axis == shift_ax;
        ax_row_idx = find(axis_mask);
        n_ax_s     = numel(ax_row_idx);
        base_col   = sensitivity_axis_colors(shift_ax, :);
        if n_ax_s == 0; continue; end

        for sens_ax = 1:n_axes
            for ori_idx = 1:numel(orientation_labels)
                ori_label = orientation_labels{ori_idx};

                fig = figure('Color', 'w', 'Position', [100, 100, 1100, 650]);
                hold on;
                leg_h = gobjects(n_ax_s, 1);

                for i = 1:n_ax_s
                    mag_scale  = 1 - (mod(i-1, 3) * 0.35);
                    col_scaled = min(1, base_col + (1-base_col) * (1-mag_scale));
                    rsq_row    = ax_row_idx(i);

                    leg_h(i) = plot(distances, ...
                        squeeze(rsq_store.(ori_label)(rsq_row, :, sens_ax)), ...
                        'LineStyle', valid_styles{rsq_row}, ...
                        'Color', col_scaled, 'LineWidth', pub_line_width, ...
                        'Marker', valid_markers{rsq_row}, ...
                        'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col_scaled, 'MarkerEdgeColor', col_scaled);
                end

                yline(1.00, '--k', 'LineWidth', 1.2, 'Alpha', 0.5, 'Label', 'r²=1.00', ...
                    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
                yline(0.99, ':', 'LineWidth', 1.2, 'Alpha', 0.5, 'Color', [0.4 0.4 0.4], ...
                    'Label', 'r²=0.99', 'LabelHorizontalAlignment', 'left', 'FontSize', 10);
                yline(0.95, ':', 'LineWidth', 1.2, 'Alpha', 0.5, 'Color', [0.6 0.6 0.6], ...
                    'Label', 'r²=0.95', 'LabelHorizontalAlignment', 'left', 'FontSize', 10);

                xlim([distances(1), distances(end)]);
                xticks(0:20:ceil(distances(end)));
                ylim([0, 1.05]);

                title(sprintf('Source shift along %s — %s  |  Sensor axis %d of %d', ...
                    shift_axis_short{shift_ax}, orientation_display{ori_idx}, ...
                    sens_ax, n_axes), 'FontSize', 14, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
                ylabel('r²  (shifted vs original leadfield)', 'FontSize', 13);

                lgd     = legend(leg_h, {'+2mm','+4mm','+6mm','-2mm','-4mm','-6mm'}, ...
                    'Location', 'eastoutside', 'FontSize', 12);
                lgd.Box = 'off';
                title(lgd, sprintf('Shift\n(%s axis)', shift_axis_labels{shift_ax}));

                grid on;
                set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');

                fname = sprintf('source_%sshift_sensorax%d_%s', ...
                    shift_axis_labels{shift_ax}, sens_ax, ori_label);
                exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
                saveas(fig, fullfile(save_dir, [fname '.fig']));
                close(fig);
                fprintf('    Saved: %s\n', fname);
            end
        end
    end

    % Overview figures — all three shift axes side by side
    fprintf('  Generating overview figures...\n');
    for sens_ax = 1:n_axes
        for ori_idx = 1:numel(orientation_labels)
            ori_label = orientation_labels{ori_idx};

            fig = figure('Color', 'w', 'Position', [100, 100, 1900, 650]);
            tl  = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
            title(tl, sprintf('Source Shift Perturbation — %s  |  Sensor axis %d of %d', ...
                orientation_display{ori_idx}, sens_ax, n_axes), ...
                'FontSize', 13, 'FontWeight', 'bold');

            for shift_ax = 1:3
                axis_mask  = valid_shift_axis == shift_ax;
                ax_row_idx = find(axis_mask);
                n_ax_s     = numel(ax_row_idx);
                base_col   = sensitivity_axis_colors(shift_ax, :);

                nexttile(tl);
                hold on;
                leg_h = gobjects(n_ax_s, 1);

                for i = 1:n_ax_s
                    mag_scale  = 1 - (mod(i-1, 3) * 0.35);
                    col_scaled = min(1, base_col + (1-base_col) * (1-mag_scale));
                    rsq_row    = ax_row_idx(i);

                    leg_h(i) = plot(distances, ...
                        squeeze(rsq_store.(ori_label)(rsq_row, :, sens_ax)), ...
                        'LineStyle', valid_styles{rsq_row}, ...
                        'Color', col_scaled, 'LineWidth', pub_line_width, ...
                        'Marker', valid_markers{rsq_row}, ...
                        'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col_scaled, 'MarkerEdgeColor', col_scaled);
                end

                yline(1.00, '--k', 'LineWidth', 1.0, 'Alpha', 0.5);
                yline(0.99, ':', 'LineWidth', 1.0, 'Alpha', 0.5, 'Color', [0.4 0.4 0.4]);
                yline(0.95, ':', 'LineWidth', 1.0, 'Alpha', 0.5, 'Color', [0.6 0.6 0.6]);

                xlim([distances(1), distances(end)]);
                xticks(0:20:ceil(distances(end)));
                ylim([0, 1.05]);

                title(sprintf('%s axis shift\n(±2, ±4, ±6 mm)', ...
                    shift_axis_short{shift_ax}), 'FontSize', 13, 'FontWeight', 'bold');
                xlabel('Distance along cord (mm)', 'FontSize', 12);
                if shift_ax == 1
                    ylabel({'r² (shifted vs original)'; '1.0 = no effect'}, 'FontSize', 12);
                end

                lgd     = legend(leg_h, {'+2mm','+4mm','+6mm','-2mm','-4mm','-6mm'}, ...
                    'Location', 'eastoutside', 'FontSize', 10);
                lgd.Box = 'off';
                grid on;
                set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickDir', 'out');
                hold off;
            end

            fname = sprintf('source_overview_sensorax%d_%s', sens_ax, ori_label);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end
    fprintf('Source curve figures complete.\n\n');
end


%% SENSOR MODE

if run_sensor
    fprintf('SENSOR PERTURBATION CURVE FIGURES\n');

    sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
    if ~isfile(sen_file)
        error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
    end
    load(sen_file);

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    fprintf('  Generating bundle figures...\n');
    for b = 1:n_sensor_bundles
        bundle_mask = valid_bundle_idx == b;
        bund_rows   = find(bundle_mask);
        n_in_bundle = numel(bund_rows);
        base_col    = sensor_bundle_colors(b, :);
        if n_in_bundle == 0; continue; end

        shift_colors = zeros(n_in_bundle, 3);
        for i = 1:n_in_bundle
            t = (i-1) / max(n_in_bundle-1, 1);
            shift_colors(i,:) = min(1, base_col + (1-base_col) * t * 0.6);
        end

        for sens_ax = 1:n_axes
            for ori_idx = 1:numel(orientation_labels)
                ori_label = orientation_labels{ori_idx};

                fig = figure('Color', 'w', 'Position', [100, 100, 1100, 650]);
                hold on;
                leg_h = gobjects(n_in_bundle, 1);

                for i = 1:n_in_bundle
                    rsq_row = bund_rows(i);
                    col     = shift_colors(i, :);
                    leg_h(i) = plot(distances, ...
                        squeeze(rsq_store.(ori_label)(rsq_row, :, sens_ax)), ...
                        '-', 'Color', col, 'LineWidth', pub_line_width, ...
                        'Marker', 'o', 'MarkerIndices', marker_idx, ...
                        'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                end

                yline(1.00, '--k', 'LineWidth', 1.2, 'Alpha', 0.5, 'Label', 'r²=1.00', ...
                    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
                yline(0.99, ':', 'LineWidth', 1.2, 'Alpha', 0.5, 'Color', [0.4 0.4 0.4], ...
                    'Label', 'r²=0.99', 'LabelHorizontalAlignment', 'left', 'FontSize', 10);
                yline(0.95, ':', 'LineWidth', 1.2, 'Alpha', 0.5, 'Color', [0.6 0.6 0.6], ...
                    'Label', 'r²=0.95', 'LabelHorizontalAlignment', 'left', 'FontSize', 10);

                xlim([distances(1), distances(end)]); ylim([0, 1.05]);
                xticks(0:20:ceil(distances(end)));

                title(sprintf('Sensor shift — %s  |  %s  |  Sensor axis %d', ...
                    sensor_bundle_display{b}, orientation_display{ori_idx}, sens_ax), ...
                    'FontSize', 14, 'FontWeight', 'bold');
                xlabel('Distance along spinal cord (mm)', 'FontSize', 14);
                ylabel('r²  (shifted vs original)', 'FontSize', 13);

                lgd     = legend(leg_h, valid_labels(bundle_mask), ...
                    'Location', 'eastoutside', 'FontSize', 12);
                lgd.Box = 'off';
                title(lgd, sensor_bundle_display{b});
                grid on;
                set(gca, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out');

                fname = sprintf('sensor_bundle%d_sensorax%d_%s', b, sens_ax, ori_label);
                exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
                saveas(fig, fullfile(save_dir, [fname '.fig']));
                close(fig);
                fprintf('    Saved: %s\n', fname);
            end
        end
    end
    fprintf('Sensor curve figures complete.\n\n');
end

fprintf('pt_plot_curves complete.\n');

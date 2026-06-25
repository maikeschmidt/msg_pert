% pt_plot_heatmaps - Heatmap summaries for perturbation analysis
%
% Loads pre-computed r² data (from pt_compute_rsq) and produces heatmap
% figures that summarise perturbation sensitivity across shift types,
% orientations, and forward model methods.
%
% FIGURES PRODUCED
%   Source — within-model (one per method, per sensor axis):
%     source_heatmap_within_<method>_sensorax<N>.png
%       18-shift × 3-orientation heatmap; rows grouped by shift axis.
%       Cell = median r² across source positions.
%
%   Source — between-models (one per sensor axis):
%     source_heatmap_between_sensorax<N>.png
%       methods × 9-shift-groups (axis × magnitude, ± averaged).
%       Cell = mean over orientations of median r².
%
%   Sensor — within-model (one per method, per sensor axis):
%     sensor_heatmap_within_<method>_sensorax<N>.png
%       3-bundle × 3-orientation heatmap.
%       Cell = mean over bundle shifts of median r² across sources.
%
%   Sensor — between-models (one per sensor axis):
%     sensor_heatmap_between_sensorax<N>.png
%       methods × 3-bundles heatmap.
%       Cell = mean over orientations and shifts of median r².
%
% USAGE:
%   pt_plot_heatmaps
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

n_ori = numel(orientation_labels);


%% =========================================================================
%% SOURCE HEATMAPS
%% =========================================================================

if run_source
    fprintf('SOURCE HEATMAPS\n');

    src_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
    if ~isfile(src_file)
        error('Source r² file not found: %s\nRun pt_compute_rsq first.', src_file);
    end
    load(src_file);   %#ok<LOAD>

    n_loaded = numel(loaded_methods);
    if n_loaded == 0
        warning('No method results in pert_source_rsq.mat — skipping source heatmaps.');
    else

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'source');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    method_label_map = containers.Map(fwd_methods, fwd_method_labels);
    n_shifts = numel(valid_labels);

    % Extract magnitude (mm) per shift row from label (e.g. "X+2mm" -> 2)
    mags_per_row = zeros(1, n_shifts);
    for i = 1:n_shifts
        tok = regexp(valid_labels{i}, '[+-](\d+)mm', 'tokens');
        if ~isempty(tok)
            mags_per_row(i) = str2double(tok{1}{1});
        end
    end

    % ----------------------------------------------------------------
    % SOURCE WITHIN-MODEL: 18-shift × 3-orientation, one per method
    % ----------------------------------------------------------------
    fprintf('  Within-model source heatmaps...\n');
    for m_idx = 1:n_loaded
        method    = loaded_methods{m_idx};
        rsq_store = rsq_by_method.(method);
        mlabel    = method_label_map(method);

        for sens_ax = 1:n_axes
            % Build [n_shifts x n_ori] matrix: median r² over source positions
            hmap = zeros(n_shifts, n_ori);
            for i = 1:n_shifts
                for j = 1:n_ori
                    ori_label  = orientation_labels{j};
                    hmap(i, j) = median(squeeze(rsq_store.(ori_label)(i, :, sens_ax)));
                end
            end

            fig_h = max(450, n_shifts * 26 + 180);
            fig   = figure('Color', 'w', 'Position', [50, 50, 720, fig_h]);
            hax   = axes(fig);   %#ok<LAXES>

            imagesc(hax, hmap);
            colormap(hax, flipud(cool));
            clim(hax, [max(0, floor(min(hmap(:)) * 100)/100 - 0.01), 1]);
            cb = colorbar(hax);
            cb.Label.String   = 'median r²';
            cb.Label.FontSize = 11;

            xticks(hax, 1:n_ori);
            xticklabels(hax, orientation_display);
            yticks(hax, 1:n_shifts);
            yticklabels(hax, valid_labels);
            set(hax, 'FontSize', 11, 'TickDir', 'out', 'YDir', 'normal');
            title(hax, sprintf('[%s]  Source shifts — median r²  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');

            % Annotate cells
            clo = clim(hax); crange = clo(2) - clo(1);
            for i = 1:n_shifts
                for j = 1:n_ori
                    val  = hmap(i, j);
                    tcol = pick_text_color(val, clo(1), crange);
                    text(hax, j, i, sprintf('%.3f', val), ...
                        'HorizontalAlignment', 'center', 'FontSize', 9, ...
                        'FontWeight', 'bold', 'Color', tcol);
                end
            end

            % Separator lines between shift axes (after last row of each axis)
            for sh_ax = 1:2
                boundary = find(valid_shift_axis == sh_ax, 1, 'last') + 0.5;
                yline(hax, boundary, 'k-', 'LineWidth', 1.8);
            end

            fname = sprintf('source_heatmap_within_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end

    % ----------------------------------------------------------------
    % SOURCE BETWEEN-MODELS: methods × 9 shift groups (axis × magnitude, ± averaged)
    % ----------------------------------------------------------------
    fprintf('  Between-models source heatmaps...\n');

    group_labels = cell(1, 9);
    col = 0;
    for sh_ax = 1:3
        for mag = [2, 4, 6]
            col = col + 1;
            group_labels{col} = sprintf('%s ±%dmm', source_shift_axes{sh_ax}, mag);
        end
    end
    n_groups = 9;

    for sens_ax = 1:n_axes
        hmap       = zeros(n_loaded, n_groups);
        row_labels = cell(1, n_loaded);

        for m_idx = 1:n_loaded
            method    = loaded_methods{m_idx};
            rsq_store = rsq_by_method.(method);
            row_labels{m_idx} = method_label_map(method);

            col = 0;
            for sh_ax = 1:3
                for mag = [2, 4, 6]
                    col  = col + 1;
                    rows = find(valid_shift_axis == sh_ax & mags_per_row == mag);
                    vals = [];
                    for j = 1:n_ori
                        ori_label = orientation_labels{j};
                        for r = rows(:)'
                            vals(end+1) = median(squeeze( ...
                                rsq_store.(ori_label)(r, :, sens_ax))); %#ok<AGROW>
                        end
                    end
                    hmap(m_idx, col) = mean(vals);
                end
            end
        end

        fig_w = max(600, n_groups * 70 + 220);
        fig_h = max(300, n_loaded  * 65 + 180);
        fig   = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);
        hax   = axes(fig);   %#ok<LAXES>

        imagesc(hax, hmap);
        colormap(hax, flipud(cool));
        clim(hax, [max(0, floor(min(hmap(:)) * 100)/100 - 0.01), 1]);
        cb = colorbar(hax);
        cb.Label.String   = 'mean median r² (± averaged)';
        cb.Label.FontSize = 11;

        xticks(hax, 1:n_groups); xticklabels(hax, group_labels); xtickangle(hax, 40);
        yticks(hax, 1:n_loaded); yticklabels(hax, row_labels);
        set(hax, 'FontSize', 11, 'TickDir', 'out');
        title(hax, sprintf('Source shifts — cross-model comparison  |  Sensor axis %d of %d', ...
            sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');

        % Vertical separators between shift axes (every 3 columns)
        xline(hax, 3.5, 'k-', 'LineWidth', 1.8);
        xline(hax, 6.5, 'k-', 'LineWidth', 1.8);

        clo = clim(hax); crange = clo(2) - clo(1);
        for r = 1:n_loaded
            for c = 1:n_groups
                val  = hmap(r, c);
                tcol = pick_text_color(val, clo(1), crange);
                text(hax, c, r, sprintf('%.3f', val), ...
                    'HorizontalAlignment', 'center', 'FontSize', 10, ...
                    'FontWeight', 'bold', 'Color', tcol);
            end
        end

        fname = sprintf('source_heatmap_between_sensorax%d', sens_ax);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end

    end
    fprintf('Source heatmaps complete.\n\n');
end


%% =========================================================================
%% SENSOR HEATMAPS
%% =========================================================================

if run_sensor
    fprintf('SENSOR HEATMAPS\n');

    sen_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
    if ~isfile(sen_file)
        error('Sensor r² file not found: %s\nRun pt_compute_rsq first.', sen_file);
    end
    load(sen_file);   %#ok<LOAD>

    n_loaded = numel(loaded_methods);
    if n_loaded == 0
        warning('No method results in pert_sensor_rsq.mat — skipping sensor heatmaps.');
    else

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    method_label_map = containers.Map(fwd_methods, fwd_method_labels);

    % ----------------------------------------------------------------
    % SENSOR WITHIN-MODEL: 3-bundle × 3-orientation, one per method
    % ----------------------------------------------------------------
    fprintf('  Within-model sensor heatmaps...\n');
    for m_idx = 1:n_loaded
        method    = loaded_methods{m_idx};
        rsq_store = rsq_by_method.(method);
        mlabel    = method_label_map(method);

        for sens_ax = 1:n_axes
            hmap = zeros(n_sensor_bundles, n_ori);
            for b = 1:n_sensor_bundles
                bund_rows = find(valid_bundle_idx == b);
                for j = 1:n_ori
                    ori_label = orientation_labels{j};
                    vals = zeros(1, numel(bund_rows));
                    for ki = 1:numel(bund_rows)
                        vals(ki) = median(squeeze( ...
                            rsq_store.(ori_label)(bund_rows(ki), :, sens_ax)));
                    end
                    hmap(b, j) = mean(vals);
                end
            end

            fig   = figure('Color', 'w', 'Position', [50, 50, 620, 340]);
            hax   = axes(fig);   %#ok<LAXES>

            imagesc(hax, hmap);
            colormap(hax, flipud(cool));
            clim(hax, [max(0, floor(min(hmap(:)) * 100)/100 - 0.01), 1]);
            cb = colorbar(hax);
            cb.Label.String   = 'mean median r²';
            cb.Label.FontSize = 11;

            xticks(hax, 1:n_ori);
            xticklabels(hax, orientation_display);
            yticks(hax, 1:n_sensor_bundles);
            yticklabels(hax, sensor_bundle_display);
            set(hax, 'FontSize', 11, 'TickDir', 'out', 'YDir', 'normal');
            title(hax, sprintf('[%s]  Sensor shifts — mean median r²  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');

            clo = clim(hax); crange = clo(2) - clo(1);
            for r = 1:n_sensor_bundles
                for c = 1:n_ori
                    val  = hmap(r, c);
                    tcol = pick_text_color(val, clo(1), crange);
                    text(hax, c, r, sprintf('%.3f', val), ...
                        'HorizontalAlignment', 'center', 'FontSize', 12, ...
                        'FontWeight', 'bold', 'Color', tcol);
                end
            end

            fname = sprintf('sensor_heatmap_within_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end

    % ----------------------------------------------------------------
    % SENSOR BETWEEN-MODELS: methods × 3 bundles
    % ----------------------------------------------------------------
    fprintf('  Between-models sensor heatmaps...\n');
    for sens_ax = 1:n_axes
        hmap       = zeros(n_loaded, n_sensor_bundles);
        row_labels = cell(1, n_loaded);

        for m_idx = 1:n_loaded
            method    = loaded_methods{m_idx};
            rsq_store = rsq_by_method.(method);
            row_labels{m_idx} = method_label_map(method);

            for b = 1:n_sensor_bundles
                bund_rows = find(valid_bundle_idx == b);
                vals = [];
                for j = 1:n_ori
                    ori_label = orientation_labels{j};
                    for ki = 1:numel(bund_rows)
                        vals(end+1) = median(squeeze( ...
                            rsq_store.(ori_label)(bund_rows(ki), :, sens_ax))); %#ok<AGROW>
                    end
                end
                hmap(m_idx, b) = mean(vals);
            end
        end

        fig_h = max(300, n_loaded * 65 + 180);
        fig   = figure('Color', 'w', 'Position', [50, 50, 600, fig_h]);
        hax   = axes(fig);   %#ok<LAXES>

        imagesc(hax, hmap);
        colormap(hax, flipud(cool));
        clim(hax, [max(0, floor(min(hmap(:)) * 100)/100 - 0.01), 1]);
        cb = colorbar(hax);
        cb.Label.String   = 'mean median r²';
        cb.Label.FontSize = 11;

        xticks(hax, 1:n_sensor_bundles); xticklabels(hax, sensor_bundle_display);
        yticks(hax, 1:n_loaded); yticklabels(hax, row_labels);
        set(hax, 'FontSize', 11, 'TickDir', 'out');
        title(hax, sprintf('Sensor shifts — cross-model comparison  |  Sensor axis %d of %d', ...
            sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');

        clo = clim(hax); crange = clo(2) - clo(1);
        for r = 1:n_loaded
            for c = 1:n_sensor_bundles
                val  = hmap(r, c);
                tcol = pick_text_color(val, clo(1), crange);
                text(hax, c, r, sprintf('%.3f', val), ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, ...
                    'FontWeight', 'bold', 'Color', tcol);
            end
        end

        fname = sprintf('sensor_heatmap_between_sensorax%d', sens_ax);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('    Saved: %s\n', fname);
    end

    end
    fprintf('Sensor heatmaps complete.\n\n');
end

fprintf('pt_plot_heatmaps complete.\n');


% ---- Local functions ----

function tcol = pick_text_color(val, clo, crange)
% White text on dark (low-r²) cells, black text on light (high-r²) cells.
% The flipud(cool) map is dark (magenta) at low values and light (cyan) at high.
    if crange < 1e-9
        tcol = 'k';
        return
    end
    norm_val = (val - clo) / crange;   % 0 = low r² (dark bg), 1 = high r² (light bg)
    if norm_val < 0.45
        tcol = 'w';
    else
        tcol = 'k';
    end
end

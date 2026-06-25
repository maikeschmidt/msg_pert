% pt_plot_heatmaps - Pairwise RE and r² heatmaps for perturbation analysis
%
% Computes pairwise relative error (RE) and squared Pearson correlation (r²)
% between all shifted model variants (and the original) for each forward
% model method, then displays the results as annotated colour heatmaps.
% Mirrors the approach in msg_fwd/plot_pairwise_heatmaps.m but:
%   (1) separates dipole orientations into individual subfigures instead of
%       concatenating them, and
%   (2) the items being compared are the shifted model variants rather than
%       different forward model methods.
%
% FIGURES PRODUCED
%   Source — within-method (one per method, per sensor axis):
%     source_within_<method>_sensorax<N>.png
%       2 × n_ori tiled layout (RE top, r² bottom).
%       Items: original + 18 source shifts (19 × 19 pairwise).
%       Rows/cols grouped by shift axis (X/Y/Z) with separator lines.
%
%   Source — cross-method (one per shift axis, per sensor axis):
%     source_cross_<axis>_sensorax<N>.png
%       2 × n_ori tiled layout.
%       Items: for each method: original + that axis's 6 shifts.
%       Rows/cols grouped by method with separator lines.
%
%   Sensor — within-method (one per method, per sensor axis):
%     sensor_within_<method>_sensorax<N>.png
%       2 × n_ori tiled layout.
%       Items: original + 24 sensor shifts (25 × 25 pairwise).
%       Rows/cols grouped by bundle with separator lines.
%
%   Sensor — cross-method (one per bundle, per sensor axis):
%     sensor_cross_bundle<B>_sensorax<N>.png
%       2 × n_ori tiled layout.
%       Items: for each method: original + that bundle's 8 shifts.
%       Rows/cols grouped by method with separator lines.
%
% USAGE:
%   pt_plot_heatmaps
%
% DEPENDENCIES:
%   config_pert, leadfields_organised.mat, pert_source_rsq.mat,
%   pert_sensor_rsq.mat, compare_results() (via msg_fwd/functions/)
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

% Load raw leadfields (needed for pairwise comparison)
lf_file = fullfile(forward_fields_base, 'leadfields_organised.mat');
if ~isfile(lf_file)
    error('leadfields_organised.mat not found: %s\nRun pt_load_leadfields first.', lf_file);
end
load(lf_file, 'leadfields');   %#ok<LOAD>

% Method abbreviations for cross-method labels
abbrev_map = containers.Map( ...
    {'bslaw', 'sphere', 'bem', 'fem'}, ...
    {'BS',    'SP',     'BEM', 'FEM'});
method_label_map = containers.Map(fwd_methods, fwd_method_labels);


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
        warning('No method results found — skipping source heatmaps.');
    else

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'source');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    n_shifts = numel(valid_labels);   % 18 (filtered)

    % Short labels: extract sign+magnitude+axis from valid_labels e.g. 'X+2mm' -> '+2X'
    src_shift_short = cell(1, n_shifts);
    for i = 1:n_shifts
        tok = regexp(valid_labels{i}, '([XYZ])([+-]\d+)mm', 'tokens');
        if ~isempty(tok)
            src_shift_short{i} = [tok{1}{2} tok{1}{1}];   % e.g. '+2X'
        else
            src_shift_short{i} = valid_labels{i};
        end
    end

    % Separator positions for within-method: after orig (1), after X group (1+6),
    % after Y group (1+12). Computed from valid_shift_axis.
    n_x = sum(valid_shift_axis == 1);
    n_y = sum(valid_shift_axis == 2);
    src_within_seps = [1, 1+n_x, 1+n_x+n_y];   % insert line AFTER these item indices

    % ----------------------------------------------------------------
    % SOURCE WITHIN-METHOD
    % ----------------------------------------------------------------
    fprintf('  Source within-method heatmaps...\n');
    for m_idx = 1:n_loaded
        method = loaded_methods{m_idx};
        mlabel = method_label_map(method);

        ref_key    = [method '_' sensitivity_ref_key];
        shift_keys = cellfun(@(k) [method '_' k], valid_keys_geom, ...
            'UniformOutput', false);

        item_keys   = [ref_key, shift_keys];
        item_labels = [{'orig'}, src_shift_short];
        n_items     = numel(item_keys);

        % Validate all keys exist
        item_ok = cellfun(@(k) isfield(leadfields, k), item_keys);
        if ~item_ok(1)
            warning('[%s] Reference not found — skipping within-method heatmap.', method);
            continue
        end

        for sens_ax = 1:n_axes
            min_s = get_min_sensors_lf(leadfields, item_keys(item_ok), ...
                orientation_labels{1}, sens_ax);

            [re_cell, cc_cell] = compute_pairwise_per_ori( ...
                leadfields, item_keys, item_ok, orientation_labels, ...
                sens_ax, src_range, min_s);

            n_eff   = n_items;
            fig_w   = max(1200, n_eff * 32 * n_ori + 300);
            fig_h   = max(700,  n_eff * 32 * 2 + 200);
            fig = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);
            tl  = tiledlayout(2, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('[%s]  Source shifts — within-method pairwise  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');

            for ori_idx = 1:n_ori
                ori   = orientation_labels{ori_idx};
                odsp  = orientation_display{ori_idx};

                % RE (top row)
                hax = nexttile(tl, ori_idx);
                draw_heatmap(hax, re_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  RE (%)'], cool, ...
                    [0, max_nonnan(re_cell.(ori)) * 100], 'RE (%)', ...
                    src_within_seps, '%.1f');

                % r² (bottom row)
                hax = nexttile(tl, n_ori + ori_idx);
                draw_heatmap(hax, cc_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  r² (%)'], flipud(cool), ...
                    [min_nonnan(cc_cell.(ori)) * 100, 100], 'r² (%)', ...
                    src_within_seps, '%.1f');
            end

            fname = sprintf('source_within_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
        fprintf('  [%s] Within-method done\n', method);
    end

    % ----------------------------------------------------------------
    % SOURCE CROSS-METHOD (one figure per shift axis)
    % ----------------------------------------------------------------
    fprintf('  Source cross-method heatmaps...\n');
    for sh_ax = 1:3
        ax_label = source_shift_axes{sh_ax};   % 'X', 'Y', or 'Z'
        ax_mask  = valid_shift_axis == sh_ax;
        ax_rows  = find(ax_mask);              % 6 row indices into valid_keys_geom

        for sens_ax = 1:n_axes
            item_keys   = {};
            item_labels = {};
            method_sep  = [];   % separator after these item indices

            pos = 0;
            for m_idx = 1:n_loaded
                method = loaded_methods{m_idx};
                abbrev = get_abbrev(method, abbrev_map);

                ref_key    = [method '_' sensitivity_ref_key];
                shift_keys = cellfun(@(k) [method '_' k], ...
                    valid_keys_geom(ax_rows), 'UniformOutput', false);

                item_keys   = [item_keys,   ref_key,          shift_keys]; %#ok<AGROW>
                item_labels = [item_labels, {[abbrev ':orig']}, ...
                    cellfun(@(l) [abbrev ':' l], src_shift_short(ax_rows), ...
                    'UniformOutput', false)];                              %#ok<AGROW>

                pos = pos + 1 + numel(ax_rows);
                if m_idx < n_loaded
                    method_sep(end+1) = pos;   %#ok<AGROW>
                end
            end

            n_items = numel(item_keys);
            item_ok = cellfun(@(k) isfield(leadfields, k), item_keys);

            if sum(item_ok) < 2
                fprintf('    Skipping cross-method %s axis, sensor axis %d — insufficient models.\n', ...
                    ax_label, sens_ax);
                continue
            end

            min_s = get_min_sensors_lf(leadfields, item_keys(item_ok), ...
                orientation_labels{1}, sens_ax);

            [re_cell, cc_cell] = compute_pairwise_per_ori( ...
                leadfields, item_keys, item_ok, orientation_labels, ...
                sens_ax, src_range, min_s);

            fig_w = max(1200, n_items * 32 * n_ori + 300);
            fig_h = max(700,  n_items * 32 * 2 + 200);
            fig = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);
            tl  = tiledlayout(2, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('Source shifts — cross-method pairwise  |  %s-axis shifts  |  Sensor axis %d of %d', ...
                ax_label, sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');

            for ori_idx = 1:n_ori
                ori  = orientation_labels{ori_idx};
                odsp = orientation_display{ori_idx};

                hax = nexttile(tl, ori_idx);
                draw_heatmap(hax, re_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  RE (%)'], cool, ...
                    [0, max_nonnan(re_cell.(ori)) * 100], 'RE (%)', ...
                    method_sep, '%.1f');

                hax = nexttile(tl, n_ori + ori_idx);
                draw_heatmap(hax, cc_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  r² (%)'], flipud(cool), ...
                    [min_nonnan(cc_cell.(ori)) * 100, 100], 'r² (%)', ...
                    method_sep, '%.1f');
            end

            fname = sprintf('source_cross_%s_sensorax%d', ax_label, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
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
        warning('No method results found — skipping sensor heatmaps.');
    else

    save_dir = fullfile(save_base_dir, 'perturbation_analysis', 'sensor');
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    n_sen_shifts = numel(valid_labels);   % 24 (filtered)

    % Short labels for sensor shifts: "B<b>-S<s>"
    sen_shift_short = cell(1, n_sen_shifts);
    for i = 1:n_sen_shifts
        sen_shift_short{i} = sprintf('B%d-S%d', valid_bundle_idx(i), valid_shift_idx(i));
    end

    % Separator positions for within-method: after orig (1), after bundle 1, after bundle 2
    b1_count = sum(valid_bundle_idx == 1);
    b2_count = sum(valid_bundle_idx == 2);
    sen_within_seps = [1, 1+b1_count, 1+b1_count+b2_count];

    % ----------------------------------------------------------------
    % SENSOR WITHIN-METHOD
    % ----------------------------------------------------------------
    fprintf('  Sensor within-method heatmaps...\n');
    for m_idx = 1:n_loaded
        method = loaded_methods{m_idx};
        mlabel = method_label_map(method);

        ref_key    = [method '_' sensor_sensitivity_ref_key];
        shift_keys = cellfun(@(k) [method '_' k], valid_keys_geom, ...
            'UniformOutput', false);

        item_keys   = [ref_key, shift_keys];
        item_labels = [{'orig'}, sen_shift_short];
        n_items     = numel(item_keys);

        item_ok = cellfun(@(k) isfield(leadfields, k), item_keys);
        if ~item_ok(1)
            warning('[%s] Sensor reference not found — skipping within-method heatmap.', method);
            continue
        end

        for sens_ax = 1:n_axes
            min_s = get_min_sensors_lf(leadfields, item_keys(item_ok), ...
                orientation_labels{1}, sens_ax);

            [re_cell, cc_cell] = compute_pairwise_per_ori( ...
                leadfields, item_keys, item_ok, orientation_labels, ...
                sens_ax, src_range, min_s);

            fig_w = max(1200, n_items * 30 * n_ori + 300);
            fig_h = max(700,  n_items * 30 * 2 + 200);
            fig = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);
            tl  = tiledlayout(2, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('[%s]  Sensor shifts — within-method pairwise  |  Sensor axis %d of %d', ...
                mlabel, sens_ax, n_axes), 'FontSize', 13, 'FontWeight', 'bold');

            for ori_idx = 1:n_ori
                ori  = orientation_labels{ori_idx};
                odsp = orientation_display{ori_idx};

                hax = nexttile(tl, ori_idx);
                draw_heatmap(hax, re_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  RE (%)'], cool, ...
                    [0, max_nonnan(re_cell.(ori)) * 100], 'RE (%)', ...
                    sen_within_seps, '%.1f');

                hax = nexttile(tl, n_ori + ori_idx);
                draw_heatmap(hax, cc_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  r² (%)'], flipud(cool), ...
                    [min_nonnan(cc_cell.(ori)) * 100, 100], 'r² (%)', ...
                    sen_within_seps, '%.1f');
            end

            fname = sprintf('sensor_within_%s_sensorax%d', method, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
        fprintf('  [%s] Within-method done\n', method);
    end

    % ----------------------------------------------------------------
    % SENSOR CROSS-METHOD (one figure per bundle)
    % ----------------------------------------------------------------
    fprintf('  Sensor cross-method heatmaps...\n');
    for b = 1:n_sensor_bundles
        b_mask = valid_bundle_idx == b;
        b_rows = find(b_mask);   % shift indices for this bundle

        for sens_ax = 1:n_axes
            item_keys   = {};
            item_labels = {};
            method_sep  = [];

            pos = 0;
            for m_idx = 1:n_loaded
                method = loaded_methods{m_idx};
                abbrev = get_abbrev(method, abbrev_map);

                ref_key    = [method '_' sensor_sensitivity_ref_key];
                shift_keys = cellfun(@(k) [method '_' k], ...
                    valid_keys_geom(b_rows), 'UniformOutput', false);

                item_keys   = [item_keys,   ref_key,          shift_keys]; %#ok<AGROW>
                item_labels = [item_labels, {[abbrev ':orig']}, ...
                    cellfun(@(l) [abbrev ':' l], sen_shift_short(b_rows), ...
                    'UniformOutput', false)];                              %#ok<AGROW>

                pos = pos + 1 + numel(b_rows);
                if m_idx < n_loaded
                    method_sep(end+1) = pos;   %#ok<AGROW>
                end
            end

            n_items = numel(item_keys);
            item_ok = cellfun(@(k) isfield(leadfields, k), item_keys);

            if sum(item_ok) < 2
                fprintf('    Skipping cross-method bundle %d, sensor axis %d — insufficient models.\n', ...
                    b, sens_ax);
                continue
            end

            min_s = get_min_sensors_lf(leadfields, item_keys(item_ok), ...
                orientation_labels{1}, sens_ax);

            [re_cell, cc_cell] = compute_pairwise_per_ori( ...
                leadfields, item_keys, item_ok, orientation_labels, ...
                sens_ax, src_range, min_s);

            fig_w = max(1200, n_items * 30 * n_ori + 300);
            fig_h = max(700,  n_items * 30 * 2 + 200);
            fig = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);
            tl  = tiledlayout(2, n_ori, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf('Sensor shifts — cross-method pairwise  |  %s  |  Sensor axis %d of %d', ...
                sensor_bundle_display{b}, sens_ax, n_axes), ...
                'FontSize', 13, 'FontWeight', 'bold');

            for ori_idx = 1:n_ori
                ori  = orientation_labels{ori_idx};
                odsp = orientation_display{ori_idx};

                hax = nexttile(tl, ori_idx);
                draw_heatmap(hax, re_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  RE (%)'], cool, ...
                    [0, max_nonnan(re_cell.(ori)) * 100], 'RE (%)', ...
                    method_sep, '%.1f');

                hax = nexttile(tl, n_ori + ori_idx);
                draw_heatmap(hax, cc_cell.(ori) * 100, item_labels, ...
                    [odsp '  —  r² (%)'], flipud(cool), ...
                    [min_nonnan(cc_cell.(ori)) * 100, 100], 'r² (%)', ...
                    method_sep, '%.1f');
            end

            fname = sprintf('sensor_cross_bundle%d_sensorax%d', b, sens_ax);
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('    Saved: %s\n', fname);
        end
    end

    end
    fprintf('Sensor heatmaps complete.\n\n');
end

fprintf('pt_plot_heatmaps complete.\n');


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [re_cell, cc_cell] = compute_pairwise_per_ori( ...
        leadfields, item_keys, item_ok, orientation_labels, ...
        sens_ax, src_range, min_s)
% Compute pairwise RE and r² for each orientation separately.
% Missing items (item_ok=false) produce NaN rows/cols.

    n       = numel(item_keys);
    n_ori   = numel(orientation_labels);
    re_cell = struct();
    cc_cell = struct();

    for oi = 1:n_ori
        ori = orientation_labels{oi};

        % Build L: one matrix per valid item, [min_s x n_sources]
        L = cell(1, n);
        for ki = 1:n
            if item_ok(ki)
                L{ki} = build_L_matrix(leadfields, item_keys{ki}, ori, ...
                    sens_ax, src_range, min_s);
            else
                L{ki} = [];   % placeholder — handled below
            end
        end

        % Compute pairwise metrics; skip missing items
        re_mat = nan(n, n);
        cc_mat = nan(n, n);
        for ii = 1:n
            if isempty(L{ii}); continue; end
            for jj = 1:n
                if isempty(L{jj}); continue; end
                [re_ij, cc_ij] = pairwise_re_cc(L{ii}, L{jj});
                re_mat(ii, jj) = re_ij;
                cc_mat(ii, jj) = cc_ij;
            end
        end

        re_cell.(ori) = re_mat;
        cc_cell.(ori) = cc_mat;
    end
end


function M = build_L_matrix(leadfields, key, ori_label, sens_ax, src_range, min_s)
% Build [min_s x n_src_plot] leadfield matrix for one model/orientation/axis.
    n_src = numel(src_range);
    M     = zeros(min_s, n_src);
    for si = 1:n_src
        src    = src_range(si);
        v      = leadfields.(key).(ori_label){sens_ax, src};
        M(:, si) = v(1:min_s);
    end
end


function [re_val, cc_val] = pairwise_re_cc(La, Lb)
% Median RE and r² between two [n_sensors x n_sources] matrices.
    n_src = size(La, 2);
    e     = zeros(1, n_src);
    c     = zeros(1, n_src);
    for s = 1:n_src
        va   = La(:, s);
        vb   = Lb(:, s);
        e(s) = norm(vb - va, 1) / (norm(va, 1) + norm(vb, 1));
        tmp  = corrcoef(va, vb);
        c(s) = tmp(1, 2)^2;
    end
    re_val = median(e, 'omitnan');
    cc_val = median(c, 'omitnan');
end


function min_s = get_min_sensors_lf(leadfields, keys, ori_label, sens_ax)
% Minimum sensor count across a set of valid leadfield keys.
    min_s = inf;
    for ki = 1:numel(keys)
        v     = leadfields.(keys{ki}).(ori_label){sens_ax, 1};
        min_s = min(min_s, numel(v));
    end
end


function draw_heatmap(hax, data, labels, ttl, cmap, clo, cbar_lbl, sep_lines, fmt)
% Draw annotated pairwise heatmap on axes hax.
    n = size(data, 1);

    imagesc(hax, data);
    colormap(hax, cmap);
    if clo(1) < clo(2)
        clim(hax, clo);
    end
    cb = colorbar(hax);
    cb.Label.String   = cbar_lbl;
    cb.Label.FontSize = 10;

    xticks(hax, 1:n); xticklabels(hax, labels); xtickangle(hax, 45);
    yticks(hax, 1:n); yticklabels(hax, labels);
    set(hax, 'FontSize', 8, 'TickDir', 'out');
    title(hax, ttl, 'FontSize', 10, 'FontWeight', 'bold');

    % Separator lines between groups
    for s = sep_lines
        if s >= 1 && s < n
            xline(hax, s + 0.5, 'k-', 'LineWidth', 1.8);
            yline(hax, s + 0.5, 'k-', 'LineWidth', 1.8);
        end
    end

    % Cell annotations
    crange = clo(2) - clo(1);
    for r = 1:n
        for c = 1:n
            val = data(r, c);
            if isnan(val); continue; end
            if crange > 0 && (val - clo(1)) / crange < 0.45
                tcol = 'w';
            else
                tcol = 'k';
            end
            text(hax, c, r, sprintf(fmt, val), ...
                'HorizontalAlignment', 'center', 'FontSize', 7, ...
                'FontWeight', 'bold', 'Color', tcol);
        end
    end
end


function v = max_nonnan(M)
    v = max(M(~isnan(M(:))));
    if isempty(v); v = 1; end
end

function v = min_nonnan(M)
    v = min(M(~isnan(M(:))));
    if isempty(v); v = 0; end
end

function abbrev = get_abbrev(method, abbrev_map)
    if isKey(abbrev_map, method)
        abbrev = abbrev_map(method);
    else
        abbrev = upper(method(1:min(3, numel(method))));
    end
end

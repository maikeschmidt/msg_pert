% sim_plot_noise_curves - r-squared vs noise level for SQUID MSG, OP-MSG, ESG
%
% Step 4 (plotting) of the realistic-measurement analysis. Reads the r-squared
% values produced by sim_simulate_noise and plots how each system degrades as
% its sensor noise floor is swept up and down around its published baseline.
%
% USAGE:
%   sim_plot_noise_curves
%
% OUTPUT (to <sim_save_dir>/noise_curves/):
%   noise_curve_average_<array>.png / .fig
%       r^2 vs noise level, averaged across every source on the cord.
%       One panel per dipole orientation. Shaded band = inter-quartile range
%       across sources, so the spread along the cord is visible.
%
%   noise_curve_src<MM>mm_<array>.png / .fig
%       Same, for the single source at sim_focus_src_mm.
%
%   noise_curve_vs_position_<array>.png / .fig
%       r^2 vs cord position, one line per noise level, one panel per system.
%       Shows WHERE along the cord each system holds up.
%
%   noise_curve_table.tsv
%       Mean r^2 (across cord) per system, orientation, and noise level.
%
% X-AXIS CONVENTION:
%   The x-axis is the noise level as a MULTIPLE of each system's own baseline,
%   not an absolute noise value. MSG is measured in fT and ESG in uV, so their
%   absolute noise floors cannot share an axis. "x baseline" can: it asks how
%   far each system is from its own real-world operating point. The absolute
%   baseline for each system is spelled out in the legend.
%
% DEPENDENCIES:
%   config_sim, sim_noise_rsq.mat (from sim_simulate_noise)
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_pert
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

clearvars
close all
clc

config_sim;

rsq_file = fullfile(sim_out_dir, 'sim_noise_rsq.mat');
if ~isfile(rsq_file)
    error('Not found: %s\nRun sim_simulate_noise first.', rsq_file);
end
load(rsq_file);   %#ok<LOAD>

save_dir = fullfile(sim_save_dir, 'noise_curves');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

n_sys = numel(sys_labels);
n_ori = numel(sim_orientations);
n_lev = numel(sim_noise_factors);

% Legend entries carry the absolute baseline so the "x baseline" axis stays
% interpretable despite the unit mismatch between MSG and ESG.
leg_labels = cell(1, n_sys);
for k = 1:n_sys
    leg_labels{k} = sprintf('%s  (1x = %g %s, %.0f Hz)', ...
        sys_labels{k}, sys_baseline(k), sys_units{k}, bw_eff(k));
end

fprintf('sim_plot_noise_curves\n');
fprintf('  %d systems, %d orientations, %d noise levels\n\n', n_sys, n_ori, n_lev);


% =========================================================================
% FIGURE 1: r^2 vs noise level, averaged across the cord
% =========================================================================

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 5*n_ori, 4.5]);
tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');

for o = 1:n_ori
    ax = nexttile(tl, o);
    hold(ax, 'on');

    h = gobjects(1, n_sys);
    for k = 1:n_sys
        % [n_src x n_lev] for this system + orientation
        R = squeeze(rsq_mean(k, o, :, :));

        med = median(R, 1, 'omitnan');
        q1  = prctile(R, 25, 1);
        q3  = prctile(R, 75, 1);

        col = sys_colors(k, :);

        % IQR band across sources
        fill(ax, [sim_noise_factors, fliplr(sim_noise_factors)], ...
                 [q1, fliplr(q3)], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        h(k) = plot(ax, sim_noise_factors, med, '-', ...
            'Color', col, 'LineWidth', pub_line_width, ...
            'Marker', 'o', 'MarkerSize', pub_marker_size, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', col);
    end

    % Baseline marker — the real-world operating point of every system
    xline(ax, 1, '--', 'baseline', 'Color', [0.4 0.4 0.4], ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    set(ax, 'XScale', 'log');
    set(ax, 'XTick', sim_noise_factors);
    set(ax, 'XTickLabel', arrayfun(@(f) sprintf('%gx', f), ...
        sim_noise_factors, 'UniformOutput', false));
    xlabel(ax, 'Sensor noise (\times system baseline)');
    ylabel(ax, 'r^2  vs noise-free field');
    ylim(ax, [0, 1.02]);
    grid(ax, 'on');
    box(ax, 'on');
    title(ax, sim_ori_display{o}, 'FontWeight', 'bold');

    if o == 1
        legend(ax, h, leg_labels, 'Location', 'southwest', ...
            'Interpreter', 'tex', 'FontSize', 8);
    end
end

title(tl, sprintf(['Sensor noise vs recoverable field  |  %s array  |  ' ...
                   'median across cord (band = IQR)'], sim_array), ...
    'FontSize', 13, 'FontWeight', 'bold');
subtitle(tl, sprintf(['evoked burst %g nA\\cdotm @ %g Hz, averaged over %d trials; ' ...
                      'white noise over each system''s own bandwidth (%s Hz)'], ...
    sim_dipole_nAm, sim_freq, sim_n_trials, ...
    strjoin(arrayfun(@(b) sprintf('%.0f', b), bw_eff, 'UniformOutput', false), '/')), ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);

fname = sprintf('noise_curve_average_%s', sim_array);
exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
saveas(fig, fullfile(save_dir, [fname '.fig']));
close(fig);
fprintf('  Saved: %s\n', fname);


% =========================================================================
% FIGURE 2: r^2 vs noise level at one chosen source
% =========================================================================

[~, focus_idx] = min(abs(src_mm - sim_focus_src_mm));
focus_mm       = src_mm(focus_idx);

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 5*n_ori, 4.5]);
tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');

for o = 1:n_ori
    ax = nexttile(tl, o);
    hold(ax, 'on');

    h = gobjects(1, n_sys);
    for k = 1:n_sys
        vals = squeeze(rsq_mean(k, o, focus_idx, :))';
        sds  = squeeze(rsq_sd(k, o, focus_idx, :))';
        col  = sys_colors(k, :);

        % Spread across noise realisations at this one source
        fill(ax, [sim_noise_factors, fliplr(sim_noise_factors)], ...
                 [vals - sds, fliplr(vals + sds)], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        h(k) = plot(ax, sim_noise_factors, vals, '-', ...
            'Color', col, 'LineWidth', pub_line_width, ...
            'Marker', 'o', 'MarkerSize', pub_marker_size, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', col);
    end

    xline(ax, 1, '--', 'baseline', 'Color', [0.4 0.4 0.4], ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    set(ax, 'XScale', 'log');
    set(ax, 'XTick', sim_noise_factors);
    set(ax, 'XTickLabel', arrayfun(@(f) sprintf('%gx', f), ...
        sim_noise_factors, 'UniformOutput', false));
    xlabel(ax, 'Sensor noise (\times system baseline)');
    ylabel(ax, 'r^2  vs noise-free field');
    ylim(ax, [0, 1.02]);
    grid(ax, 'on');
    box(ax, 'on');
    title(ax, sim_ori_display{o}, 'FontWeight', 'bold');

    if o == 1
        legend(ax, h, leg_labels, 'Location', 'southwest', ...
            'Interpreter', 'tex', 'FontSize', 8);
    end
end

title(tl, sprintf(['Sensor noise vs recoverable field  |  %s array  |  ' ...
                   'source at %.0f mm'], sim_array, focus_mm), ...
    'FontSize', 13, 'FontWeight', 'bold');
subtitle(tl, sprintf('band = \\pm1 s.d. across %d noise realisations', ...
    sim_n_realisations), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

fname = sprintf('noise_curve_src%03dmm_%s', round(focus_mm), sim_array);
exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
saveas(fig, fullfile(save_dir, [fname '.fig']));
close(fig);
fprintf('  Saved: %s\n', fname);


% =========================================================================
% FIGURE 3: r^2 vs cord position, one line per noise level
% =========================================================================
% Rows = orientation, columns = system. Shows where along the cord each system
% still recovers the field once noise is added.

lev_cols = parula(n_lev);

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [1, 1, 4.5*n_sys, 3.4*n_ori]);
tl = tiledlayout(n_ori, n_sys, 'TileSpacing', 'compact', 'Padding', 'tight');

for o = 1:n_ori
    for k = 1:n_sys
        ax = nexttile(tl, (o-1)*n_sys + k);
        hold(ax, 'on');

        for L = 1:n_lev
            plot(ax, src_mm, squeeze(rsq_mean(k, o, :, L)), '-', ...
                'Color', lev_cols(L, :), 'LineWidth', pub_line_width - 0.5);
        end

        xlabel(ax, 'Distance along cord (mm)');
        ylabel(ax, 'r^2');
        ylim(ax, [0, 1.02]);
        xlim(ax, [src_mm(1), src_mm(end)]);
        grid(ax, 'on');
        box(ax, 'on');
        title(ax, sprintf('%s — %s', sys_labels{k}, sim_ori_display{o}), ...
            'FontWeight', 'bold', 'FontSize', 10);
    end
end

cb = colorbar(ax);
colormap(ax, lev_cols);
caxis(ax, [0.5, n_lev + 0.5]);
cb.Ticks      = 1:n_lev;
cb.TickLabels = arrayfun(@(f) sprintf('%gx', f), ...
    sim_noise_factors, 'UniformOutput', false);
cb.Label.String = 'Noise (\times baseline)';
cb.Layout.Tile  = 'east';

title(tl, sprintf('r^2 along the cord at each noise level  |  %s array', sim_array), ...
    'FontSize', 13, 'FontWeight', 'bold');

fname = sprintf('noise_curve_vs_position_%s', sim_array);
exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
saveas(fig, fullfile(save_dir, [fname '.fig']));
close(fig);
fprintf('  Saved: %s\n', fname);


% =========================================================================
% TABLE: mean r^2 across the cord, per system / orientation / noise level
% =========================================================================

rows = {};
rows{end+1} = sprintf('%-12s\t%-16s\t%-10s\t%-14s\t%-10s\t%-10s', ...
    'System', 'Orientation', 'NoiseFactor', 'NoiseAbs', 'MeanRsq', 'MedianSNR');

for k = 1:n_sys
    for o = 1:n_ori
        for L = 1:n_lev
            R = squeeze(rsq_mean(k, o, :, L));
            S = squeeze(snr(k, o, :, L));
            rows{end+1} = sprintf('%-12s\t%-16s\t%-10g\t%-14.3g\t%-10.4f\t%-10.3g', ...
                sys_labels{k}, sim_ori_display{o}, sim_noise_factors(L), ...
                sigma_abs(k, L), mean(R, 'omitnan'), median(S, 'omitnan'));   %#ok<AGROW>
        end
    end
end

tbl_file = fullfile(save_dir, 'noise_curve_table.tsv');
fid = fopen(tbl_file, 'w');
fprintf(fid, '%s\n', rows{:});
fclose(fid);
fprintf('  Saved: noise_curve_table.tsv\n');

fprintf('\nsim_plot_noise_curves complete.\nFigures: %s\n', save_dir);

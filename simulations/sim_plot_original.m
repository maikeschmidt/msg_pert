% sim_plot_original - The per-system noise-curve figures for the ORIGINAL
%                     (unperturbed) geometry
%
% Reproduces the three "base" figures for the original geometry, overlaying all
% sensor systems (SQUID MSG, OP-MSG, ESG). The geometry-variant comparison
% (sim_plot_comparison) builds on top of this; this script is the base case.
%
% Reads the per-system results sim_run_geometries wrote for the baseline variant
% (<sim_out_dir>/<baseline>/sim_noise_<system>.mat), so run sim_run_geometries
% first.
%
% USAGE:
%   sim_plot_original
%
% OUTPUT (to <sim_save_dir>/original/):
%   original_curve_average_<array>.png/.fig
%       r^2 vs noise, median across the cord, band = inter-quartile range over
%       sources. One panel per orientation, all systems overlaid.
%   original_curve_src<MM>mm_<array>.png/.fig
%       Same, for the single source at sim_focus_src_mm; band = +/-1 s.d. across
%       noise realisations.
%   original_curve_vs_position_<array>.png/.fig
%       r^2 vs cord position, one line per noise level; rows = orientation,
%       columns = system.
%
% X-AXIS: noise as a multiple of each system's OWN baseline floor (MSG fT and
% ESG uV are not comparable in absolute terms). Absolute baselines are in the
% legend.
%
% DEPENDENCIES:
%   config_sim, per-system sim_noise_<system>.mat from sim_run_geometries
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

% Baseline (unperturbed) variant subfolder
gsel = find(strcmp({sim_geometries.group}, 'baseline'), 1);
if isempty(gsel); gsel = 1; end
base_name = sim_geometries(gsel).name;
in_sub    = fullfile(sim_out_dir, base_name);

save_dir = fullfile(sim_save_dir, 'original');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

nf    = sim_noise_factors;
n_ori = numel(sim_orientations);
n_lev = numel(nf);

% -------------------------------------------------------------------------
% Load each system's baseline-variant result
% -------------------------------------------------------------------------
S = struct('label', {}, 'short', {}, 'color', {}, 'leg', {}, ...
           'rsq_mean', {}, 'rsq_sd', {}, 'src_mm', {});
for k = 1:numel(sim_systems)
    f = fullfile(in_sub, sprintf('sim_noise_%s.mat', sim_systems(k).short));
    if ~isfile(f)
        fprintf('  %s: no baseline result (%s) — skipped\n', sim_systems(k).label, f);
        continue
    end
    d = load(f, 'rsq_mean', 'rsq_sd', 'src_mm', 'bw_eff');
    bw = d.bw_eff(min(k, numel(d.bw_eff)));
    S(end+1) = struct( ...
        'label', sim_systems(k).label, 'short', sim_systems(k).short, ...
        'color', sim_systems(k).color, ...
        'leg', sprintf('%s  (1x = %g %s, %.0f Hz)', sim_systems(k).label, ...
                       sim_systems(k).noise_baseline, sim_systems(k).noise_unit, bw), ...
        'rsq_mean', d.rsq_mean, 'rsq_sd', d.rsq_sd, 'src_mm', d.src_mm); %#ok<SAGROW>
end
if isempty(S)
    error('No baseline results found in %s. Run sim_run_geometries first.', in_sub);
end
n_sys  = numel(S);
src_mm = S(1).src_mm;

fprintf('sim_plot_original\n  %d systems on the "%s" geometry\n\n', n_sys, base_name);


% =========================================================================
% FIGURE 1: median across cord (band = IQR)
% =========================================================================
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 5*n_ori, 4.5]);
tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');
for o = 1:n_ori
    ax = nexttile(tl, o); hold(ax, 'on');
    h = gobjects(1, n_sys);
    for k = 1:n_sys
        R   = squeeze(S(k).rsq_mean(o, :, :));   % [n_src x n_lev]
        med = median(R, 1, 'omitnan');
        q1  = prctile(R, 25, 1);  q3 = prctile(R, 75, 1);
        col = S(k).color;
        fill(ax, [nf, fliplr(nf)], [q1, fliplr(q3)], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        h(k) = plot(ax, nf, med, '-', 'Color', col, 'LineWidth', pub_line_width, ...
            'Marker', 'o', 'MarkerSize', pub_marker_size, 'MarkerFaceColor', col);
    end
    xline(ax, 1, '--', 'baseline', 'Color', [0.4 0.4 0.4], ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    format_noise_axis(ax, nf);
    title(ax, sim_ori_display{o}, 'FontWeight', 'bold');
    if o == 1; legend(ax, h, {S.leg}, 'Location', 'southwest', 'FontSize', 8); end
end
title(tl, sprintf('Sensor noise vs recoverable field  |  %s array  |  median across cord (band = IQR)', ...
    sim_array), 'FontSize', 13, 'FontWeight', 'bold');
subtitle(tl, sprintf('evoked burst %g nA\\cdotm @ %g Hz, %d trials; original geometry', ...
    sim_dipole_nAm, sim_freq, sim_n_trials), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
save_fig(fig, save_dir, sprintf('original_curve_average_%s', sim_array));


% =========================================================================
% FIGURE 2: single source (band = +/-1 s.d. over realisations)
% =========================================================================
[~, fidx]  = min(abs(src_mm - sim_focus_src_mm));
focus_mm   = src_mm(fidx);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 5*n_ori, 4.5]);
tl  = tiledlayout(1, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');
for o = 1:n_ori
    ax = nexttile(tl, o); hold(ax, 'on');
    h = gobjects(1, n_sys);
    for k = 1:n_sys
        vals = squeeze(S(k).rsq_mean(o, fidx, :))';
        sds  = squeeze(S(k).rsq_sd(o, fidx, :))';
        col  = S(k).color;
        fill(ax, [nf, fliplr(nf)], [max(vals-sds,0), fliplr(min(vals+sds,1))], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        h(k) = plot(ax, nf, vals, '-', 'Color', col, 'LineWidth', pub_line_width, ...
            'Marker', 'o', 'MarkerSize', pub_marker_size, 'MarkerFaceColor', col);
    end
    xline(ax, 1, '--', 'baseline', 'Color', [0.4 0.4 0.4], ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    format_noise_axis(ax, nf);
    title(ax, sim_ori_display{o}, 'FontWeight', 'bold');
    if o == 1; legend(ax, h, {S.leg}, 'Location', 'southwest', 'FontSize', 8); end
end
title(tl, sprintf('Sensor noise vs recoverable field  |  %s array  |  source at %.0f mm', ...
    sim_array, focus_mm), 'FontSize', 13, 'FontWeight', 'bold');
subtitle(tl, sprintf('band = \\pm1 s.d. across %d noise realisations; original geometry', ...
    sim_n_realisations), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
save_fig(fig, save_dir, sprintf('original_curve_src%03dmm_%s', round(focus_mm), sim_array));


% =========================================================================
% FIGURE 3: r^2 along cord, one line per noise level
% =========================================================================
lev_cols = parula(n_lev);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 4.5*n_sys, 3.4*n_ori]);
tl = tiledlayout(n_ori, n_sys, 'TileSpacing', 'compact', 'Padding', 'tight');
for o = 1:n_ori
    for k = 1:n_sys
        ax = nexttile(tl, (o-1)*n_sys + k); hold(ax, 'on');
        for L = 1:n_lev
            plot(ax, src_mm, squeeze(S(k).rsq_mean(o, :, L)), '-', ...
                'Color', lev_cols(L, :), 'LineWidth', pub_line_width - 0.5);
        end
        xlabel(ax, 'Distance along cord (mm)'); ylabel(ax, 'r^2');
        ylim(ax, [0, 1.02]); xlim(ax, [src_mm(1), src_mm(end)]);
        grid(ax, 'on'); box(ax, 'on');
        title(ax, sprintf('%s — %s', S(k).label, sim_ori_display{o}), ...
            'FontWeight', 'bold', 'FontSize', 10);
    end
end
cb = colorbar(ax); colormap(ax, lev_cols); caxis(ax, [0.5, n_lev + 0.5]);
cb.Ticks = 1:n_lev;
cb.TickLabels = arrayfun(@(f) sprintf('%gx', f), nf, 'UniformOutput', false);
cb.Label.String = 'Noise (\times baseline)'; cb.Layout.Tile = 'east';
title(tl, sprintf('r^2 along the cord at each noise level  |  %s array  |  original geometry', ...
    sim_array), 'FontSize', 13, 'FontWeight', 'bold');
save_fig(fig, save_dir, sprintf('original_curve_vs_position_%s', sim_array));

fprintf('\nsim_plot_original complete.\nFigures: %s\n', save_dir);


% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------
function format_noise_axis(ax, nf)
    set(ax, 'XScale', 'log', 'XTick', nf, ...
        'XTickLabel', arrayfun(@(f) sprintf('%gx', f), nf, 'UniformOutput', false));
    xlabel(ax, 'Sensor noise (\times system baseline)');
    ylabel(ax, 'r^2  vs noise-free field');
    ylim(ax, [0, 1.02]); grid(ax, 'on'); box(ax, 'on');
end

function save_fig(fig, save_dir, name)
    exportgraphics(fig, fullfile(save_dir, [name '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [name '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', name);
end

% sim_plot_perturbation_noise - r^2 vs noise curves per perturbation, by bundle
%
% Plots the output of sim_perturbation_noise: for each perturbation family
% (source, sensor, conductivity) and each applicable sensor system, how the
% recoverable field degrades with sensor noise, grouped into small / medium /
% large perturbation bundles.
%
% USAGE:
%   sim_plot_perturbation_noise
%
% OUTPUT (to <sim_save_dir>/perturbation_noise/):
%   pertnoise_<pert>_<system>.png / .fig
%     One figure per perturbation type per system.
%     Rows    = the two r^2 references:
%                 top    "vs perfect"  — model error + noise together
%                 bottom "vs self"     — noise only, per perturbed config
%     Columns = dipole orientation.
%     Lines   = perturbation bundle (small / medium / large), each the mean
%               across that bundle's 8 shifts and across the whole cord.
%
%   pertnoise_summary_<system>.png / .fig
%     One compact figure per system: rows = perturbation type, columns =
%     orientation, showing only the "vs perfect" (combined) curves — the most
%     directly interpretable view for comparing perturbation families.
%
% X-AXIS:
%   Noise as a multiple of each system's own baseline floor (see config_sim for
%   why absolute noise cannot be shared across MSG and ESG).
%
% DEPENDENCIES:
%   config_sim, sim_pert_noise.mat (from sim_perturbation_noise)
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

nf_file = fullfile(sim_out_dir, 'sim_pert_noise.mat');
if ~isfile(nf_file)
    error('Not found: %s\nRun sim_perturbation_noise first.', nf_file);
end
load(nf_file);   %#ok<LOAD>

save_dir = fullfile(sim_save_dir, 'perturbation_noise');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

n_sys = numel(sys_labels);
n_ori = numel(sim_orientations);
nf    = sim_noise_factors;

pert_title = containers.Map( ...
    {'source', 'sensor', 'cond'}, ...
    {'Source-position perturbation', 'Sensor-position perturbation', ...
     'Tissue-conductivity perturbation'});

fprintf('sim_plot_perturbation_noise\n');
fprintf('  %d perturbation families, %d systems\n\n', numel(pert_names), n_sys);


% =========================================================================
% DETAILED FIGURE: one per (perturbation type, system)
% =========================================================================

for p = 1:numel(pert_names)
    pname = pert_names{p};
    R     = results.(pname);
    nb    = R.n_bundle;

    if pert_title.isKey(pname); ptitle = pert_title(pname); else; ptitle = pname; end

    for si = 1:n_sys

        fig = figure('Color', 'w', 'Units', 'inches', ...
            'Position', [1, 1, 4.6*n_ori, 8]);
        tl = tiledlayout(2, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');

        % Two rows: perfect reference (combined) then self reference (noise only)
        ref_data   = {R.cordmean_perf, R.cordmean_self};
        ref_titles = {'vs perfect field  (model error + noise)', ...
                      'vs same perturbed field  (noise only)'};

        for rr = 1:2
            D = ref_data{rr};   % [sys ori bundle lev]

            for o = 1:n_ori
                ax = nexttile(tl, (rr-1)*n_ori + o);
                hold(ax, 'on');

                h = gobjects(1, nb);
                for b = 1:nb
                    curve = squeeze(D(si, o, b, :))';
                    col   = R.bundle_colors(min(b, size(R.bundle_colors,1)), :);
                    h(b)  = plot(ax, nf, curve, '-', ...
                        'Color', col, 'LineWidth', pub_line_width, ...
                        'Marker', 'o', 'MarkerSize', pub_marker_size, ...
                        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
                end

                xline(ax, 1, '--', 'baseline', 'Color', [0.4 0.4 0.4], ...
                    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

                set(ax, 'XScale', 'log', 'XTick', nf, ...
                    'XTickLabel', arrayfun(@(f) sprintf('%gx', f), nf, ...
                                           'UniformOutput', false));
                xlabel(ax, 'Sensor noise (\times baseline)');
                ylabel(ax, 'r^2');
                ylim(ax, [0, 1.02]);
                grid(ax, 'on'); box(ax, 'on');

                if rr == 1
                    title(ax, sim_ori_display{o}, 'FontWeight', 'bold');
                end
                if o == 1
                    text(-0.28, 0.5, ref_titles{rr}, 'Units', 'normalized', ...
                        'Rotation', 90, 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', 'FontWeight', 'bold', ...
                        'FontSize', 9);
                    legend(ax, h, R.bundle_display, 'Location', 'southwest', ...
                        'FontSize', 8);
                end
            end
        end

        title(tl, sprintf('%s  |  %s', ptitle, sys_labels{si}), ...
            'FontSize', 13, 'FontWeight', 'bold');
        subtitle(tl, sprintf(['mean across bundle shifts and cord; evoked %g nA\\cdotm ' ...
                              '@ %g Hz, %d trials'], ...
            sim_dipole_nAm, sim_freq, sim_n_trials), ...
            'FontSize', 9, 'Color', [0.4 0.4 0.4]);

        fname = sprintf('pertnoise_%s_%s', pname, sys_shorts{si});
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end


% =========================================================================
% SUMMARY FIGURE: one per system, all perturbation families, combined ref
% =========================================================================

for si = 1:n_sys

    fig = figure('Color', 'w', 'Units', 'inches', ...
        'Position', [1, 1, 4.6*n_ori, 3.2*numel(pert_names)]);
    tl = tiledlayout(numel(pert_names), n_ori, ...
        'TileSpacing', 'compact', 'Padding', 'tight');

    for p = 1:numel(pert_names)
        pname = pert_names{p};
        R     = results.(pname);
        nb    = R.n_bundle;
        if pert_title.isKey(pname); ptitle = pert_title(pname); else; ptitle = pname; end

        for o = 1:n_ori
            ax = nexttile(tl, (p-1)*n_ori + o);
            hold(ax, 'on');

            h = gobjects(1, nb);
            for b = 1:nb
                curve = squeeze(R.cordmean_perf(si, o, b, :))';
                col   = R.bundle_colors(min(b, size(R.bundle_colors,1)), :);
                h(b)  = plot(ax, nf, curve, '-', ...
                    'Color', col, 'LineWidth', pub_line_width, ...
                    'Marker', 'o', 'MarkerSize', pub_marker_size - 1, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col);
            end

            xline(ax, 1, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
            set(ax, 'XScale', 'log', 'XTick', nf, ...
                'XTickLabel', arrayfun(@(f) sprintf('%gx', f), nf, ...
                                       'UniformOutput', false));
            xlabel(ax, 'Noise (\times baseline)');
            ylabel(ax, 'r^2');
            ylim(ax, [0, 1.02]);
            grid(ax, 'on'); box(ax, 'on');

            if p == 1
                title(ax, sim_ori_display{o}, 'FontWeight', 'bold');
            end
            if o == 1
                text(-0.30, 0.5, ptitle, 'Units', 'normalized', 'Rotation', 90, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontWeight', 'bold', 'FontSize', 9);
                legend(ax, h, R.bundle_display, 'Location', 'southwest', 'FontSize', 7);
            end
        end
    end

    title(tl, sprintf('Perturbation vs sensor noise (r^2 vs perfect field)  |  %s', ...
        sys_labels{si}), 'FontSize', 13, 'FontWeight', 'bold');

    fname = sprintf('pertnoise_summary_%s', sys_shorts{si});
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

fprintf('\nsim_plot_perturbation_noise complete.\nFigures: %s\n', save_dir);

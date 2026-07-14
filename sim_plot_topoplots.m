% sim_plot_topoplots - Perfect-forward-field topoplots across models
%
% Plots the noise-free ("perfect") forward field of a single source for every
% model in sim_models, so the spatial character of each can be compared
% directly.
%
% This is the figure that tests the smoothness hypothesis: if BEM MSG really
% is more individualised and spatially sharper than Biot-Savart MSG and BEM
% ESG, it should be visible here as tighter, more structured field patterns —
% and that sharpness is the candidate explanation for why BEM MSG turned out
% MORE sensitive to perturbation, not less.
%
% USAGE:
%   sim_plot_topoplots
%
% OUTPUT (to <sim_save_dir>/topoplots/):
%   perfect_topoplot_<array>_sensorax<N>_src<MM>mm.png / .fig
%   One figure per array (front/back) per sensor axis.
%   Rows = forward models, columns = dipole orientations.
%
% COLOUR LIMITS:
%   Shared per ROW (per model), not globally. MSG rows are fT/nAm and the ESG
%   row is uV/nAm — a shared colour axis across rows would be meaningless. A
%   per-row limit keeps the three orientation panels within a model directly
%   comparable, which is the comparison that matters.
%
% SENSOR AXES:
%   MSG (triaxial) has 3 sensor axes; ESG electrodes have 2. On the sensor-axis
%   3 figure the ESG row is therefore drawn as an explicit "no 3rd axis" panel
%   rather than being silently dropped.
%
% DEPENDENCIES:
%   config_sim, pt_add_functions
%   sim_load_leadfield(), sim_sensor_positions()   — msg_pert/functions/
%   plot_topoplot_publication()                    — msg_fwd/functions/
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
pt_add_functions;

save_dir = fullfile(sim_save_dir, 'topoplots');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

arrays    = {'front', 'back'};
n_models  = numel(sim_models);
n_ori     = numel(sim_orientations);

% Cord distance -> source index (sources are evenly spaced along the cord)
src_idx = round(sim_topo_src_mm / src_spacing_mm) + 1;

fprintf('sim_plot_topoplots\n');
fprintf('  Source: %.0f mm along cord (index %d)\n', sim_topo_src_mm, src_idx);
fprintf('  Models: %d   Arrays: front, back\n\n', n_models);


% =========================================================================
% LOAD ALL MODELS (both arrays)
% =========================================================================

lf_all  = cell(n_models, numel(arrays));   % {model, array}
pos_all = cell(n_models, numel(arrays));

for m = 1:n_models
    for a = 1:numel(arrays)
        arr  = arrays{a};
        file = sim_models(m).(arr);

        fprintf('  Loading %-20s %-5s ... ', sim_models(m).label, arr);
        lf_all{m, a} = sim_load_leadfield(file, sim_models(m).var, ...
            sim_models(m).scale, sim_models(m).is_meg);
        pos_all{m, a} = sim_sensor_positions(sim_models(m).geom_file, ...
            arr, sim_models(m).is_meg);
        fprintf('%d sources, %d axes, %d sensors/axis\n', ...
            lf_all{m, a}.n_sources, lf_all{m, a}.n_sensor_axes, ...
            lf_all{m, a}.n_sensors_per_axis);

        if src_idx > lf_all{m, a}.n_sources
            error(['Source index %d exceeds the %d sources in %s (%s).\n' ...
                   'Check sim_topo_src_mm (%.0f mm) and src_spacing_mm (%.1f mm).'], ...
                   src_idx, lf_all{m, a}.n_sources, sim_models(m).label, arr, ...
                   sim_topo_src_mm, src_spacing_mm);
        end
    end
end
fprintf('\n');


% =========================================================================
% ONE FIGURE PER ARRAY PER SENSOR AXIS
% =========================================================================
% Sensor axis count is taken as the maximum across models (3 for MSG). Models
% with fewer axes get an explicit placeholder panel on the extra figures.

max_axes = max(cellfun(@(l) l.n_sensor_axes, lf_all(:)));
axis_labels = {'X-axis', 'Y-axis', 'Z-axis'};

for a = 1:numel(arrays)
    arr = arrays{a};

    for ax = 1:max_axes

        fig = figure('Color', 'w', 'Units', 'inches', ...
            'Position', [1, 1, 3*n_ori + 1.5, 2.8*n_models]);
        tl = tiledlayout(n_models, n_ori, ...
            'TileSpacing', 'compact', 'Padding', 'tight');

        for m = 1:n_models
            lf  = lf_all{m, a};
            pos = pos_all{m, a};

            % ── Model has no such sensor axis (ESG on axis 3) ─────────────
            if ax > lf.n_sensor_axes
                for ori = 1:n_ori
                    nexttile(tl, (m-1)*n_ori + ori);
                    text(0.5, 0.5, sprintf('%s\nhas no sensor axis %d', ...
                        sim_models(m).label, ax), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 10, 'Color', [0.45 0.45 0.45], ...
                        'Interpreter', 'none');
                    axis off
                end
                continue
            end

            % ── Row-shared colour limit across the three orientations ─────
            row_max = 0;
            for ori = 1:n_ori
                vals    = lf.(sim_orientations{ori}){ax, src_idx};
                row_max = max(row_max, max(abs(vals)));
            end
            if row_max == 0; row_max = 1; end   % guard against an all-zero panel
            row_clim = [-row_max, row_max];

            for ori = 1:n_ori
                nexttile(tl, (m-1)*n_ori + ori);

                vals = lf.(sim_orientations{ori}){ax, src_idx};
                plot_topoplot_publication(pos{ax}, vals, row_clim, lf.is_meg);

                if m == 1
                    title(sim_ori_display{ori}, ...
                        'FontWeight', 'bold', 'FontSize', 11);
                end

                % Row label on the leftmost panel, with the unit spelled out
                if ori == 1
                    if lf.is_meg
                        unit_str = 'fT/nAm';
                    else
                        unit_str = '\muV/nAm';
                    end
                    ylabel_txt = sprintf('%s\n(%s)', sim_models(m).label, unit_str);
                    text(-0.12, 0.5, ylabel_txt, ...
                        'Units', 'normalized', ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 10);
                end
            end
        end

        title(tl, sprintf(['Perfect forward field  |  %s array  |  %s  |  ' ...
                           'source %.0f mm along cord'], ...
                upper(arr(1)), axis_labels{ax}, sim_topo_src_mm), ...
            'FontSize', 13, 'FontWeight', 'bold');
        subtitle(tl, sprintf(['%s array — colour limits shared within each row; ' ...
                              'MSG rows are fT/nAm, ESG row is \\muV/nAm'], arr), ...
            'FontSize', 9, 'Color', [0.4 0.4 0.4]);

        fname = sprintf('perfect_topoplot_%s_sensorax%d_src%03dmm', ...
            arr, ax, round(sim_topo_src_mm));
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('  Saved: %s\n', fname);
    end
end

fprintf('\nsim_plot_topoplots complete.\nFigures: %s\n', save_dir);

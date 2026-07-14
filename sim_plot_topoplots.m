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
            sim_models(m).scale, sim_models(m).is_meg, sim_models(m).n_axes);
        pos_all{m, a} = sim_sensor_positions(sim_models(m).geom_file, ...
            arr, sim_models(m).is_meg, lf_all{m, a}.n_sensor_axes);
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
% AMPLITUDE DIAGNOSTIC
% =========================================================================
% Peak |leadfield| at this source, per model. Cross-check these against the
% peak-amplitude figures from msg_fwd (plot_absmax_curves / plot_sm_absmax):
% the two MSG models should land in the same ballpark, since Biot-Savart and
% BEM agree closely for MSG. If they differ by a round factor here (1e3, 1e6,
% 1e9, 1e15) the unit scale in config_sim is wrong for that model, not the
% physics.
%
% Expected scales depend on whether the FieldTrip dipoleunit patch is active
% (see bem_patched in config_sim):
%   Patched   — BEM MSG x 1e15,  BEM ESG x 1e6
%   Unpatched — BEM MSG x 1e6,   BEM ESG x 1e-3   (leadfield is per A*m, not nA*m)
%   Biot-Savart x 1 either way (it writes fT/nAm directly, bypassing FieldTrip)

fprintf('  Peak |leadfield| at source %d (compare with msg_fwd absmax figures):\n', src_idx);
for m = 1:n_models
    if sim_models(m).is_meg
        unit_txt = 'fT/nAm';
    else
        unit_txt = 'uV/nAm';
    end
    fprintf('    %-20s scale x%-8g ', sim_models(m).label, sim_models(m).scale);
    for a = 1:numel(arrays)
        lf = lf_all{m, a};
        pk = 0;
        for ax = 1:lf.n_sensor_axes
            for ori = 1:n_ori
                pk = max(pk, max(abs(lf.(sim_orientations{ori}){ax, src_idx})));
            end
        end
        fprintf('%s: %9.3g %s   ', arrays{a}, pk, unit_txt);
    end
    fprintf('\n');
end

% The two MSG models should agree to within a factor of a few. A round-number
% ratio is the signature of a unit-scale error.
msg_idx = find([sim_models.is_meg]);
if numel(msg_idx) == 2
    pk = zeros(1, 2);
    for i = 1:2
        lf = lf_all{msg_idx(i), 2};   % back array
        for ax = 1:lf.n_sensor_axes
            for ori = 1:n_ori
                pk(i) = max(pk(i), max(abs(lf.(sim_orientations{ori}){ax, src_idx})));
            end
        end
    end
    ratio = pk(1) / pk(2);
    fprintf('    MSG %s / %s peak ratio = %.4g\n', ...
        sim_models(msg_idx(1)).label, sim_models(msg_idx(2)).label, ratio);
    if ratio > 20 || ratio < 0.05
        warning(['The two MSG models differ by %.4g x. msg_fwd shows Biot-Savart ' ...
                 'and BEM agreeing closely for MSG, so a ratio this large points ' ...
                 'to a unit-scale mismatch in config_sim (check .scale on each ' ...
                 'model against msg_fwd/load_and_organise_leadfields).'], ratio);
    end
end
fprintf('\n');


% =========================================================================
% ONE FIGURE PER ARRAY PER COMPARISON SLOT
% =========================================================================
% Figures are organised by SLOT, not by raw axis index, because the axis index
% is not a shared quantity across models. Slots follow the MSG convention
% (1 = X, 2 = Y, 3 = Z); each model declares which slot each of its own axes
% belongs in via .axis_slot.
%
% For ESG that mapping is:  tangential -> slot 1 (X),  radial -> slot 3 (Z).
% Note this is NOT the ESG channel order — radial is ESG's SECOND axis but
% belongs in the THIRD slot. Indexing by raw axis number would have put ESG
% radial next to MSG Y, comparing two unrelated measurements. ESG has no
% counterpart to MSG Y, so the slot-2 figure shows a placeholder for its row.

slot_names = {'X-axis', 'Y-axis', 'Z-axis'};
n_slots    = numel(slot_names);

% Every axis a model actually has must be named and assigned to a slot.
for m = 1:n_models
    for a = 1:numel(arrays)
        n_ax_m = lf_all{m, a}.n_sensor_axes;
        if numel(sim_models(m).axis_names) < n_ax_m || ...
           numel(sim_models(m).axis_slot)  < n_ax_m
            error(['Model "%s" has %d sensor axes in its %s leadfield but only ' ...
                   '%d axis_names and %d axis_slot entries.\nSet both in ' ...
                   'config_sim.m so every axis is named and placed.'], ...
                   sim_models(m).label, n_ax_m, arrays{a}, ...
                   numel(sim_models(m).axis_names), numel(sim_models(m).axis_slot), m);
        end
    end
end

for a = 1:numel(arrays)
    arr = arrays{a};

    for slot = 1:n_slots

        fig = figure('Color', 'w', 'Units', 'inches', ...
            'Position', [1, 1, 3*n_ori + 1.5, 2.8*n_models]);
        tl = tiledlayout(n_models, n_ori, ...
            'TileSpacing', 'compact', 'Padding', 'tight');

        for m = 1:n_models
            lf  = lf_all{m, a};
            pos = pos_all{m, a};

            % Which of THIS model's axes (if any) sits in this slot?
            ax = find(sim_models(m).axis_slot(1:lf.n_sensor_axes) == slot, 1);

            % ── Model has no axis in this slot (ESG in the Y slot) ────────
            if isempty(ax)
                for ori = 1:n_ori
                    nexttile(tl, (m-1)*n_ori + ori);
                    text(0.5, 0.5, sprintf('%s\nhas no %s equivalent\n(%s only)', ...
                        sim_models(m).label, slot_names{slot}, ...
                        strjoin(sim_models(m).axis_names(1:lf.n_sensor_axes), ' + ')), ...
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

                % Row label on the leftmost panel: model, ITS name for this
                % axis, and the unit — so the row is self-describing and is
                % never mistaken for the same measurement as the row above.
                if ori == 1
                    if lf.is_meg
                        unit_str = 'fT/nAm';
                    else
                        unit_str = '\muV/nAm';
                    end
                    ylabel_txt = sprintf('%s\n%s (%s)', sim_models(m).label, ...
                        sim_models(m).axis_names{ax}, unit_str);
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
                arr, slot_names{slot}, sim_topo_src_mm), ...
            'FontSize', 13, 'FontWeight', 'bold');
        subtitle(tl, ['MSG X/Y/Z are field components; ESG tangential aligns with X ' ...
                      'and radial with Z. Colour limits shared within a row.'], ...
            'FontSize', 9, 'Color', [0.4 0.4 0.4]);

        fname = sprintf('perfect_topoplot_%s_%s_src%03dmm', ...
            arr, lower(strrep(slot_names{slot}, '-', '')), round(sim_topo_src_mm));
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('  Saved: %s\n', fname);
    end
end

fprintf('\nsim_plot_topoplots complete.\nFigures: %s\n', save_dir);

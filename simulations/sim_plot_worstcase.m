% sim_plot_worstcase - Systems compared under the WORST-CASE perturbation
%
% One figure that puts the three sensor systems (SQUID MSG, OP-MSG, ESG) on the
% same axes for the largest perturbation of each type, scored against the
% ORIGINAL (unperturbed) field. This directly answers: for the biggest shift,
% how differently do the systems hold up?
%
% Layout: rows = perturbation family (source / sensor / conductivity),
%         columns = dipole orientation,
%         lines  = the three systems (each its large-shift variant).
%
% r^2 REFERENCE — the original, not self:
%   Unlike the comparison curves (which score each variant against its OWN
%   noise-free field), this figure scores the noisy WORST-CASE field against the
%   unperturbed original. So r^2 begins below 1 even at zero noise, by the amount
%   the large shift alone already changed the field — which is exactly "the
%   difference the biggest shift has". It is computed fresh here (perturbed field
%   vs original field, under noise) rather than read from the per-variant files,
%   which store the self-referenced r^2.
%
% USAGE:
%   sim_plot_worstcase
%
% OUTPUT (to <sim_save_dir>/comparison/):
%   worstcase_systems_<array>.png / .fig
%
% DEPENDENCIES:
%   config_sim, sim_lf_path(), sim_load_leadfield(), sim_evoked_noise_rsq()
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

worst_bundle = 'large';   % SET THIS: which bundle counts as "worst case"

save_dir = fullfile(sim_save_dir, 'comparison');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

nf    = sim_noise_factors;
n_ori = numel(sim_orientations);
n_lev = numel(nf);
n_sys = numel(sim_systems);

group_order = {'source', 'sensor', 'cond'};
group_title = containers.Map( ...
    {'source', 'sensor', 'cond'}, ...
    {'Source shift (large)', 'Sensor shift (large)', 'Conductivity change (large)'});

% Noise s.d. per system per level (bandwidth + trial averaging), as elsewhere
nyquist    = sim_fs / 2;
trial_gain = sqrt(sim_n_trials);
sigma_abs  = zeros(n_sys, n_lev);
for k = 1:n_sys
    bw = min(sim_systems(k).bandwidth_hz, nyquist);
    sigma_abs(k, :) = sim_systems(k).noise_baseline * nf * sqrt(bw) / trial_gain;
end
w = sim_waveform(:);

% Baseline variant (the original the worst case is scored against)
bsel = find(strcmp({sim_geometries.group}, 'baseline'), 1);
if isempty(bsel)
    error('No baseline variant in sim_geometries to use as the original reference.');
end
base_geo = sim_geometries(bsel);

fprintf('sim_plot_worstcase\n');
fprintf('  worst-case bundle = "%s", array = %s\n\n', worst_bundle, sim_array);


% -------------------------------------------------------------------------
% Precompute the ORIGINAL field vectors per system (the reference), once
% -------------------------------------------------------------------------
orig_V = cell(1, n_sys);   % orig_V{k}{o} = [n_ch x n_src]
for k = 1:n_sys
    model = sim_models(sim_systems(k).model);
    p = sim_lf_path(model, base_geo, sim_array);
    if ~isfile(p)
        fprintf('  %s: original leadfield missing — system skipped\n', sim_systems(k).label);
        continue
    end
    lf = sim_load_leadfield(p, model.var, model.scale, model.is_meg, model.n_axes);
    orig_V{k} = flatten_orientations(lf, sim_orientations);
end


% -------------------------------------------------------------------------
% Compute worst-case r^2 vs original, per group / system / orientation
% -------------------------------------------------------------------------
% curves.(group){k} = struct with cm/lo/hi [n_ori x n_lev], or [] if absent
curves = struct();
present_groups = {};

for gi = 1:numel(group_order)
    grp = group_order{gi};

    % Find the large-bundle variant for this group
    vsel = find(strcmp({sim_geometries.group}, grp) & ...
                strcmp({sim_geometries.bundle}, worst_bundle), 1);
    if isempty(vsel)
        fprintf('  [%s] no "%s" variant — row skipped\n', grp, worst_bundle);
        continue
    end
    geo = sim_geometries(vsel);
    present_groups{end+1} = grp; %#ok<SAGROW>

    per_sys = cell(1, n_sys);
    for k = 1:n_sys
        if isempty(orig_V{k}); continue; end
        model = sim_models(sim_systems(k).model);

        p = sim_lf_path(model, geo, sim_array);
        if ~isfile(p)
            fprintf('  [%s] %s: no leadfield — line skipped\n', grp, sim_systems(k).label);
            continue
        end
        if strcmp(geo.kind, 'cond'); sc = model.cond_scale; else; sc = model.scale; end
        lf = sim_load_leadfield(p, model.var, sc, model.is_meg, model.n_axes);

        V = flatten_orientations(lf, sim_orientations);
        if ~isequal(size(V{1}), size(orig_V{k}{1}))
            fprintf('  [%s] %s: size mismatch vs original — skipped\n', grp, sim_systems(k).label);
            continue
        end

        cm = nan(n_ori, n_lev); lo = nan(n_ori, n_lev); hi = nan(n_ori, n_lev);
        for o = 1:n_ori
            rm = sim_evoked_noise_rsq(V{o}, w, sigma_abs(k, :), ...
                sim_n_realisations, orig_V{k}{o});   % r^2 vs ORIGINAL
            cm(o, :) = mean(rm, 1, 'omitnan');
            lo(o, :) = prctile(rm, 25, 1);
            hi(o, :) = prctile(rm, 75, 1);
        end
        per_sys{k} = struct('cm', cm, 'lo', lo, 'hi', hi);
        fprintf('  [%s] %s done\n', grp, sim_systems(k).label);
    end
    curves.(grp) = per_sys;
end

if isempty(present_groups)
    error('No worst-case variants found. Check sim_geometries has "%s" bundles.', worst_bundle);
end


% -------------------------------------------------------------------------
% Plot: rows = groups, cols = orientations, lines = systems
% -------------------------------------------------------------------------
n_grp = numel(present_groups);
fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [1, 1, 4.6*n_ori, 3.2*n_grp]);
tl = tiledlayout(n_grp, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');

for gi = 1:n_grp
    grp     = present_groups{gi};
    per_sys = curves.(grp);
    if group_title.isKey(grp); gtitle = group_title(grp); else; gtitle = grp; end

    for o = 1:n_ori
        ax = nexttile(tl, (gi-1)*n_ori + o); hold(ax, 'on');
        leg_h = gobjects(0); leg_t = {};

        for k = 1:n_sys
            if isempty(per_sys{k}); continue; end
            col = sim_systems(k).color;
            fill(ax, [nf, fliplr(nf)], ...
                 [per_sys{k}.lo(o, :), fliplr(per_sys{k}.hi(o, :))], col, ...
                'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            hh = plot(ax, nf, per_sys{k}.cm(o, :), '-', 'Color', col, ...
                'LineWidth', pub_line_width, 'Marker', 'o', ...
                'MarkerSize', pub_marker_size, 'MarkerFaceColor', col);
            leg_h(end+1) = hh; leg_t{end+1} = sim_systems(k).label; %#ok<SAGROW>
        end

        xline(ax, 1, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        set(ax, 'XScale', 'log', 'XTick', nf, ...
            'XTickLabel', arrayfun(@(f) sprintf('%gx', f), nf, 'UniformOutput', false));
        xlabel(ax, 'Sensor noise (\times system baseline)');
        ylabel(ax, 'r^2  vs original field');
        ylim(ax, [0, 1.02]); grid(ax, 'on'); box(ax, 'on');

        if gi == 1; title(ax, sim_ori_display{o}, 'FontWeight', 'bold'); end
        if o == 1
            text(-0.30, 0.5, gtitle, 'Units', 'normalized', 'Rotation', 90, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', 'FontSize', 10);
            if ~isempty(leg_h)
                legend(ax, leg_h, leg_t, 'Location', 'southwest', 'FontSize', 8);
            end
        end
    end
end

title(tl, sprintf('Worst-case perturbation vs sensor noise, by system  |  %s array', sim_array), ...
    'FontSize', 13, 'FontWeight', 'bold');
subtitle(tl, sprintf(['largest shift per family; evoked %g nA\\cdotm @ %g Hz, %d trials; ' ...
                      'line = mean across cord, band = IQR;  r^2 vs ORIGINAL field'], ...
    sim_dipole_nAm, sim_freq, sim_n_trials), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

fname = sprintf('worstcase_systems_%s', sim_array);
exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
saveas(fig, fullfile(save_dir, [fname '.fig']));
close(fig);
fprintf('\nsim_plot_worstcase complete.\nSaved: %s\n', fullfile(save_dir, fname));


% -------------------------------------------------------------------------
% Local function: flatten a leadfield to {o} = [n_ch x n_src] per orientation
% -------------------------------------------------------------------------
function V = flatten_orientations(lf, orientations)
    n_ori = numel(orientations);
    n_src = lf.n_sources;
    n_ch  = lf.n_sensor_axes * lf.n_sensors_per_axis;
    V = cell(1, n_ori);
    for o = 1:n_ori
        M = zeros(n_ch, n_src);
        for s = 1:n_src
            M(:, s) = vertcat(lf.(orientations{o}){:, s});
        end
        V{o} = M;
    end
end

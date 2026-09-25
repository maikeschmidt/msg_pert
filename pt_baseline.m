% pt_baseline - Stage 1: the unperturbed, noise-free forward fields
%
% Establishes what the three forward models predict before anything is
% perturbed or any noise is added — the baseline every later stage is read
% against. Three analyses:
%
%   1. Topoplots of the unperturbed field at chosen cord positions, for
%      Biot-Savart MSG, BEM MSG and BEM ESG side by side.
%   2. Source amplitude along the cord: peak and RMS field per source, and a
%      focality index (the fraction of sensors seeing at least half the peak;
%      smaller = more spatially concentrated field).
%   3. Forward-model-type decomposition: RE, gain, RDM and r2 along the cord
%      between model types on the same geometry (Biot-Savart vs BEM for MSG).
%      This sets the scale of "how different are two reasonable models",
%      which later stages compare the perturbation effects against.
%
% USAGE:
%   pt_baseline
%
% INPUT:
%   <forward_fields_base>/leadfields_organised.mat for each modality used in
%   baseline_models (config_pert)
%
% OUTPUT (to <stage_results_dir>/1_baseline/):
%   topoplots/topo_src<MM>mm_slot<K>.png/.fig   rows = models, cols = orientations
%   amplitude_profiles.png/.fig                 RMS field along cord, per model
%   amplitude_normalised.png/.fig               models overlaid, each / own max
%   field_focality.png/.fig                     focality index along cord
%   decomposition_<pair>.png/.fig               RE / gain / RDM / r2 along cord
%   baseline_amplitude.csv/.txt/.tex            amplitude summary per model
%   baseline_model_type.csv/.txt/.tex           model-type metrics with CIs
%   baseline_per_source.csv                     every per-source value
%   stage1_long.csv                             rows for pt_summary_table
%
% NOTES:
%   - MSG is in fT/nAm and ESG in uV/nAm, so absolute amplitudes are shown in
%     separate panels; only the normalised and focality figures overlay them.
%   - Model-type CIs resample SOURCE POSITIONS: they say how much the median
%     would move had the cord been sampled at different positions.
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

config_pert;
pt_add_functions;
rng(stats_seed);

fprintf('Stage 1 — noise-free baseline\n\n');

out_dir  = fullfile(stage_results_dir, '1_baseline');
topo_dir = fullfile(out_dir, 'topoplots');
if ~exist(topo_dir, 'dir'); mkdir(topo_dir); end

ref_geom = sensitivity_ref_key;           % unperturbed source geometry
oris3    = orientation_labels;            % VD RC LR (topoplots, amplitude)
mopts    = metric_defaults();


% =========================================================================
% LOAD
% =========================================================================

LF = struct();
mods_needed = unique([baseline_models(:, 1); baseline_pairs(:, 1)])';
for m = mods_needed
    f = fullfile(mods_cfg.(m{1}).forward_fields_base, 'leadfields_organised.mat');
    if ~isfile(f)
        warning('Stage 1: %s not found — %s models skipped.', f, upper(m{1}));
        continue
    end
    tmp = load(f, 'leadfields');
    LF.(m{1}) = tmp.leadfields;
end

n_mod  = size(baseline_models, 1);
models = struct('mod', {}, 'method', {}, 'label', {}, 'E', {}, 'color', {});
for i = 1:n_mod
    md  = baseline_models{i, 1};
    key = [baseline_models{i, 2} '_' ref_geom];
    if ~isfield(LF, md) || ~isfield(LF.(md), key)
        warning('Stage 1: %s not loaded — skipped.', key);
        continue
    end
    models(end+1) = struct('mod', md, 'method', baseline_models{i, 2}, ...
        'label', baseline_models{i, 3}, 'E', LF.(md).(key), ...
        'color', baseline_model_colors(i, :)); %#ok<SAGROW>
end
if isempty(models)
    error('Stage 1: none of baseline_models could be loaded.');
end

n_sources = models(1).E.n_sources;
src_range = 2:(n_sources - 1);
dist_mm   = src_range * src_spacing_mm;


% =========================================================================
% 1. TOPOPLOTS
% =========================================================================

fprintf('[1/3] Topoplots\n');
slot_names = {'X / Tangential', 'Y', 'Z / Radial'};

for src_mm = baseline_topo_src_mm
    si = min(max(round(src_mm / src_spacing_mm), 1), n_sources);
    for slot = 1:3
        fig = figure('Color', 'w', 'Units', 'inches', ...
            'Position', [1 1 3.2*numel(oris3) + 1, 2.6*numel(models)]);
        tl = tiledlayout(numel(models), numel(oris3) + 1, ...
            'TileSpacing', 'compact', 'Padding', 'compact');

        for r = 1:numel(models)
            M  = models(r);
            mc = mods_cfg.(M.mod);
            ax_idx = find(mc.axis_slot == slot, 1);

            nexttile(tl, (r-1)*(numel(oris3)+1) + 1);
            if isempty(ax_idx)
                row_lbl = sprintf('%s\n(no %s axis)', M.label, slot_names{slot});
            else
                row_lbl = sprintf('%s\n%s axis', M.label, mc.axis_names{ax_idx});
            end
            text(0.5, 0.5, row_lbl, 'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', 'FontSize', 10, 'Rotation', 90);
            axis off

            pos = [];
            if ~isempty(ax_idx)
                gfile = fullfile(mc.geoms_path, ['geometries_' ref_geom '.mat']);
                try
                    pba = sim_sensor_positions(gfile, analysis_array, ...
                        mc.sensor_is_meg, mc.sensor_n_axes);
                    pos = pba{ax_idx};
                    if size(pos, 1) ~= M.E.n_sensors_per_axis
                        warning('Stage 1: %s has %d sensor positions but %d leadfield rows — topoplot skipped.', ...
                            M.label, size(pos, 1), M.E.n_sensors_per_axis);
                        pos = [];
                    end
                catch err
                    warning('Stage 1: sensor positions for %s: %s', M.label, err.message);
                end
            end

            lim = 0;
            if ~isempty(pos)
                for o = 1:numel(oris3)
                    lim = max(lim, max(abs(M.E.(oris3{o}){ax_idx, si})));
                end
            end

            for o = 1:numel(oris3)
                nexttile(tl, (r-1)*(numel(oris3)+1) + 1 + o);
                if isempty(pos) || lim == 0
                    text(0.5, 0.5, 'n/a', 'HorizontalAlignment', 'center', ...
                        'Color', [0.5 0.5 0.5]);
                    axis off
                else
                    plot_topoplot_publication(pos, M.E.(oris3{o}){ax_idx, si}, ...
                        [-lim lim], mc.sensor_is_meg);
                end
                if r == 1
                    title(ori_titles.(oris3{o}), 'FontSize', 11);
                end
            end
        end
        title(tl, sprintf('Unperturbed field, source at %d mm — %s slot (%s array)', ...
            si * src_spacing_mm, slot_names{slot}, analysis_array), ...
            'FontSize', 13, 'FontWeight', 'bold');
        pt_save_figure(fig, topo_dir, sprintf('topo_src%dmm_slot%d', ...
            si * src_spacing_mm, slot), fig_resolution);
    end
end


% =========================================================================
% 2. SOURCE AMPLITUDE ALONG THE CORD
% =========================================================================

fprintf('[2/3] Source amplitude\n');

amp = struct();   % amp(r).rms / .peak / .focal : {n_cfg} of [n_ori3 x n_src]
per_source_rows = {};
amp_rows = {};

for r = 1:numel(models)
    E      = models(r).E;
    n_ax   = E.n_sensor_axes;
    cfgs   = [1:n_ax, 0];
    cnames = [mods_cfg.(models(r).mod).axis_names, {'Array'}];
    amp(r).cfg_names = cnames;
    for c = 1:numel(cfgs)
        R = nan(numel(oris3), numel(src_range));
        P = R; F = R;
        for o = 1:numel(oris3)
            L = pt_lf_matrix(E, oris3{o}, cfgs(c), E.n_sensors_per_axis, src_range);
            pk = max(abs(L), [], 1);
            R(o, :) = sqrt(mean(L.^2, 1));
            P(o, :) = pk;
            F(o, :) = mean(abs(L) >= 0.5 * pk, 1);
            for s = 1:numel(src_range)
                per_source_rows(end+1, :) = {models(r).label, models(r).mod, ...
                    models(r).method, oris3{o}, cnames{c}, dist_mm(s), ...
                    P(o, s), R(o, s), F(o, s)}; %#ok<SAGROW>
            end
            amp_rows(end+1, :) = {models(r).label, mods_cfg.(models(r).mod).unit, ...
                oris3{o}, cnames{c}, median(R(o, :)), min(R(o, :)), max(R(o, :)), ...
                median(P(o, :)), median(F(o, :))}; %#ok<SAGROW>
        end
        amp(r).rms{c}   = R;
        amp(r).peak{c}  = P;
        amp(r).focal{c} = F;
    end
end

% 2a. Absolute RMS: one row per model (its own unit), one column per orientation
fig = figure('Color', 'w', 'Position', [50 50 1400 330*numel(models)]);
tl  = tiledlayout(numel(models), numel(oris3), 'TileSpacing', 'compact');
for r = 1:numel(models)
    n_c  = numel(amp(r).cfg_names);
    cols = lines(n_c);
    cols(end, :) = [0 0 0];                     % whole array in black
    for o = 1:numel(oris3)
        ax = nexttile(tl); hold(ax, 'on');
        for c = 1:n_c
            plot(ax, dist_mm, amp(r).rms{c}(o, :), 'Color', cols(c, :), ...
                'LineWidth', 1.6 + 0.6*(c == n_c), 'DisplayName', amp(r).cfg_names{c});
        end
        grid(ax, 'on'); set(ax, 'TickDir', 'out');
        if r == 1, title(ax, ori_titles.(oris3{o})); end
        if o == 1
            ylabel(ax, sprintf('%s\nRMS (%s)', models(r).label, ...
                mods_cfg.(models(r).mod).unit));
        end
        if r == numel(models), xlabel(ax, 'Distance along cord (mm)'); end
        if o == numel(oris3), legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 8); end
    end
end
title(tl, 'Unperturbed field amplitude along the cord (RMS over sensors)', ...
    'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'amplitude_profiles', fig_resolution);

% 2b/2c. Models overlaid: whole-array RMS normalised to own maximum, and focality
for what = {'amplitude_normalised', 'field_focality'}
    fig = figure('Color', 'w', 'Position', [50 50 1400 380]);
    tl  = tiledlayout(1, numel(oris3), 'TileSpacing', 'compact');
    for o = 1:numel(oris3)
        ax = nexttile(tl); hold(ax, 'on');
        for r = 1:numel(models)
            if strcmp(what{1}, 'amplitude_normalised')
                y = amp(r).rms{end}(o, :) / max(amp(r).rms{end}(:));
            else
                y = amp(r).focal{end}(o, :);
            end
            plot(ax, dist_mm, y, 'Color', models(r).color, 'LineWidth', 1.8, ...
                'DisplayName', models(r).label);
        end
        grid(ax, 'on'); set(ax, 'TickDir', 'out');
        title(ax, ori_titles.(oris3{o}));
        xlabel(ax, 'Distance along cord (mm)');
        if o == 1
            if strcmp(what{1}, 'amplitude_normalised')
                ylabel(ax, 'RMS / model maximum');
            else
                ylabel(ax, 'Fraction of sensors \geq half peak');
            end
        end
        if o == numel(oris3), legend(ax, 'Location', 'best', 'Box', 'off'); end
    end
    if strcmp(what{1}, 'amplitude_normalised')
        title(tl, 'Relative amplitude along the cord (whole array, each model / its own maximum)', 'FontWeight', 'bold');
    else
        title(tl, 'Field focality (whole array; smaller = more concentrated field)', 'FontWeight', 'bold');
    end
    pt_save_figure(fig, out_dir, what{1}, fig_resolution);
end

Ta = cell2table(amp_rows, 'VariableNames', {'model', 'unit', 'orientation', ...
    'axis_set', 'rms_median', 'rms_min', 'rms_max', 'peak_median', 'focality_median'});
pt_write_table(Ta, fullfile(out_dir, 'baseline_amplitude'), struct( ...
    'cols', {{ 'model', 'Model', 'Model', '%s'; ...
               'orientation', 'Ori', 'Ori.', '%s'; ...
               'axis_set', 'Axes', 'Axes', '%s'; ...
               'unit', 'Unit', 'Unit', '%s'; ...
               'rms_median', 'RMS med', 'RMS (median)', '%.3g'; ...
               'rms_min', 'RMS min', 'RMS (min)', '%.3g'; ...
               'rms_max', 'RMS max', 'RMS (max)', '%.3g'; ...
               'peak_median', 'Peak med', 'Peak (median)', '%.3g'; ...
               'focality_median', 'Focality', 'Focality', '%.3f' }}, ...
    'title', 'STAGE 1 — UNPERTURBED FIELD AMPLITUDE', ...
    'notes', {{'Values over cord positions (end sources excluded).', ...
               'Focality = fraction of sensors with |L| >= half the peak; smaller is more concentrated.'}}, ...
    'caption', ['Amplitude of the unperturbed forward field along the cord. ' ...
                'RMS and peak are over sensors, then summarised over cord ' ...
                'positions. Focality is the fraction of sensors reaching half ' ...
                'the peak value.'], ...
    'label', 'tab:baseline-amplitude', 'group_by', 'model', ...
    'tex_rows', strcmp(Ta.axis_set, 'Array')));

Tps = cell2table(per_source_rows, 'VariableNames', {'model', 'modality', ...
    'method', 'orientation', 'axis_set', 'distance_mm', 'peak', 'rms', 'focality'});


% =========================================================================
% 3. FORWARD-MODEL-TYPE DECOMPOSITION
% =========================================================================

fprintf('[3/3] Forward-model-type decomposition\n');

long_rows = pt_long_row([]);
mt_rows   = {};

for p = 1:size(baseline_pairs, 1)
    md = baseline_pairs{p, 1};
    ka = [baseline_pairs{p, 2} '_' ref_geom];
    kb = [baseline_pairs{p, 3} '_' ref_geom];
    lbl = baseline_pairs{p, 4};
    if ~isfield(LF, md) || ~isfield(LF.(md), ka) || ~isfield(LF.(md), kb)
        warning('Stage 1: pair "%s" not available — skipped.', lbl);
        continue
    end
    EA = LF.(md).(ka); EB = LF.(md).(kb);
    n_ax   = EA.n_sensor_axes;
    cfgs   = [1:n_ax, 0];
    cnames = [mods_cfg.(md).axis_names, {'Array'}];
    n_tr   = min(EA.n_sensors_per_axis, EB.n_sensors_per_axis);

    S = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
    for c = 1:numel(cfgs)
        S(c).label = sprintf('%s axes', cnames{c});
        for o = 1:numel(metric_orientations)
            ori = metric_orientations{o};
            A = pt_lf_matrix(EA, ori, cfgs(c), n_tr, src_range);
            B = pt_lf_matrix(EB, ori, cfgs(c), n_tr, src_range);
            Mo = pt_metrics_block(A, B, mopts);
            gain = (exp(Mo.lnmag) - 1) * 100;
            S(c).re(o, :)   = Mo.re;
            S(c).gain(o, :) = gain;
            S(c).rdm(o, :)  = Mo.rdm;
            S(c).rsq(o, :)  = Mo.rsq;

            D = pt_describe(Mo.re, stats_n_boot);
            mt_rows(end+1, :) = {lbl, mods_cfg.(md).display, ori, cnames{c}, ...
                numel(Mo.re), D.median, D.ci_lo, D.ci_hi, max(Mo.re), ...
                median(Mo.rsq, 'omitnan'), min(Mo.rsq), median(Mo.rdm, 'omitnan'), ...
                median(gain, 'omitnan'), median(abs(gain), 'omitnan'), ...
                median(sqrt(gain.^2 + (Mo.rdm*100).^2), 'omitnan')}; %#ok<SAGROW>

            meta = struct('stage', '1_baseline', 'modality', md, 'system', 'noise-free', ...
                'method', [baseline_pairs{p, 3} ' vs ' baseline_pairs{p, 2}], ...
                'factor', 'model_type', 'level', lbl, 'orientation', ori, ...
                'axis_set', cnames{c}, 'unit', 'source');
            long_rows = [long_rows; pt_long_row(meta, struct('re', Mo.re, ...
                'rsq', Mo.rsq, 'rdm', Mo.rdm, 'gain', gain, 'absgain', abs(gain)), ...
                stats_n_boot)]; %#ok<AGROW>
        end
    end

    cols = lines(numel(cfgs)); cols(end, :) = [0 0 0];
    plot_metric_decomposition(S, struct( ...
        'dist', dist_mm, 'orientation_labels', {metric_orientations}, ...
        'ori_titles', ori_titles, 'colors', cols, ...
        'title', sprintf('%s — %s, unperturbed geometry (reference: %s)', ...
            lbl, mods_cfg.(md).display, upper(baseline_pairs{p, 2})), ...
        'save_dir', out_dir, ...
        'save_name', ['decomposition_' matlab.lang.makeValidName(lbl)]));
end

if ~isempty(mt_rows)
    Tm = cell2table(mt_rows, 'VariableNames', {'comparison', 'modality', ...
        'orientation', 'axis_set', 'n_sources', 're', 're_ci_lo', 're_ci_hi', ...
        're_max', 'rsq_median', 'rsq_min', 'rdm_median', 'gain_median', ...
        'absgain_median', 're_from_decomposition'});
    pt_write_table(Tm, fullfile(out_dir, 'baseline_model_type'), struct( ...
        'cols', {{ 'comparison', 'Comparison', 'Comparison', '%s'; ...
                   'orientation', 'Ori', 'Ori.', '%s'; ...
                   'axis_set', 'Axes', 'Axes', '%s'; ...
                   're', 'RE % [95% CI]', 'RE (\%) [95\% CI]', 'ci'; ...
                   're_max', 'RE max', 'RE max', '%.2f'; ...
                   'rsq_median', 'r2', '$r^2$', '%.4f'; ...
                   'rdm_median', 'RDM', 'RDM', '%.4f'; ...
                   'gain_median', 'gain %', 'Gain (\%)', '%+.2f'; ...
                   're_from_decomposition', 'sqrt(g2+RDM2)', '$\sqrt{g^2+\mathrm{RDM}^2}$', '%.2f' }}, ...
        'title', 'STAGE 1 — FORWARD-MODEL-TYPE COMPARISON (unperturbed, noise-free)', ...
        'notes', {{'RE is asymmetric: the reference model is the denominator.', ...
                   'Medians over cord positions; CI resamples source positions.', ...
                   'RE ~ sqrt(gain^2 + (RDM*100)^2): compare the last column with RE.'}}, ...
        'caption', ['Difference between forward model types on the unperturbed ' ...
                    'geometry. RE is relative to the reference model; $r^2$ and ' ...
                    'RDM describe field shape only, gain describes amplitude ' ...
                    'only. Medians over cord positions with 95\% bootstrap CIs ' ...
                    'over source positions.'], ...
        'label', 'tab:baseline-model-type', 'group_by', 'axis_set'));
end

writetable(Tps, fullfile(out_dir, 'baseline_per_source.csv'));
writetable(long_rows, fullfile(out_dir, 'stage1_long.csv'));

fprintf('\nStage 1 complete: %s\n', out_dir);

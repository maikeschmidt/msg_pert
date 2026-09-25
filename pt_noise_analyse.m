% pt_noise_analyse - Stage 5b: how sensor noise changes the perturbation results
%
% Reads the simulations of pt_noise_simulate and repeats the questions of
% stages 2-4 on the NOISY measurements, at every noise level:
%
%   A. Noise alone          how far does noise move the measured field from
%                           the true one, with a perfect model?
%   B. Perturbation + noise the measured error when the geometry was wrong,
%                           against the noise-only error at the same level
%   C. Detectability        in what fraction of realisations does each
%                           perturbation push the error beyond the range noise
%                           alone produces, and at what noise level does it
%                           stop doing so (the critical noise level)?
%   D. Error budget         share of the squared error due to the
%                           perturbation rather than the noise
%   E. Systems compared     SQUID MSG vs ESG, OP-MSG vs ESG, SQUID vs OP-MSG,
%                           per perturbation and bundle, at every level: do
%                           the stage-4 differences survive realistic noise?
%   F. Within-modality      the stage 2/3 magnitude and type tests, repeated
%      tests repeated       at the reference noise level (noise_reference_factor)
%
% STATISTICAL UNIT
%   As in stages 2-4: one perturbation realisation, summarised by its cord
%   median, here averaged over the noise realisations. Comparisons between
%   two MSG systems are always paired (same geometry, same shifts); MSG vs
%   ESG follows msg_esg_paired. System comparisons use the whole array.
%
% USAGE:
%   pt_noise_analyse
%
% OUTPUT (to <stage_results_dir>/5_noise/):
%   noise_only_<system>.png             A: metric vs noise level
%   noise_pert_<system>_<type>.png      B: pert + noise vs noise alone
%   detectability_<system>.png          C: detection rate vs noise level
%   error_budget.png                    D: perturbation share of squared RE
%   systems_<a>_vs_<b>.png              E: effect size vs noise level
%   noise_only.csv/.txt/.tex            A (with amplitude SNR)
%   noise_combined.csv/.txt/.tex        B, C, D per type x bundle x level
%   noise_critical_levels.csv/.txt/.tex C
%   noise_system_stats.csv/.txt/.tex    E
%   noise_within_stats.csv/.txt/.tex    F
%   stage5_long.csv                     rows for pt_summary_table
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

pt_modality('clear');
config_pert;
pt_add_functions;
rng(stats_seed);

fprintf('Stage 5b — analysis of noisy measurements\n\n');

out_dir  = fullfile(stage_results_dir, '5_noise');
data_dir = fullfile(out_dir, 'data');

SYS = struct();
for k = 1:numel(noise_systems)
    f = fullfile(data_dir, sprintf('noise_%s.mat', noise_systems(k).short));
    if isfile(f)
        tmp = load(f, 'NS');
        SYS.(noise_systems(k).short) = tmp.NS;
    else
        warning('%s not found — %s skipped. Run pt_noise_simulate.', f, noise_systems(k).label);
    end
end
sys_ids = fieldnames(SYS)';
if isempty(sys_ids), error('No noise simulations found in %s.', data_dir); end

info   = pt_metric_info();
tested = info([info.tested]);
nb     = numel(bundle_names);
levels = SYS.(sys_ids{1}).levels;
n_lev  = numel(levels);
L_ref  = find(abs(levels - noise_reference_factor) < 1e-12, 1);
if isempty(L_ref)
    error('noise_reference_factor = %g is not one of noise_factors.', noise_reference_factor);
end
lev_lbl = arrayfun(@(x) sprintf('%gx', x), levels, 'UniformOutput', false);
xpos    = 1:n_lev;                     % categorical positions (level 0 has no log)
plot_metrics = {'re', 'absgain', 'rdm', 'rsq'};
mlabel  = @(nm) info(strcmp({info.name}, nm)).label;
mtxt    = @(nm) info(strcmp({info.name}, nm)).txt;
bcol    = struct('source', source_bundle_colors, 'sensor', sensor_bundle_colors, ...
                 'cond', cond_bundle_colors);

% Per-realisation value for one type/metric/axis set/orientation at level L:
% mean over noise realisations of the cord median. Row 1 is noise only.
shiftval = @(NS, ty, nm, L, c, o) mean(NS.pert.(ty).cm.(nm)(:, L, :, c, o), 3, 'omitnan');

long_rows = pt_long_row([]);


% =========================================================================
% A. NOISE ALONE
% =========================================================================

fprintf('[A] Noise alone\n');
no_rows = {};
for s = sys_ids
    NS = SYS.(s{1}); ty0 = first_type(NS, pert_types);
    n_cfg = numel(NS.axcfg_names); n_ori = numel(NS.ori_labels);

    fig = figure('Color', 'w', 'Position', [40 40 380*n_ori, 1050]);
    tl  = tiledlayout(numel(plot_metrics), n_ori, 'TileSpacing', 'compact');
    cols = lines(n_cfg); cols(end, :) = [0 0 0];
    for r = 1:numel(plot_metrics)
        for o = 1:n_ori
            ax = nexttile(tl); hold(ax, 'on');
            for c = 1:n_cfg
                Y = squeeze(NS.pert.(ty0).cm.(plot_metrics{r})(1, :, :, c, o))';   % [n_real x n_lev]
                pt_plot_band(ax, xpos, Y, cols(c, :), '-o', NS.axcfg_names{c});
            end
            style_level_axis(ax, xpos, lev_lbl, r == numel(plot_metrics));
            if r == 1, title(ax, ori_titles.(NS.ori_labels{o})); end
            if o == 1, ylabel(ax, mlabel(plot_metrics{r})); end
            if r == 1 && o == n_ori, legend(ax, 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, sprintf(['%s — measured vs true field, unperturbed model (noise only)\n' ...
        'median and IQR over noise realisations of the cord median'], NS.system.label), ...
        'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['noise_only_' s{1}], fig_resolution);

    for c = 1:n_cfg
        for o = 1:n_ori
            for L = 1:n_lev
                x = struct();
                for k = 1:numel(info)
                    nm = info(k).name;
                    if strcmp(nm, 'gain'), src = 'lnmag'; else, src = nm; end
                    v = squeeze(NS.pert.(ty0).cm.(src)(1, L, :, c, o));
                    if strcmp(nm, 'gain'), v = (exp(v) - 1) * 100; end
                    x.(nm) = v;
                end
                no_rows(end+1, :) = {NS.system.label, levels(L), NS.ori_labels{o}, ...
                    NS.axcfg_names{c}, NS.snr(L, c, o), median(x.re), ...
                    median(x.rsq), median(x.rdm), median(x.absgain)}; %#ok<SAGROW>
                meta = struct('stage', '5_noise', 'modality', NS.modality, ...
                    'system', sprintf('%s @ %s', NS.system.label, lev_lbl{L}), ...
                    'method', NS.system.method, 'factor', 'noise_only', ...
                    'level', lev_lbl{L}, 'orientation', NS.ori_labels{o}, ...
                    'axis_set', NS.axcfg_names{c}, 'unit', 'noise realisation');
                long_rows = [long_rows; pt_long_row(meta, x, stats_n_boot)]; %#ok<AGROW>
            end
        end
    end
end
Tno = cell2table(no_rows, 'VariableNames', {'system', 'noise_factor', 'orientation', ...
    'axis_set', 'snr_median', 're_median', 'rsq_median', 'rdm_median', 'absgain_median'});
pt_write_table(Tno, fullfile(out_dir, 'noise_only'), struct( ...
    'cols', {{ 'noise_factor', 'x floor', '$\times$ floor', '%g'; ...
               'orientation', 'Ori', 'Ori.', '%s'; ...
               'axis_set', 'Axes', 'Axes', '%s'; ...
               'snr_median', 'SNR', 'SNR', '%.3g'; ...
               're_median', 'RE %', 'RE (\%)', '%.2f'; ...
               'rsq_median', 'r2', '$r^2$', '%.4f'; ...
               'rdm_median', 'RDM', 'RDM', '%.4f'; ...
               'absgain_median', '|gain| %', '$|$Gain$|$ (\%)', '%.2f' }}, ...
    'title', 'STAGE 5 — NOISE ALONE (unperturbed model, measured vs true field)', ...
    'notes', {{'SNR = median over cord of rms(true field) / s.d. of the noise on the recovered field.', ...
               'Metrics: median over noise realisations of the cord median.'}}, ...
    'caption', ['Effect of sensor noise alone: the recovered field pattern of the ' ...
               'unperturbed model compared with the true field, per system and ' ...
               'noise level (multiples of each system''s own noise floor).'], ...
    'label', 'tab:noise-only', 'group_by', 'system', ...
    'tex_rows', strcmp(Tno.orientation, headline_orientation) & strcmp(Tno.axis_set, 'Array')));


% =========================================================================
% B, C, D. PERTURBATION + NOISE, DETECTABILITY, ERROR BUDGET
% =========================================================================

fprintf('[B/C/D] Perturbation + noise, detectability, error budget\n');
cb_rows = {};
crit_rows = {};
det_metrics = {'re', 'rdm', 'rsq'};

for s = sys_ids
    NS = SYS.(s{1});
    c = numel(NS.axcfg_names);                       % whole array
    o = find(strcmp(NS.ori_labels, headline_orientation));
    types = pert_types(isfield(NS.pert, pert_types));

    for t = 1:numel(types)
        ty = types{t}; P = NS.pert.(ty);

        % B. figure: rows = metrics, cols = bundles
        fig = figure('Color', 'w', 'Position', [40 40 1400 1050]);
        tl  = tiledlayout(numel(plot_metrics), nb, 'TileSpacing', 'compact');
        for r = 1:numel(plot_metrics)
            for b = 1:nb
                ax = nexttile(tl); hold(ax, 'on');
                Y0 = squeeze(P.cm.(plot_metrics{r})(1, :, :, c, o))';
                pt_plot_band(ax, xpos, Y0, [0.5 0.5 0.5], '--', 'noise only');
                V  = squeeze(mean(P.cm.(plot_metrics{r})(2:end, :, :, c, o), 3, 'omitnan'));
                pt_plot_band(ax, xpos, V(P.bundle_idx == b, :), bcol.(ty)(b, :), '-o', ...
                    sprintf('%s + noise', bundle_names{b}));
                style_level_axis(ax, xpos, lev_lbl, r == numel(plot_metrics));
                if r == 1, title(ax, sprintf('%s bundle', bundle_names{b})); end
                if b == 1, ylabel(ax, mlabel(plot_metrics{r})); end
                if r == 1, legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 7); end
            end
        end
        title(tl, sprintf(['%s — %s perturbation measured through noise (%s, whole array)\n' ...
            'grey = unperturbed model + noise; colour = perturbed + noise (IQR over realisations)'], ...
            NS.system.label, lower(pert_type_display.(ty)), ori_titles.(headline_orientation)), ...
            'FontWeight', 'bold');
        pt_save_figure(fig, out_dir, sprintf('noise_pert_%s_%s', s{1}, ty), fig_resolution);

        % C + D numbers, every orientation and axis set in the CSV
        for cc = 1:numel(NS.axcfg_names)
            for oo = 1:numel(NS.ori_labels)
                for b = 1:nb
                    sel = find(P.bundle_idx == b) + 1;     % +1: row 1 is noise only
                    rate = struct();
                    for dm = det_metrics
                        wr  = info(strcmp({info.name}, dm{1})).worse;
                        rate.(dm{1}) = nan(1, n_lev);
                        for L = 1:n_lev
                            null = squeeze(P.cm.(dm{1})(1, L, :, cc, oo));
                            X    = squeeze(P.cm.(dm{1})(sel, L, :, cc, oo));
                            if wr > 0
                                thr = pt_quantile(null, noise_detect_pct, 1);
                                det = X > thr;
                            else
                                thr = pt_quantile(null, 100 - noise_detect_pct, 1);
                                det = X < thr;
                            end
                            if levels(L) == 0
                                % Noise-free: any real difference is detected
                                det = abs(X - null(1)) > 1e-12;
                            end
                            det(isnan(X)) = false;
                            rate.(dm{1})(L) = mean(det(:));
                        end
                        [crit, crit_txt] = critical_level(levels, rate.(dm{1}), noise_detect_rate);
                        if cc == numel(NS.axcfg_names)
                            crit_rows(end+1, :) = {NS.system.label, pert_type_display.(ty), ...
                                bundle_names{b}, NS.ori_labels{oo}, mtxt(dm{1}), crit, crit_txt}; %#ok<SAGROW>
                        end
                    end
                    for L = 1:n_lev
                        re_i  = squeeze(P.cm.re(sel, L, :, cc, oo));                 % [n_b x n_real]
                        re_0  = squeeze(P.cm.re(sel, 1, 1, cc, oo));                 % noise-free
                        share = median(re_0.^2 ./ mean(re_i.^2, 2, 'omitnan'), 'omitnan') * 100;
                        x = struct();
                        for k = 1:numel(info)
                            nm = info(k).name;
                            if strcmp(nm, 'gain'), src = 'lnmag'; else, src = nm; end
                            v = mean(P.cm.(src)(sel, L, :, cc, oo), 3, 'omitnan');
                            if strcmp(nm, 'gain'), v = (exp(v) - 1) * 100; end
                            x.(nm) = v(:);
                        end
                        D  = pt_describe(x.re, stats_n_boot);
                        re_noise = median(squeeze(P.cm.re(1, L, :, cc, oo)), 'omitnan');
                        cb_rows(end+1, :) = {NS.system.label, pert_type_display.(ty), ...
                            bundle_names{b}, levels(L), NS.ori_labels{oo}, NS.axcfg_names{cc}, ...
                            re_noise, D.median, D.ci_lo, D.ci_hi, median(x.rsq, 'omitnan'), ...
                            median(x.rdm, 'omitnan'), median(x.absgain, 'omitnan'), share, ...
                            rate.re(L), rate.rdm(L), rate.rsq(L)}; %#ok<SAGROW>
                        meta = struct('stage', '5_noise', 'modality', NS.modality, ...
                            'system', sprintf('%s @ %s', NS.system.label, lev_lbl{L}), ...
                            'method', NS.system.method, 'factor', ty, ...
                            'level', bundle_names{b}, 'orientation', NS.ori_labels{oo}, ...
                            'axis_set', NS.axcfg_names{cc}, 'unit', 'realisation');
                        long_rows = [long_rows; pt_long_row(meta, x, stats_n_boot)]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % C. detectability figure: rows = detection metrics, cols = types
    fig = figure('Color', 'w', 'Position', [40 40 400*numel(types) + 100, 280*numel(det_metrics) + 120]);
    tl  = tiledlayout(numel(det_metrics), numel(types), 'TileSpacing', 'compact');
    Tcb_sys = cell2table(cb_rows(strcmp(cb_rows(:, 1), NS.system.label), :), 'VariableNames', ...
        {'system', 'perturbation', 'bundle', 'noise_factor', 'orientation', 'axis_set', ...
         're_noise_only', 're', 're_ci_lo', 're_ci_hi', 'rsq', 'rdm', 'absgain', ...
         'pert_share_pct', 'detect_re', 'detect_rdm', 'detect_rsq'});
    for d = 1:numel(det_metrics)
        for t = 1:numel(types)
            ax = nexttile(tl); hold(ax, 'on');
            for b = 1:nb
                sel = strcmp(Tcb_sys.perturbation, pert_type_display.(types{t})) & ...
                      strcmp(Tcb_sys.bundle, bundle_names{b}) & ...
                      strcmp(Tcb_sys.orientation, headline_orientation) & ...
                      strcmp(Tcb_sys.axis_set, 'Array');
                y = Tcb_sys.(['detect_' det_metrics{d}])(sel);
                plot(ax, xpos(1:numel(y)), y, '-o', 'Color', bcol.(types{t})(b, :), ...
                    'MarkerFaceColor', bcol.(types{t})(b, :), 'LineWidth', 1.6, ...
                    'DisplayName', bundle_names{b});
            end
            yline(ax, noise_detect_rate, ':k', 'HandleVisibility', 'off');
            ylim(ax, [-0.02 1.02]);
            style_level_axis(ax, xpos, lev_lbl, d == numel(det_metrics));
            if d == 1, title(ax, pert_type_display.(types{t})); end
            if t == 1, ylabel(ax, sprintf('Detection rate (%s)', mlabel(det_metrics{d}))); end
            if d == 1 && t == numel(types), legend(ax, 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, sprintf(['%s — is the perturbation visible above the noise? (%s, whole array)\n' ...
        'fraction of realisations beyond the %gth percentile of the noise-only error'], ...
        NS.system.label, ori_titles.(headline_orientation), noise_detect_pct), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['detectability_' s{1}], fig_resolution);
end

Tcb = cell2table(cb_rows, 'VariableNames', {'system', 'perturbation', 'bundle', ...
    'noise_factor', 'orientation', 'axis_set', 're_noise_only', 're', 're_ci_lo', ...
    're_ci_hi', 'rsq', 'rdm', 'absgain', 'pert_share_pct', 'detect_re', ...
    'detect_rdm', 'detect_rsq'});
Tcb.group = strcat(Tcb.system, {' — '}, Tcb.perturbation);
pt_write_table(Tcb, fullfile(out_dir, 'noise_combined'), struct( ...
    'cols', {{ 'bundle', 'Bundle', 'Bundle', '%s'; ...
               'noise_factor', 'x floor', '$\times$ floor', '%g'; ...
               're_noise_only', 'RE noise', 'RE noise only', '%.2f'; ...
               're', 'RE pert+noise [CI]', 'RE pert.+noise [95\% CI]', 'ci'; ...
               'rsq', 'r2', '$r^2$', '%.4f'; ...
               'rdm', 'RDM', 'RDM', '%.4f'; ...
               'pert_share_pct', 'pert share %', 'Pert. share (\%)', '%.1f'; ...
               'detect_re', 'detect RE', 'Detect (RE)', '%.2f'; ...
               'detect_rdm', 'detect RDM', 'Detect (RDM)', '%.2f' }}, ...
    'title', 'STAGE 5 — PERTURBATION MEASURED THROUGH NOISE', ...
    'notes', {{'RE noise = unperturbed model + noise vs true field; RE pert+noise = perturbed + noise vs true field.', ...
               'pert share = perturbation part of the squared RE: RE_noise-free^2 / mean(RE_measured^2).', ...
               sprintf('detect = fraction of realisations beyond the %gth percentile of the noise-only error.', noise_detect_pct)}}, ...
    'caption', ['Perturbation effects as measured through sensor noise. The ' ...
               'perturbed field plus noise is compared with the true, unperturbed ' ...
               'field. The perturbation share is the fraction of the squared ' ...
               'relative error attributable to the perturbation; detection is ' ...
               'the fraction of realisations whose error exceeds the noise-only ' ...
               'range.'], ...
    'label', 'tab:noise-combined', 'group_by', 'group', ...
    'tex_rows', strcmp(Tcb.orientation, headline_orientation) & strcmp(Tcb.axis_set, 'Array') & ...
                (Tcb.noise_factor == 0 | abs(Tcb.noise_factor - noise_reference_factor) < 1e-12 | ...
                 Tcb.noise_factor == max(levels))));

Tcr = cell2table(crit_rows, 'VariableNames', {'system', 'perturbation', 'bundle', ...
    'orientation', 'metric', 'critical_factor', 'critical_text'});
Tcr.group = strcat(Tcr.system, {' — '}, Tcr.metric);
pt_write_table(Tcr, fullfile(out_dir, 'noise_critical_levels'), struct( ...
    'cols', {{ 'perturbation', 'Perturbation', 'Perturbation', '%s'; ...
               'bundle', 'Bundle', 'Bundle', '%s'; ...
               'orientation', 'Ori', 'Ori.', '%s'; ...
               'critical_text', 'critical x floor', 'Critical $\times$ floor', '%s' }}, ...
    'title', 'STAGE 5 — CRITICAL NOISE LEVEL (whole array)', ...
    'notes', {{sprintf('Highest noise level (x each system''s floor) at which >= %.0f%% of realisations are still detected; log-interpolated.', noise_detect_rate*100), ...
               '">" = still detected at the largest level simulated; "<" = not detected even at the smallest non-zero level.'}}, ...
    'caption', sprintf(['Critical noise level: the largest noise level, as a multiple ' ...
               'of each system''s own floor, at which at least %.0f\\%% of ' ...
               'perturbation realisations remain distinguishable from noise alone.'], ...
               noise_detect_rate*100), ...
    'label', 'tab:noise-critical', 'group_by', 'group', ...
    'tex_rows', strcmp(Tcr.orientation, headline_orientation)));

% D. error budget figure: pert share vs noise level, systems overlaid
types_all = pert_types(ismember(pert_types, fieldnames(SYS.(sys_ids{1}).pert)));
fig = figure('Color', 'w', 'Position', [40 40 400*numel(types_all) + 100, 300*nb + 120]);
tl  = tiledlayout(nb, numel(types_all), 'TileSpacing', 'compact');
for b = 1:nb
    for t = 1:numel(types_all)
        ax = nexttile(tl); hold(ax, 'on');
        for s = sys_ids
            sel = strcmp(Tcb.system, SYS.(s{1}).system.label) & ...
                  strcmp(Tcb.perturbation, pert_type_display.(types_all{t})) & ...
                  strcmp(Tcb.bundle, bundle_names{b}) & ...
                  strcmp(Tcb.orientation, headline_orientation) & strcmp(Tcb.axis_set, 'Array');
            if ~any(sel), continue; end
            plot(ax, xpos(1:sum(sel)), Tcb.pert_share_pct(sel), '-o', ...
                'Color', SYS.(s{1}).system.color, 'MarkerFaceColor', SYS.(s{1}).system.color, ...
                'LineWidth', 1.6, 'DisplayName', SYS.(s{1}).system.label);
        end
        ylim(ax, [0 105]);
        style_level_axis(ax, xpos, lev_lbl, b == nb);
        if b == 1, title(ax, pert_type_display.(types_all{t})); end
        if t == 1, ylabel(ax, sprintf('%s bundle\nPerturbation share (%%)', bundle_names{b})); end
        if b == 1 && t == numel(types_all), legend(ax, 'Location', 'best', 'Box', 'off'); end
    end
end
title(tl, sprintf(['Error budget — share of the squared RE due to the perturbation ' ...
    '(%s, whole array)\n100%% = noise negligible; 0%% = perturbation invisible'], ...
    ori_titles.(headline_orientation)), 'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'error_budget', fig_resolution);


% =========================================================================
% E. SYSTEMS COMPARED AT EVERY NOISE LEVEL
% =========================================================================

fprintf('[E] Systems compared across noise levels\n');
sopts = struct('bundle_names', {bundle_names}, 'n_perm', stats_n_perm);
tabs = {};
for q = 1:size(noise_compare_pairs, 1)
    ia = noise_compare_pairs{q, 1}; ib = noise_compare_pairs{q, 2};
    if ~isfield(SYS, ia) || ~isfield(SYS, ib), continue; end
    NA = SYS.(ia); NB = SYS.(ib);
    cmp = sprintf('%s vs %s', NA.system.label, NB.system.label);
    types = pert_types(isfield(NA.pert, pert_types) & isfield(NB.pert, pert_types));
    ca = numel(NA.axcfg_names); cb = numel(NB.axcfg_names);   % whole array
    for t = 1:numel(types)
        ty = types{t};
        if strcmp(NA.modality, NB.modality), paired = true;
        else, paired = msg_esg_paired.(ty); end
        for o = 1:numel(NA.ori_labels)
            ob = find(strcmp(NB.ori_labels, NA.ori_labels{o}));
            for L = 1:n_lev
                for k = 1:numel(tested)
                    va = shiftval(NA, ty, tested(k).name, L, ca, o);
                    vb = shiftval(NB, ty, tested(k).name, L, cb, ob);
                    T = pt_stats_between(va(2:end), vb(2:end), NA.pert.(ty).bundle_idx, ...
                        NB.pert.(ty).bundle_idx, paired, sopts, cmp);
                    T.perturbation = repmat({pert_type_display.(ty)}, height(T), 1);
                    T.noise_factor = repmat(levels(L), height(T), 1);
                    T.orientation  = repmat(NA.ori_labels(o), height(T), 1);
                    T.metric       = repmat({tested(k).txt}, height(T), 1);
                    tabs{end+1} = T; %#ok<SAGROW>
                end
            end
        end
    end
end
Tsys = table();
if ~isempty(tabs)
    Tsys = vertcat(tabs{:});
    Tsys = pt_fdr_by_group(Tsys, {'comparison', 'metric'}, stats_alpha);
    Tsys = movevars(Tsys, {'perturbation', 'metric', 'noise_factor', 'orientation'}, 'Before', 1);
    Tsys.group = strcat(Tsys.comparison, {' — '}, Tsys.metric);
    pt_write_table(Tsys, fullfile(out_dir, 'noise_system_stats'), struct( ...
        'cols', {{ 'perturbation', 'Perturbation', 'Perturbation', '%s'; ...
                   'subset', 'Bundle', 'Bundle', '%s'; ...
                   'noise_factor', 'x floor', '$\times$ floor', '%g'; ...
                   'median_a', 'A', 'A', '%.3g'; ...
                   'median_b', 'B', 'B', '%.3g'; ...
                   'effect', 'effect', 'Effect', '%+.2f'; ...
                   'p_fdr', 'p (FDR)', '$p_\mathrm{FDR}$', 'p' }}, ...
        'title', 'STAGE 5 — SENSOR SYSTEMS COMPARED UNDER NOISE (whole array)', ...
        'notes', {{'Row 0x floor = the noise-free result, so each block shows how a stage-4 difference changes as noise grows.', ...
                   'Positive effect = system A shows the larger value of the metric.', ...
                   'p_fdr: Benjamini-Hochberg within comparison x metric, across types, bundles, levels, orientations.'}}, ...
        'caption', ['Sensor systems compared on the same perturbations at each ' ...
                   'noise level. Each observation is one perturbation realisation ' ...
                   '(cord median, averaged over noise realisations).'], ...
        'label', 'tab:noise-systems', 'group_by', 'group', ...
        'tex_rows', strcmp(Tsys.orientation, headline_orientation) & ...
                    strcmp(Tsys.subset, 'pooled') & ...
                    (Tsys.noise_factor == 0 | abs(Tsys.noise_factor - noise_reference_factor) < 1e-12 | ...
                     Tsys.noise_factor == max(levels))));

    % Figure per pair: effect size vs noise level, filled = FDR-significant
    for q = 1:size(noise_compare_pairs, 1)
        ia = noise_compare_pairs{q, 1}; ib = noise_compare_pairs{q, 2};
        if ~isfield(SYS, ia) || ~isfield(SYS, ib), continue; end
        cmp = sprintf('%s vs %s', SYS.(ia).system.label, SYS.(ib).system.label);
        Tq  = Tsys(strcmp(Tsys.comparison, cmp) & strcmp(Tsys.orientation, headline_orientation), :);
        types = unique(Tq.perturbation, 'stable');
        fig = figure('Color', 'w', 'Position', [40 40 400*numel(types) + 100, 1000]);
        tl  = tiledlayout(numel(tested), numel(types), 'TileSpacing', 'compact');
        for k = 1:numel(tested)
            for t = 1:numel(types)
                ax = nexttile(tl); hold(ax, 'on');
                ty = pert_types{strcmp(struct2cell(pert_type_display), types{t})};
                subs = [bundle_names, {'pooled'}];
                for b = 1:numel(subs)
                    sel = strcmp(Tq.perturbation, types{t}) & strcmp(Tq.subset, subs{b}) & ...
                          strcmp(Tq.metric, tested(k).txt);
                    if b <= nb, col = bcol.(ty)(b, :); else, col = [0 0 0]; end
                    e = Tq.effect(sel); sg = Tq.sig(sel);
                    plot(ax, xpos(1:numel(e)), e, '-', 'Color', col, 'LineWidth', 1.4, ...
                        'DisplayName', subs{b});
                    scatter(ax, xpos(sg), e(sg), 36, col, 'filled', 'HandleVisibility', 'off');
                    scatter(ax, xpos(~sg), e(~sg), 36, col, 'HandleVisibility', 'off');
                end
                yline(ax, 0, ':k', 'HandleVisibility', 'off');
                ylim(ax, [-1.05 1.05]);
                style_level_axis(ax, xpos, lev_lbl, k == numel(tested));
                if k == 1, title(ax, types{t}); end
                if t == 1, ylabel(ax, sprintf('Effect on %s', tested(k).label)); end
                if k == 1 && t == numel(types), legend(ax, 'Location', 'best', 'Box', 'off'); end
            end
        end
        title(tl, sprintf(['%s — effect size vs noise level (%s, whole array)\n' ...
            '+ = %s larger; filled = FDR-significant'], cmp, ...
            ori_titles.(headline_orientation), SYS.(ia).system.label), 'FontWeight', 'bold');
        pt_save_figure(fig, out_dir, sprintf('systems_%s_vs_%s', ia, ib), fig_resolution);
    end
end


% =========================================================================
% F. WITHIN-MODALITY TESTS REPEATED UNDER NOISE
% =========================================================================

fprintf('[F] Within-modality tests at the noise-free and reference levels\n');
tabs = {};
for s = sys_ids
    NS = SYS.(s{1});
    types = pert_types(isfield(NS.pert, pert_types));
    so = struct('types', {types}, 'type_display', pert_type_display, ...
        'bundle_names', {bundle_names}, 'n_perm', stats_n_perm);
    if source_sensor_paired, so.paired_types = {'source', 'sensor'};
    else, so.paired_types = cell(0, 2); end
    c = numel(NS.axcfg_names);
    for L = unique([1, L_ref])
        for o = 1:numel(NS.ori_labels)
            for k = 1:numel(tested)
                vals = struct(); bund = struct(); mags = struct();
                for t = 1:numel(types)
                    v = shiftval(NS, types{t}, tested(k).name, L, c, o);
                    vals.(types{t}) = v(2:end);
                    bund.(types{t}) = NS.pert.(types{t}).bundle_idx;
                    mags.(types{t}) = NS.pert.(types{t}).magnitude;
                end
                T = pt_stats_within(vals, bund, mags, so);
                T.system       = repmat({NS.system.label}, height(T), 1);
                T.noise_factor = repmat(levels(L), height(T), 1);
                T.orientation  = repmat(NS.ori_labels(o), height(T), 1);
                T.metric       = repmat({tested(k).txt}, height(T), 1);
                tabs{end+1} = T; %#ok<SAGROW>
            end
        end
    end
end
Tw = vertcat(tabs{:});
Tw = pt_fdr_by_group(Tw, {'system', 'noise_factor', 'family', 'metric'}, stats_alpha);
Tw = movevars(Tw, {'system', 'noise_factor', 'metric', 'orientation'}, 'Before', 1);
Tw.group = strcat(Tw.system, {' — '}, Tw.family, {' — '}, Tw.metric);
Tw = sortrows(Tw, {'system', 'family', 'metric'});
pt_write_table(Tw, fullfile(out_dir, 'noise_within_stats'), struct( ...
    'cols', {{ 'noise_factor', 'x floor', '$\times$ floor', '%g'; ...
               'subset', 'Subset', 'Subset', '%s'; ...
               'comparison', 'Comparison', 'Comparison', '%s'; ...
               'median_diff', 'median diff', 'Median diff.', '%+.3g'; ...
               'effect', 'effect', 'Effect', '%+.2f'; ...
               'p_fdr', 'p (FDR)', '$p_\mathrm{FDR}$', 'p' }}, ...
    'title', sprintf('STAGE 5 — WITHIN-MODALITY TESTS, NOISE-FREE AND AT %gx FLOOR (whole array)', noise_reference_factor), ...
    'notes', {{'The stage 2/3 tests, repeated on measured (noisy) values. Compare each 0x row with its reference-level row.', ...
               'p_fdr: Benjamini-Hochberg within system x noise level x family x metric.'}}, ...
    'caption', sprintf(['Within-modality perturbation tests on the noise-free field ' ...
               'and at %g$\\times$ each system''s noise floor.'], noise_reference_factor), ...
    'label', 'tab:noise-within', 'group_by', 'group', ...
    'tex_rows', strcmp(Tw.orientation, headline_orientation) & ...
                ~strcmp(Tw.family, 'bundle_pairs')));

writetable(long_rows, fullfile(out_dir, 'stage5_long.csv'));
fprintf('\nStage 5 complete: %s\n', out_dir);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function ty = first_type(NS, order)
    f  = order(isfield(NS.pert, order));
    ty = f{1};
end

function style_level_axis(ax, xpos, lbl, show_x)
    set(ax, 'XTick', xpos, 'XTickLabel', lbl, 'XLim', [xpos(1) - 0.3, xpos(end) + 0.3], ...
        'TickDir', 'out', 'FontSize', 9);
    grid(ax, 'on');
    if show_x, xlabel(ax, 'Noise (\times system floor)'); end
end

function [crit, txt] = critical_level(levels, rate, target)
% Highest noise level at which the detection rate still reaches target,
% log-interpolated between simulated levels. Level 0 is excluded.
    nz = levels > 0;
    lv = levels(nz); rt = rate(nz);
    if all(rt >= target)
        crit = Inf; txt = sprintf('> %g', lv(end));
    elseif rt(1) < target
        crit = 0;   txt = sprintf('< %g', lv(1));
    else
        i = find(rt >= target, 1, 'last');
        if i == numel(lv)
            crit = lv(end);
        else
            f = (rt(i) - target) / (rt(i) - rt(i+1));
            crit = exp(log(lv(i)) + f * (log(lv(i+1)) - log(lv(i))));
        end
        txt = sprintf('%.3g', crit);
    end
end

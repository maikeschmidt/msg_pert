% pt_compare_modalities - Stage 4: MSG versus ESG under identical perturbations
%
% Takes every analysis of stages 2 and 3 and asks whether the two modalities
% respond differently to the same perturbation, on the noise-free forward
% fields and with the same forward model type for both (config_pert:
% modality_compare_method). Five analyses:
%
%   A. Side-by-side description     every metric, type, bundle, axis set
%   B. MSG vs ESG tests             per type x bundle (and pooled), paired by
%                                   realisation where the two modalities share
%                                   the same geometric change (msg_esg_paired)
%   C. Interaction                  is the source-vs-sensor difference itself
%                                   different between MSG and ESG?
%   D. Dose-response slopes         metric change per mm (or per %) in each
%                                   modality, with a bootstrap CI of the
%                                   difference
%   E. Along the cord               where MSG and ESG differ, per-position
%                                   tests FDR-corrected across positions
%
% AXIS MATCHING
%   MSG is triaxial, ESG has tangential and radial electrodes. Comparisons
%   are made on the whole array (all axes stacked; needs no matching, and is
%   the headline) and on matched axes by comparison slot (config_pert:
%   axis_slot), i.e. ESG tangential vs MSG X and ESG radial vs MSG Z.
%
%   The metrics are all ratios (RE, gain) or shape measures (r2, RDM), so
%   they are dimensionless and comparable although MSG is in fT/nAm and ESG
%   in uV/nAm.
%
% USAGE:
%   pt_compare_modalities
%
% INPUT:
%   pert_metrics.mat from pt_compute_metrics, for both modalities
%
% OUTPUT (to <stage_results_dir>/4_msg_vs_esg/):
%   headline_<axes>.png                   metrics x types, MSG vs ESG per bundle
%   along_cord_<type>.png                 cord profiles with significant positions
%   effect_size_map.png                   effect sizes, all types/bundles/axes
%   effect_size_map_orientations.png      the same across dipole orientations
%   gain_vs_rdm.png                       amplitude vs topography, both modalities
%   dose_response.png                     metric vs perturbation size, both
%   modality_descriptive.csv/.txt/.tex    section A
%   modality_stats.csv/.txt/.tex          section B
%   modality_interaction.csv/.txt/.tex    section C
%   modality_slopes.csv/.txt/.tex         section D
%   modality_per_source.csv               section E
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

fprintf('Stage 4 — MSG vs ESG (noise-free)\n\n');

out_dir = fullfile(stage_results_dir, '4_msg_vs_esg');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

A = load(fullfile(mods_cfg.msg.forward_fields_base, 'pert_metrics.mat'), 'PM'); A = A.PM;
B = load(fullfile(mods_cfg.esg.forward_fields_base, 'pert_metrics.mat'), 'PM'); B = B.PM;

meth  = modality_compare_method;
info  = pt_metric_info();
tested = info([info.tested]);
oris  = A.ori_labels;
o_head = find(strcmp(oris, headline_orientation));
nb    = numel(bundle_names);
dist  = A.dist_mm;
if ~isequal(A.dist_mm, B.dist_mm)
    error('MSG and ESG source grids differ; per-position comparison is not possible.');
end

types = pert_types(cellfun(@(ty) isfield(A.pert, ty) && isfield(B.pert, ty) && ...
    isfield(A.pert.(ty).M, meth) && isfield(B.pert.(ty).M, meth), pert_types));
fprintf('  Method: %s   Types in both: %s\n', pt_method_label(meth), strjoin(types, ', '));

% Axis matches: {name, MSG axis-set index, ESG axis-set index}
amatch = {'Array', numel(A.axcfg_names), numel(B.axcfg_names)};
for k = 1:numel(B.axcfg_slot)
    if B.axcfg_slot(k) == 0, continue; end
    ia = find(A.axcfg_slot == B.axcfg_slot(k), 1);
    if isempty(ia), continue; end
    amatch(end+1, :) = {sprintf('%s / %s', B.axcfg_names{k}, A.axcfg_names{ia}), ia, k}; %#ok<SAGROW>
end
fprintf('  Axis sets compared: %s\n\n', strjoin(amatch(:, 1)', ', '));

sopts = struct('bundle_names', {bundle_names}, 'n_perm', stats_n_perm);
colA = modality_colors.msg; colB = modality_colors.esg;
vals = @(PM, ty, metric, c, o) pt_shift_summary(PM.pert.(ty).M.(meth), metric, c, o);


% =========================================================================
% A + B. DESCRIPTION AND TESTS
% =========================================================================

fprintf('[A/B] Description and MSG vs ESG tests\n');
desc_rows = {};    % long: one row per metric
wide_rows = {};    % wide: one row per type x bundle x ori x axes
tabs = {};
for t = 1:numel(types)
    ty = types{t};
    ba = A.pert.(ty).bundle_idx; bb = B.pert.(ty).bundle_idx;
    paired = msg_esg_paired.(ty);
    for a = 1:size(amatch, 1)
        for o = 1:numel(oris)
            VA = cell(1, numel(info)); VB = VA;
            for k = 1:numel(info)
                VA{k} = vals(A, ty, info(k).name, amatch{a, 2}, o);
                VB{k} = vals(B, ty, info(k).name, amatch{a, 3}, o);
                if ~info(k).tested, continue; end
                T = pt_stats_between(VA{k}, VB{k}, ba, bb, paired, sopts, 'MSG vs ESG');
                T.perturbation = repmat({pert_type_display.(ty)}, height(T), 1);
                T.orientation  = repmat(oris(o), height(T), 1);
                T.axis_set     = repmat(amatch(a, 1), height(T), 1);
                T.metric       = repmat({info(k).txt}, height(T), 1);
                tabs{end+1} = T; %#ok<SAGROW>
            end
            for b = [1:nb, 0]
                if b == 0, sa = true(size(ba)); sb = true(size(bb)); sub = 'pooled';
                else, sa = ba == b; sb = bb == b; sub = bundle_names{b}; end
                wr = {pert_type_display.(ty), sub, oris{o}, amatch{a, 1}};
                for k = 1:numel(info)
                    ma = median(VA{k}(sa), 'omitnan'); mb = median(VB{k}(sb), 'omitnan');
                    desc_rows(end+1, :) = {pert_type_display.(ty), sub, oris{o}, ...
                        amatch{a, 1}, info(k).txt, ma, mb, mb / ma}; %#ok<SAGROW>
                    wr = [wr, {ma, mb}]; %#ok<AGROW>
                end
                wide_rows(end+1, :) = wr; %#ok<SAGROW>
            end
        end
    end
end

Tdesc = cell2table(desc_rows, 'VariableNames', {'perturbation', 'bundle', ...
    'orientation', 'axis_set', 'metric', 'msg_median', 'esg_median', 'esg_over_msg'});
writetable(Tdesc, fullfile(out_dir, 'modality_descriptive_long.csv'));

wvars = {'perturbation', 'bundle', 'orientation', 'axis_set'};
for k = 1:numel(info)
    wvars = [wvars, {['msg_' info(k).name], ['esg_' info(k).name]}]; %#ok<AGROW>
end
Tdw = cell2table(wide_rows, 'VariableNames', wvars);
Tdw.group = strcat(Tdw.perturbation, {' — '}, Tdw.axis_set);
pt_write_table(Tdw, fullfile(out_dir, 'modality_descriptive'), struct( ...
    'cols', {{ 'bundle', 'Bundle', 'Bundle', '%s'; ...
               'orientation', 'Ori', 'Ori.', '%s'; ...
               'msg_re', 'RE MSG', 'RE MSG', '%.2f'; ...
               'esg_re', 'RE ESG', 'RE ESG', '%.2f'; ...
               'msg_rsq', 'r2 MSG', '$r^2$ MSG', '%.4f'; ...
               'esg_rsq', 'r2 ESG', '$r^2$ ESG', '%.4f'; ...
               'msg_rdm', 'RDM MSG', 'RDM MSG', '%.4f'; ...
               'esg_rdm', 'RDM ESG', 'RDM ESG', '%.4f'; ...
               'msg_gain', 'gain MSG', 'Gain MSG', '%+.2f'; ...
               'esg_gain', 'gain ESG', 'Gain ESG', '%+.2f' }}, ...
    'title', 'STAGE 4 — MSG AND ESG SIDE BY SIDE (noise-free)', ...
    'notes', {{sprintf('Forward model: %s for both. Medians over realisations of the cord median.', pt_method_label(meth)), ...
               'modality_descriptive_long.csv also carries the ESG/MSG ratio of every metric.'}}, ...
    'caption', sprintf(['Perturbation effects on MSG and ESG side by side (%s, ' ...
               'noise-free). Values are medians over realisations of each ' ...
               'realisation''s median over the cord.'], pt_method_label(meth)), ...
    'label', 'tab:msg-esg-descriptive', 'group_by', 'group', ...
    'tex_rows', strcmp(Tdw.orientation, headline_orientation) & ...
                strcmp(Tdw.axis_set, 'Array')));

Ts = vertcat(tabs{:});
Ts = pt_fdr_by_group(Ts, {'metric'}, stats_alpha);
Ts = movevars(Ts, {'perturbation', 'metric', 'orientation', 'axis_set'}, 'Before', 1);
Ts.group = strcat(Ts.metric, {' — '}, Ts.axis_set);
Ts = sortrows(Ts, {'metric', 'axis_set'});
pt_write_table(Ts, fullfile(out_dir, 'modality_stats'), struct( ...
    'cols', {{ 'perturbation', 'Perturbation', 'Perturbation', '%s'; ...
               'subset', 'Bundle', 'Bundle', '%s'; ...
               'test', 'Test', 'Test', '%s'; ...
               'n_a', 'n', '$n$', '%d'; ...
               'median_a', 'MSG', 'MSG', '%.3g'; ...
               'median_b', 'ESG', 'ESG', '%.3g'; ...
               'median_diff', 'MSG-ESG', 'MSG$-$ESG', '%+.3g'; ...
               'effect', 'effect', 'Effect', '%+.2f'; ...
               'p_fdr', 'p (FDR)', '$p_\mathrm{FDR}$', 'p' }}, ...
    'title', 'STAGE 4 — MSG VS ESG TESTS (noise-free)', ...
    'notes', {{'Positive effect = MSG shows the larger value of the metric.', ...
               'Effect: rank-biserial when paired, Cliff delta when unpaired.', ...
               'p_fdr: Benjamini-Hochberg within each metric across types, bundles, orientations and axis sets.'}}, ...
    'caption', ['MSG versus ESG under identical perturbations (noise-free). ' ...
               'Each observation is one perturbation realisation summarised over ' ...
               'the cord; tests are paired where the two modalities share the ' ...
               'same geometric change. Positive effect sizes mean MSG shows the ' ...
               'larger value. $p_\mathrm{FDR}$ is adjusted within each metric.'], ...
    'label', 'tab:msg-esg-stats', 'group_by', 'group', ...
    'tex_rows', strcmp(Ts.orientation, headline_orientation)));


% =========================================================================
% C. INTERACTION: is (source - sensor) different between MSG and ESG?
% =========================================================================

Tint = table();
if all(ismember({'source', 'sensor'}, types)) && source_sensor_paired
    fprintf('[C] Modality x perturbation-type interaction\n');
    pr = msg_esg_paired.source && msg_esg_paired.sensor;
    tabs = {};
    for a = 1:size(amatch, 1)
        for o = 1:numel(oris)
            for k = 1:numel(tested)
                da = vals(A, 'source', tested(k).name, amatch{a, 2}, o) - ...
                     vals(A, 'sensor', tested(k).name, amatch{a, 2}, o);
                db = vals(B, 'source', tested(k).name, amatch{a, 3}, o) - ...
                     vals(B, 'sensor', tested(k).name, amatch{a, 3}, o);
                T = pt_stats_between(da, db, A.pert.source.bundle_idx, ...
                    B.pert.source.bundle_idx, pr, sopts, ...
                    '(source - sensor): MSG vs ESG');
                T.orientation = repmat(oris(o), height(T), 1);
                T.axis_set    = repmat(amatch(a, 1), height(T), 1);
                T.metric      = repmat({tested(k).txt}, height(T), 1);
                tabs{end+1} = T; %#ok<SAGROW>
            end
        end
    end
    Tint = vertcat(tabs{:});
    Tint = pt_fdr_by_group(Tint, {'metric'}, stats_alpha);
    Tint = movevars(Tint, {'metric', 'orientation', 'axis_set'}, 'Before', 1);
    Tint.group = strcat(Tint.metric, {' — '}, Tint.axis_set);
    pt_write_table(Tint, fullfile(out_dir, 'modality_interaction'), struct( ...
        'cols', {{ 'subset', 'Bundle', 'Bundle', '%s'; ...
                   'orientation', 'Ori', 'Ori.', '%s'; ...
                   'median_a', 'MSG src-sen', 'MSG (src$-$sen)', '%+.3g'; ...
                   'median_b', 'ESG src-sen', 'ESG (src$-$sen)', '%+.3g'; ...
                   'effect', 'effect', 'Effect', '%+.2f'; ...
                   'p_fdr', 'p (FDR)', '$p_\mathrm{FDR}$', 'p' }}, ...
        'title', 'STAGE 4 — INTERACTION: SOURCE-MINUS-SENSOR DIFFERENCE, MSG VS ESG', ...
        'notes', {{'Within each modality, source and sensor shift k share a shift vector, so (source - sensor) is a paired difference.', ...
                   'A significant row means the two perturbation types rank differently in the two modalities.'}}, ...
        'caption', ['Interaction between modality and perturbation type: the ' ...
                   'paired source-minus-sensor difference in each modality, ' ...
                   'compared between MSG and ESG.'], ...
        'label', 'tab:msg-esg-interaction', 'group_by', 'group', ...
        'tex_rows', strcmp(Tint.orientation, headline_orientation)));
end


% =========================================================================
% D. DOSE-RESPONSE SLOPES
% =========================================================================

fprintf('[D] Dose-response slopes\n');
sl_rows = {};
for t = 1:numel(types)
    ty = types{t};
    ma = A.pert.(ty).magnitude; mbm = B.pert.(ty).magnitude;
    for a = 1:size(amatch, 1)
        for k = 1:numel(tested)
            va = vals(A, ty, tested(k).name, amatch{a, 2}, o_head);
            vb = vals(B, ty, tested(k).name, amatch{a, 3}, o_head);
            [sa, sb, ci] = slope_diff(ma, va, mbm, vb, msg_esg_paired.(ty), stats_n_boot);
            sl_rows(end+1, :) = {pert_type_display.(ty), amatch{a, 1}, tested(k).txt, ...
                A.pert.(ty).mag_unit, sa, sb, sa - sb, ci(1), ci(2), ...
                ~isnan(ci(1)) && (ci(1) > 0 || ci(2) < 0)}; %#ok<SAGROW>
        end
    end
end
Tsl = cell2table(sl_rows, 'VariableNames', {'perturbation', 'axis_set', 'metric', ...
    'per_unit', 'slope_msg', 'slope_esg', 'difference', 'difference_ci_lo', ...
    'difference_ci_hi', 'ci_excludes_zero'});
Tsl.group = Tsl.perturbation;
pt_write_table(Tsl, fullfile(out_dir, 'modality_slopes'), struct( ...
    'cols', {{ 'axis_set', 'Axes', 'Axes', '%s'; ...
               'metric', 'Metric', 'Metric', '%s'; ...
               'slope_msg', 'MSG /unit', 'MSG per unit', '%+.4g'; ...
               'slope_esg', 'ESG /unit', 'ESG per unit', '%+.4g'; ...
               'difference', 'MSG-ESG', 'MSG$-$ESG', '%+.4g'; ...
               'difference_ci_lo', 'CI lo', 'CI low', '%+.4g'; ...
               'difference_ci_hi', 'CI hi', 'CI high', '%+.4g' }}, ...
    'title', sprintf('STAGE 4 — DOSE-RESPONSE SLOPES (%s)', ori_titles.(headline_orientation)), ...
    'notes', {{'Least-squares slope of the per-realisation metric against perturbation size (per mm, or per % conductivity change).', ...
               '95% bootstrap CI of the MSG-ESG difference, resampling realisations (paired where the modalities are paired).'}}, ...
    'caption', ['Rate at which each metric changes with perturbation size, per ' ...
               'modality, with a 95\% bootstrap CI of the difference.'], ...
    'label', 'tab:msg-esg-slopes', 'group_by', 'group'));


% =========================================================================
% E. ALONG THE CORD (per position, FDR across positions)
% =========================================================================

fprintf('[E] Along the cord\n');
plot_metrics = {'re', 'absgain', 'rdm', 'rsq'};
ps_rows = {};
for t = 1:numel(types)
    ty = types{t};
    ba = A.pert.(ty).bundle_idx; bb = B.pert.(ty).bundle_idx;
    paired = msg_esg_paired.(ty);
    fig = figure('Color', 'w', 'Position', [40 40 1500 1100]);
    tl  = tiledlayout(numel(plot_metrics), nb + 1, 'TileSpacing', 'compact');
    for r = 1:numel(plot_metrics)
        [~, Xa] = pt_shift_summary(A.pert.(ty).M.(meth), plot_metrics{r}, amatch{1, 2}, o_head);
        [~, Xb] = pt_shift_summary(B.pert.(ty).M.(meth), plot_metrics{r}, amatch{1, 3}, o_head);
        for b = [1:nb, 0]
            if b == 0, sa = true(size(ba)); sb = true(size(bb)); sub = 'pooled';
            else, sa = ba == b; sb = bb == b; sub = bundle_names{b}; end
            p = nan(1, numel(dist));
            for s = 1:numel(dist)
                R = pt_compare_two(Xa(sa, s), Xb(sb, s), paired && isequal(sa, sb), stats_n_perm);
                p(s) = R.p;
            end
            [padj, sig] = st_bh_fdr(p, stats_alpha);
            for s = 1:numel(dist)
                ps_rows(end+1, :) = {pert_type_display.(ty), plot_metrics{r}, sub, ...
                    dist(s), median(Xa(sa, s), 'omitnan'), median(Xb(sb, s), 'omitnan'), ...
                    p(s), padj(s), sig(s)}; %#ok<SAGROW>
            end
            ax = nexttile(tl); hold(ax, 'on');
            h1 = pt_plot_band(ax, dist, Xa(sa, :), colA, '-', 'MSG');
            h2 = pt_plot_band(ax, dist, Xb(sb, :), colB, '-', 'ESG');
            yl = ylim(ax);
            plot(ax, dist(sig), repmat(yl(2) - 0.04*diff(yl), 1, sum(sig)), 's', ...
                'MarkerSize', 3, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none', ...
                'HandleVisibility', 'off');
            grid(ax, 'on'); set(ax, 'TickDir', 'out', 'FontSize', 9);
            if r == 1, title(ax, sprintf('%s bundle', sub)); end
            if b == 1, ylabel(ax, info(strcmp({info.name}, plot_metrics{r})).label); end
            if r == numel(plot_metrics), xlabel(ax, 'Distance along cord (mm)'); end
            if r == 1 && b == 0, legend(ax, [h1 h2], 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, sprintf(['MSG vs ESG — %s perturbation along the cord (%s, %s, whole array)\n' ...
        'median and IQR over realisations; squares = FDR-significant positions'], ...
        lower(pert_type_display.(ty)), pt_method_label(meth), ...
        ori_titles.(headline_orientation)), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['along_cord_' ty], fig_resolution);
end
writetable(cell2table(ps_rows, 'VariableNames', {'perturbation', 'metric', 'bundle', ...
    'distance_mm', 'msg_median', 'esg_median', 'p', 'p_fdr', 'sig'}), ...
    fullfile(out_dir, 'modality_per_source.csv'));


% =========================================================================
% FIGURES: headline, effect-size maps, gain vs RDM, dose response
% =========================================================================

fprintf('[F] Summary figures\n');

% Headline: rows = metrics, cols = types, x = bundle, MSG and ESG side by side
for a = 1:size(amatch, 1)
    fig = figure('Color', 'w', 'Position', [40 40 420*numel(types) + 100, 1050]);
    tl  = tiledlayout(numel(tested), numel(types), 'TileSpacing', 'compact');
    for k = 1:numel(tested)
        for t = 1:numel(types)
            ty = types{t};
            va = vals(A, ty, tested(k).name, amatch{a, 2}, o_head);
            vb = vals(B, ty, tested(k).name, amatch{a, 3}, o_head);
            ax = nexttile(tl); hold(ax, 'on');
            stars = cell(1, nb);
            for b = 1:nb
                xa = va(A.pert.(ty).bundle_idx == b); xa = xa(~isnan(xa));
                xb = vb(B.pert.(ty).bundle_idx == b); xb = xb(~isnan(xb));
                if ~isempty(xa)
                    boxchart(ax, (b - 0.18) * ones(size(xa)), xa, 'BoxWidth', 0.3, ...
                        'BoxFaceColor', colA, 'MarkerStyle', 'none');
                end
                if ~isempty(xb)
                    boxchart(ax, (b + 0.18) * ones(size(xb)), xb, 'BoxWidth', 0.3, ...
                        'BoxFaceColor', colB, 'MarkerStyle', 'none');
                end
                row = Ts(strcmp(Ts.perturbation, pert_type_display.(ty)) & ...
                    strcmp(Ts.subset, bundle_names{b}) & strcmp(Ts.metric, tested(k).txt) & ...
                    strcmp(Ts.orientation, headline_orientation) & ...
                    strcmp(Ts.axis_set, amatch{a, 1}), :);
                if ~isempty(row), stars{b} = pt_sig_label(row.p_fdr(1)); end
            end
            % Stars last, once every box is drawn and the y-limits are final
            yl = ylim(ax);
            ylim(ax, [yl(1), yl(2) + 0.08*diff(yl)]);
            for b = 1:nb
                text(ax, b, yl(2) + 0.04*diff(yl), stars{b}, ...
                    'HorizontalAlignment', 'center', 'FontSize', 10);
            end
            set(ax, 'XTick', 1:nb, 'XTickLabel', bundle_names, 'TickDir', 'out', ...
                'FontSize', 9, 'XLim', [0.4 nb + 0.6]);
            grid(ax, 'on');
            if k == 1, title(ax, pert_type_display.(ty)); end
            if t == 1, ylabel(ax, tested(k).label); end
            if k == 1 && t == numel(types)
                hA = plot(ax, NaN, NaN, 's', 'MarkerFaceColor', colA, 'MarkerEdgeColor', 'none');
                hB = plot(ax, NaN, NaN, 's', 'MarkerFaceColor', colB, 'MarkerEdgeColor', 'none');
                legend(ax, [hA hB], {'MSG', 'ESG'}, 'Location', 'best', 'Box', 'off');
            end
        end
    end
    title(tl, sprintf(['MSG vs ESG under identical perturbations — %s axes (%s, %s)\n' ...
        'one value per realisation; stars = FDR-adjusted MSG vs ESG test per bundle'], ...
        amatch{a, 1}, pt_method_label(meth), ori_titles.(headline_orientation)), ...
        'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['headline_' matlab.lang.makeValidName(amatch{a, 1})], ...
        fig_resolution);
end

% Effect-size maps
effect_map(Ts, types, pert_type_display, bundle_names, tested, amatch(:, 1)', ...
    {headline_orientation}, out_dir, 'effect_size_map', fig_resolution);
effect_map(Ts, types, pert_type_display, bundle_names, tested, {'Array'}, ...
    oris, out_dir, 'effect_size_map_orientations', fig_resolution);

% Gain vs RDM, both modalities
mk = {'o', 's', '^'};
fig = figure('Color', 'w', 'Position', [40 40 420*numel(types), 430]);
tl  = tiledlayout(1, numel(types), 'TileSpacing', 'compact');
for t = 1:numel(types)
    ty = types{t};
    ax = nexttile(tl); hold(ax, 'on');
    for side = 1:2
        if side == 1, PMx = A; col = colA; c = amatch{1, 2}; nm = 'MSG';
        else,         PMx = B; col = colB; c = amatch{1, 3}; nm = 'ESG'; end
        g = vals(PMx, ty, 'gain', c, o_head);
        d = vals(PMx, ty, 'rdm',  c, o_head);
        for b = 1:nb
            sel = PMx.pert.(ty).bundle_idx == b;
            scatter(ax, g(sel), d(sel), 30 + 25*b, col, mk{b}, 'filled', ...
                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.75, ...
                'DisplayName', sprintf('%s, %s', nm, bundle_names{b}));
        end
    end
    xline(ax, 0, ':k', 'HandleVisibility', 'off');
    grid(ax, 'on'); set(ax, 'TickDir', 'out');
    title(ax, pert_type_display.(ty)); xlabel(ax, 'Gain (%)');
    if t == 1, ylabel(ax, 'RDM (topography change)'); end
    if t == numel(types), legend(ax, 'Location', 'bestoutside', 'FontSize', 7, 'Box', 'off'); end
end
title(tl, sprintf('Amplitude vs topography per realisation — MSG and ESG (%s, whole array)', ...
    ori_titles.(headline_orientation)), 'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'gain_vs_rdm', fig_resolution);

% Dose response, both modalities
fig = figure('Color', 'w', 'Position', [40 40 400*numel(types) + 100, 1050]);
tl  = tiledlayout(numel(tested), numel(types), 'TileSpacing', 'compact');
for k = 1:numel(tested)
    for t = 1:numel(types)
        ty = types{t};
        ax = nexttile(tl); hold(ax, 'on');
        for side = 1:2
            if side == 1, PMx = A; col = colA; c = amatch{1, 2}; nm = 'MSG';
            else,         PMx = B; col = colB; c = amatch{1, 3}; nm = 'ESG'; end
            x = PMx.pert.(ty).magnitude;
            y = vals(PMx, ty, tested(k).name, c, o_head);
            scatter(ax, x, y, 30, col, 'filled', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', 0.7, 'DisplayName', nm);
            ok = ~isnan(x) & ~isnan(y);
            if sum(ok) >= 3
                pf = polyfit(x(ok), y(ok), 1);
                xx = linspace(min(x(ok)), max(x(ok)), 20);
                plot(ax, xx, polyval(pf, xx), '-', 'Color', col, 'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
            end
        end
        grid(ax, 'on'); set(ax, 'TickDir', 'out', 'FontSize', 9);
        if k == 1, title(ax, pert_type_display.(ty)); end
        if t == 1, ylabel(ax, tested(k).label); end
        if k == numel(tested)
            xlabel(ax, sprintf('Perturbation size (%s)', A.pert.(ty).mag_unit));
        end
        if k == 1 && t == numel(types), legend(ax, 'Location', 'best', 'Box', 'off'); end
    end
end
title(tl, sprintf('Dose response — MSG and ESG (%s, %s, whole array)', ...
    pt_method_label(meth), ori_titles.(headline_orientation)), 'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'dose_response', fig_resolution);

fprintf('\nStage 4 complete: %s\n', out_dir);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [sa, sb, ci] = slope_diff(xa, ya, xb, yb, paired, n_boot)
% Least-squares slopes of y on x for two groups, and a percentile bootstrap
% CI of their difference, resampling realisations (jointly when paired).
    sa = ls_slope(xa, ya);
    sb = ls_slope(xb, yb);
    ci = [NaN; NaN];
    if isnan(sa) || isnan(sb), return; end
    na = numel(xa); nb_ = numel(xb);
    d  = nan(n_boot, 1);
    for i = 1:n_boot
        ia = randi(na, na, 1);
        if paired && na == nb_, ib = ia; else, ib = randi(nb_, nb_, 1); end
        d(i) = ls_slope(xa(ia), ya(ia)) - ls_slope(xb(ib), yb(ib));
    end
    ci = pt_quantile(d, [2.5 97.5], 1);
end

function s = ls_slope(x, y)
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3 || var(x) == 0, s = NaN; return; end
    xc = x - mean(x);
    s  = sum(xc .* (y - mean(y))) / sum(xc.^2);
end

function effect_map(Ts, types, tdisp, bnames, tested, axsets, oris, out_dir, name, res)
% Effect sizes as a heatmap: rows = type x bundle, columns = metric x axis
% set (or x orientation). Colour = effect, text = FDR stars.
    rlab = {}; rsel = {};
    for t = 1:numel(types)
        for b = [1:numel(bnames), 0]
            if b == 0, sub = 'pooled'; else, sub = bnames{b}; end
            rlab{end+1} = sprintf('%s — %s', tdisp.(types{t}), sub); %#ok<AGROW>
            rsel{end+1} = {tdisp.(types{t}), sub}; %#ok<AGROW>
        end
    end
    clab = {}; csel = {};
    for k = 1:numel(tested)
        for a = 1:numel(axsets)
            for o = 1:numel(oris)
                clab{end+1} = sprintf('%s | %s | %s', tested(k).txt, axsets{a}, oris{o}); %#ok<AGROW>
                csel{end+1} = {tested(k).txt, axsets{a}, oris{o}}; %#ok<AGROW>
            end
        end
    end
    E = nan(numel(rlab), numel(clab)); P = E;
    for i = 1:numel(rlab)
        for j = 1:numel(clab)
            row = Ts(strcmp(Ts.perturbation, rsel{i}{1}) & strcmp(Ts.subset, rsel{i}{2}) & ...
                strcmp(Ts.metric, csel{j}{1}) & strcmp(Ts.axis_set, csel{j}{2}) & ...
                strcmp(Ts.orientation, csel{j}{3}), :);
            if ~isempty(row), E(i, j) = row.effect(1); P(i, j) = row.p_fdr(1); end
        end
    end
    fig = figure('Color', 'w', 'Position', [40 40 max(700, 70*numel(clab)) 40*numel(rlab) + 250]);
    ax = axes(fig);
    imagesc(ax, E, 'AlphaData', ~isnan(E));
    n = 128; t = linspace(0, 1, n)';
    cm = [ [t, t, ones(n,1)]; [ones(n,1), flipud(t), flipud(t)] ];   % blue-white-red
    colormap(ax, cm); caxis(ax, [-1 1]);
    cb = colorbar(ax); cb.Label.String = 'Effect size (+ = MSG larger)';
    for i = 1:numel(rlab)
        for j = 1:numel(clab)
            if ~isnan(P(i, j))
                text(ax, j, i, pt_sig_label(P(i, j)), 'HorizontalAlignment', 'center', ...
                    'FontSize', 8);
            end
        end
    end
    set(ax, 'XTick', 1:numel(clab), 'XTickLabel', clab, 'XTickLabelRotation', 45, ...
        'YTick', 1:numel(rlab), 'YTickLabel', rlab, 'TickLabelInterpreter', 'none', ...
        'FontSize', 8);
    title(ax, 'MSG vs ESG effect sizes (rank-biserial if paired, Cliff \delta otherwise)');
    pt_save_figure(fig, out_dir, name, res);
end

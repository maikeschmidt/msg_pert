% pt_within_modality - Stages 2 and 3: perturbation effects within one modality
%
% For the active modality (stage 2 = MSG, stage 3 = ESG), on the noise-free
% forward fields, answers:
%
%   A. How large is each perturbation?        descriptive table, all metrics
%   B. Where along the cord, and is it gain    decomposition along the cord,
%      or topography?                          one series per bundle
%   C. Does it scale with the size of the      dose-response: metric vs |shift|
%      perturbation?                           (Spearman) and bundle-vs-bundle
%   D. Which perturbation type matters most?   source vs sensor vs conductivity,
%                                              per bundle and pooled
%   E. Amplitude or shape, per realisation?    gain vs RDM scatter
%   F. Does the forward model type change the  e.g. Biot-Savart vs BEM, paired
%      sensitivity?                            by realisation (MSG only)
%   G. Where along the cord do source and      per-position paired test,
%      sensor shifts differ?                   FDR-corrected across positions
%
% STATISTICAL UNIT
%   One perturbation realisation (one random shift), summarised by its median
%   over cord positions: n = 8 per bundle, 24 pooled. Source and sensor shift
%   k share one shift vector, so they are compared paired; conductivity is
%   compared unpaired. p-values are Benjamini-Hochberg corrected within each
%   test family and metric, across methods, orientations, axis sets and
%   subsets (column p_fdr). Tests use RE, r2, RDM and |gain|; signed gain is
%   reported descriptively (see pt_metric_info).
%
% USAGE:
%   pt_modality('set', 'msg');  pt_within_modality     % stage 2
%   pt_modality('set', 'esg');  pt_within_modality     % stage 3
%
% INPUT:
%   <forward_fields_base>/pert_metrics.mat   from pt_compute_metrics
%
% OUTPUT (to <stage_results_dir>/2_within_msg/ or 3_within_esg/):
%   decomposition_<type>_<method>_<axes>.png    RE/gain/RDM/r2 along cord, per bundle
%   types_along_cord_<method>.png               the three types overlaid, per bundle
%   dose_response_<method>.png                  metric vs magnitude
%   types_compared_<method>.png                 per-realisation values by type
%   gain_vs_rdm_<method>.png                    amplitude vs topography
%   methods_compared.png                        forward model types (if > 1)
%   source_minus_sensor_along_cord_<method>.png per-position paired test
%   within_descriptive.csv/.txt/.tex            section A
%   within_stats.csv/.txt/.tex                  sections C and D
%   within_methods.csv/.txt/.tex                section F
%   within_per_source.csv                       section G
%   stage<N>_long.csv                           rows for pt_summary_table
%   The .tex tables hold the headline orientation and the whole array; the
%   .csv files hold every orientation and axis set.
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

if strcmp(active_modality, 'msg'), stage_no = 2; else, stage_no = 3; end
stage_tag = sprintf('%d_within_%s', stage_no, active_modality);
out_dir   = fullfile(stage_results_dir, stage_tag);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

pm_file = fullfile(forward_fields_base, 'pert_metrics.mat');
if ~isfile(pm_file)
    error('%s not found. Run pt_compute_metrics first.', pm_file);
end
load(pm_file, 'PM');

MOD   = PM.display;
fprintf('Stage %d — within %s\n\n', stage_no, MOD);

info    = pt_metric_info();
tested  = info([info.tested]);
types   = pert_types(isfield(PM.pert, pert_types));
oris    = PM.ori_labels;
cfgs    = PM.axcfg_names;
c_arr   = numel(cfgs);                         % whole array is last
c_rad   = PM.radial_cfg;
o_head  = find(strcmp(oris, headline_orientation));
dist    = PM.dist_mm;
nb      = numel(bundle_names);

bcol = struct('source', source_bundle_colors, 'sensor', sensor_bundle_colors, ...
              'cond', cond_bundle_colors);

all_methods = {};
for t = 1:numel(types)
    all_methods = union(all_methods, PM.pert.(types{t}).methods, 'stable');
end

sopts = struct('types', {types}, 'type_display', pert_type_display, ...
    'bundle_names', {bundle_names}, 'n_perm', stats_n_perm);
if source_sensor_paired
    sopts.paired_types = {'source', 'sensor'};
else
    sopts.paired_types = cell(0, 2);
end

method_label = @pt_method_label;


% =========================================================================
% A. DESCRIPTIVE TABLE  (+ long rows)
% =========================================================================

fprintf('[A] Descriptive summary\n');
desc_rows = {};
long_rows = pt_long_row([]);

for t = 1:numel(types)
    ty = types{t}; P = PM.pert.(ty);
    for m = 1:numel(P.methods)
        method = P.methods{m}; S = P.M.(method);
        for c = 1:numel(cfgs)
            for o = 1:numel(oris)
                for b = 1:nb
                    sel = P.bundle_idx == b;
                    x = struct();
                    for k = 1:numel(info)
                        v = pt_shift_summary(S, info(k).name, c, o);
                        x.(info(k).name) = v(sel);
                    end
                    if all(isnan(x.re)), continue; end
                    [~, Xre]  = pt_shift_summary(S, 're', c, o);
                    [~, Xrsq] = pt_shift_summary(S, 'rsq', c, o);
                    Xre = Xre(sel, :); Xrsq = Xrsq(sel, :);
                    D = pt_describe(x.re, stats_n_boot);
                    desc_rows(end+1, :) = {pert_type_display.(ty), method_label(method), ...
                        bundle_names{b}, oris{o}, cfgs{c}, D.n, ...
                        median(P.magnitude(sel), 'omitnan'), P.mag_unit, ...
                        D.median, D.ci_lo, D.ci_hi, max(Xre(:)), ...
                        median(x.rsq, 'omitnan'), min(Xrsq(:)), ...
                        median(x.rdm, 'omitnan'), median(x.gain, 'omitnan'), ...
                        median(x.absgain, 'omitnan'), ...
                        median(mean(Xrsq < 0.99, 2, 'omitnan')) * 100, ...
                        median(mean(Xrsq < 0.95, 2, 'omitnan')) * 100}; %#ok<SAGROW>

                    meta = struct('stage', stage_tag, 'modality', PM.modality, ...
                        'system', 'noise-free', 'method', method, 'factor', ty, ...
                        'level', bundle_names{b}, 'orientation', oris{o}, ...
                        'axis_set', cfgs{c}, 'unit', 'realisation');
                    long_rows = [long_rows; pt_long_row(meta, x, stats_n_boot)]; %#ok<AGROW>
                end
            end
        end
    end
end

Td = cell2table(desc_rows, 'VariableNames', {'perturbation', 'method', 'bundle', ...
    'orientation', 'axis_set', 'n', 'magnitude_median', 'magnitude_unit', ...
    're', 're_ci_lo', 're_ci_hi', 're_worst', 'rsq_median', 'rsq_worst', ...
    'rdm_median', 'gain_median', 'absgain_median', 'pct_cord_rsq_below_099', ...
    'pct_cord_rsq_below_095'});
Td.group = strcat(Td.perturbation, {' — '}, Td.method);

pt_write_table(Td, fullfile(out_dir, 'within_descriptive'), struct( ...
    'cols', {{ 'bundle', 'Bundle', 'Bundle', '%s'; ...
               'orientation', 'Ori', 'Ori.', '%s'; ...
               'axis_set', 'Axes', 'Axes', '%s'; ...
               'magnitude_median', 'Size', 'Size', '%.2f'; ...
               're', 'RE % [95% CI]', 'RE (\%) [95\% CI]', 'ci'; ...
               're_worst', 'RE worst', 'RE worst', '%.2f'; ...
               'rsq_median', 'r2', '$r^2$', '%.4f'; ...
               'rsq_worst', 'r2 worst', '$r^2$ worst', '%.4f'; ...
               'rdm_median', 'RDM', 'RDM', '%.4f'; ...
               'gain_median', 'gain %', 'Gain (\%)', '%+.2f'; ...
               'absgain_median', '|gain| %', '$|$Gain$|$ (\%)', '%.2f'; ...
               'pct_cord_rsq_below_095', '% cord r2<0.95', '\% cord $r^2<0.95$', '%.1f' }}, ...
    'title', sprintf('STAGE %d — %s PERTURBATION EFFECTS (noise-free)', stage_no, MOD), ...
    'notes', {{'One value per realisation = median over cord positions; the table gives the median over realisations.', ...
               'CI: 95% bootstrap over realisations. "worst" = extreme over every realisation and cord position.', ...
               'Size = median |shift| (mm) or mean conductivity change (%).'}}, ...
    'caption', sprintf(['Effect of each perturbation on the noise-free %s forward ' ...
               'field. Each realisation is summarised by its median over cord ' ...
               'positions; values are medians over realisations with 95\\%% bootstrap ' ...
               'CIs. Worst-case values are taken over every realisation and cord ' ...
               'position.'], MOD), ...
    'label', sprintf('tab:within-%s-descriptive', PM.modality), 'group_by', 'group', ...
    'tex_rows', strcmp(Td.orientation, headline_orientation) & strcmp(Td.axis_set, 'Array')));


% =========================================================================
% B. DECOMPOSITION ALONG THE CORD, ONE SERIES PER BUNDLE
% =========================================================================

fprintf('[B] Decomposition along the cord\n');
for t = 1:numel(types)
    ty = types{t}; P = PM.pert.(ty);
    for m = 1:numel(P.methods)
        method = P.methods{m}; S = P.M.(method);
        for c = unique([c_arr, c_rad])
            Ss = struct('label', {}, 're', {}, 'gain', {}, 'rdm', {}, 'rsq', {});
            for b = 1:nb
                sel = P.bundle_idx == b;
                Ss(b).label = sprintf('%s (median of %d)', bundle_names{b}, sum(sel));
                for o = 1:numel(oris)
                    for f = {'re', 'gain', 'rdm', 'rsq'}
                        [~, X] = pt_shift_summary(S, f{1}, c, o);
                        Ss(b).(f{1})(o, :) = median(X(sel, :), 1, 'omitnan');
                    end
                end
            end
            plot_metric_decomposition(Ss, struct('dist', dist, ...
                'orientation_labels', {oris}, 'ori_titles', ori_titles, ...
                'colors', bcol.(ty), ...
                'title', sprintf('%s %s perturbation — %s, %s axes', MOD, ...
                    lower(pert_type_display.(ty)), method_label(method), cfgs{c}), ...
                'save_dir', out_dir, ...
                'save_name', sprintf('decomposition_%s_%s_%s', ty, method, cfgs{c})));
        end
    end
end


% =========================================================================
% B2. THE THREE TYPES OVERLAID ALONG THE CORD (IQR over realisations)
% =========================================================================

plot_metrics = {'re', 'absgain', 'rdm', 'rsq'};
for m = 1:numel(all_methods)
    method = all_methods{m};
    fig = figure('Color', 'w', 'Position', [40 40 1500 1100]);
    tl  = tiledlayout(numel(plot_metrics), nb, 'TileSpacing', 'compact');
    for r = 1:numel(plot_metrics)
        for b = 1:nb
            ax = nexttile(tl); hold(ax, 'on');
            for t = 1:numel(types)
                ty = types{t}; P = PM.pert.(ty);
                if ~isfield(P.M, method), continue; end
                [~, X] = pt_shift_summary(P.M.(method), plot_metrics{r}, c_arr, o_head);
                pt_plot_band(ax, dist, X(P.bundle_idx == b, :), ...
                    pert_type_colors.(ty), '-', pert_type_display.(ty));
            end
            grid(ax, 'on'); set(ax, 'TickDir', 'out', 'FontSize', 9);
            if r == 1, title(ax, sprintf('%s bundle', bundle_names{b})); end
            if b == 1, ylabel(ax, info(strcmp({info.name}, plot_metrics{r})).label); end
            if r == numel(plot_metrics), xlabel(ax, 'Distance along cord (mm)'); end
            if r == 1 && b == nb, legend(ax, 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, sprintf(['%s — perturbation types along the cord (%s, %s, whole array)\n' ...
        'line = median over realisations, band = IQR'], MOD, method_label(method), ...
        ori_titles.(headline_orientation)), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['types_along_cord_' method], fig_resolution);
end


% =========================================================================
% C + D. STATISTICS: magnitude and perturbation type
% =========================================================================

fprintf('[C/D] Magnitude and type statistics\n');
stat_tabs = {};
for m = 1:numel(all_methods)
    method = all_methods{m};
    ty_m = types(cellfun(@(ty) isfield(PM.pert.(ty).M, method), types));
    if isempty(ty_m), continue; end
    so = sopts; so.types = ty_m;
    for c = 1:numel(cfgs)
        for o = 1:numel(oris)
            for k = 1:numel(tested)
                vals = struct(); bund = struct(); mags = struct();
                for t = 1:numel(ty_m)
                    ty = ty_m{t}; P = PM.pert.(ty);
                    vals.(ty) = pt_shift_summary(P.M.(method), tested(k).name, c, o);
                    bund.(ty) = P.bundle_idx;
                    mags.(ty) = P.magnitude;
                end
                T = pt_stats_within(vals, bund, mags, so);
                T.method      = repmat({method_label(method)}, height(T), 1);
                T.orientation = repmat(oris(o), height(T), 1);
                T.axis_set    = repmat(cfgs(c), height(T), 1);
                T.metric      = repmat({tested(k).txt}, height(T), 1);
                stat_tabs{end+1} = T; %#ok<SAGROW>
            end
        end
    end
end
Ts = vertcat(stat_tabs{:});
Ts = pt_fdr_by_group(Ts, {'family', 'metric'}, stats_alpha);
Ts = movevars(Ts, {'method', 'metric', 'orientation', 'axis_set'}, 'Before', 1);
Ts.group = strcat(Ts.family, {' — '}, Ts.metric);
Ts = sortrows(Ts, {'family', 'metric', 'method'});   % stable: keeps subset order

pt_write_table(Ts, fullfile(out_dir, 'within_stats'), struct( ...
    'cols', {{ 'method', 'Method', 'Method', '%s'; ...
               'subset', 'Subset', 'Subset', '%s'; ...
               'comparison', 'Comparison', 'Comparison', '%s'; ...
               'n_a', 'n', '$n$', '%d'; ...
               'median_diff', 'median diff', 'Median diff.', '%+.3g'; ...
               'effect', 'effect', 'Effect', '%+.2f'; ...
               'effect_name', 'effect type', 'Effect type', '%s'; ...
               'p', 'p', '$p$', 'p'; ...
               'p_fdr', 'p (FDR)', '$p_\mathrm{FDR}$', 'p' }}, ...
    'title', sprintf('STAGE %d — %s: MAGNITUDE AND TYPE STATISTICS (noise-free)', stage_no, MOD), ...
    'notes', {{'dose_response: Spearman rho between metric and perturbation size.', ...
               'bundle_pairs: larger bundle vs smaller (unpaired); positive effect = larger bundle larger metric.', ...
               'type_pairs: A vs B; positive effect = A larger. Source vs sensor paired by shift vector.', ...
               'p_fdr: Benjamini-Hochberg within family x metric, across methods, orientations, axis sets, subsets.'}}, ...
    'caption', sprintf(['Within-%s statistics. Each observation is one perturbation ' ...
               'realisation summarised over the cord. Dose response: Spearman ' ...
               'correlation with perturbation size (permutation $p$). Bundle and ' ...
               'type comparisons: permutation tests, paired where realisations ' ...
               'share a shift vector, with rank-biserial or Cliff''s $\\delta$ effect ' ...
               'sizes. $p_\\mathrm{FDR}$ is Benjamini--Hochberg adjusted within each ' ...
               'test family and metric.'], MOD), ...
    'label', sprintf('tab:within-%s-stats', PM.modality), 'group_by', 'group', ...
    'tex_rows', strcmp(Ts.orientation, headline_orientation) & strcmp(Ts.axis_set, 'Array')));


% =========================================================================
% C (figure). DOSE RESPONSE
% =========================================================================

fprintf('[C] Dose-response figures\n');
for m = 1:numel(all_methods)
    method = all_methods{m};
    ty_m = types(cellfun(@(ty) isfield(PM.pert.(ty).M, method), types));
    fig = figure('Color', 'w', 'Position', [40 40 380*numel(ty_m) + 100, 1050]);
    tl  = tiledlayout(numel(tested), numel(ty_m), 'TileSpacing', 'compact');
    for k = 1:numel(tested)
        for t = 1:numel(ty_m)
            ty = ty_m{t}; P = PM.pert.(ty);
            v  = pt_shift_summary(P.M.(method), tested(k).name, c_arr, o_head);
            ax = nexttile(tl); hold(ax, 'on');
            for b = 1:nb
                sel = P.bundle_idx == b;
                scatter(ax, P.magnitude(sel), v(sel), 40, bcol.(ty)(b, :), 'filled', ...
                    'MarkerEdgeColor', 'k', 'DisplayName', bundle_names{b});
            end
            row = Ts(strcmp(Ts.family, 'dose_response') & strcmp(Ts.method, method_label(method)) ...
                & strcmp(Ts.subset, pert_type_display.(ty)) & strcmp(Ts.metric, tested(k).txt) ...
                & strcmp(Ts.orientation, headline_orientation) & strcmp(Ts.axis_set, 'Array'), :);
            if ~isempty(row)
                title(ax, sprintf('\\rho = %+.2f  %s', row.effect(1), pt_sig_label(row.p_fdr(1))), ...
                    'FontSize', 9, 'FontWeight', 'normal');
            end
            grid(ax, 'on'); set(ax, 'TickDir', 'out', 'FontSize', 9);
            if k == 1
                subtitle(ax, pert_type_display.(ty), 'FontWeight', 'bold');
            end
            if k == numel(tested), xlabel(ax, sprintf('Perturbation size (%s)', P.mag_unit)); end
            if t == 1, ylabel(ax, tested(k).label); end
            if k == 1 && t == numel(ty_m), legend(ax, 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, sprintf('%s — dose response (%s, %s, whole array; one point per realisation)', ...
        MOD, method_label(method), ori_titles.(headline_orientation)), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['dose_response_' method], fig_resolution);
end


% =========================================================================
% D (figure). TYPES COMPARED, per bundle and pooled
% =========================================================================

fprintf('[D] Type comparison figures\n');
for m = 1:numel(all_methods)
    method = all_methods{m};
    ty_m = types(cellfun(@(ty) isfield(PM.pert.(ty).M, method), types));
    if numel(ty_m) < 2, continue; end
    fig = figure('Color', 'w', 'Position', [40 40 1500 1050]);
    tl  = tiledlayout(numel(tested), nb + 1, 'TileSpacing', 'compact');
    for k = 1:numel(tested)
        for b = [1:nb, 0]
            ax = nexttile(tl); hold(ax, 'on');
            for t = 1:numel(ty_m)
                ty = ty_m{t}; P = PM.pert.(ty);
                v  = pt_shift_summary(P.M.(method), tested(k).name, c_arr, o_head);
                if b > 0, v = v(P.bundle_idx == b); end
                v  = v(~isnan(v));
                if isempty(v), continue; end
                boxchart(ax, t * ones(size(v)), v, 'BoxFaceColor', pert_type_colors.(ty), ...
                    'MarkerStyle', 'none', 'BoxWidth', 0.5);
                scatter(ax, t + 0.12*randn(size(v)), v, 14, pert_type_colors.(ty), ...
                    'filled', 'MarkerFaceAlpha', 0.6);
            end
            if b == 0, sub = 'pooled'; else, sub = bundle_names{b}; end
            rows = Ts(strcmp(Ts.family, 'type_pairs') & strcmp(Ts.method, method_label(method)) ...
                & strcmp(Ts.subset, sub) & strcmp(Ts.metric, tested(k).txt) ...
                & strcmp(Ts.orientation, headline_orientation) & strcmp(Ts.axis_set, 'Array'), :);
            lbl = cell(height(rows), 1);
            for i = 1:height(rows)
                sides  = strsplit(rows.comparison{i}, ' vs ');
                lbl{i} = sprintf('%s-%s %s', sides{1}(1:min(4, end)), ...
                    sides{2}(1:min(4, end)), pt_sig_label(rows.p_fdr(i)));
            end
            title(ax, strjoin(lbl, '   '), 'FontSize', 8, 'FontWeight', 'normal');
            set(ax, 'XTick', 1:numel(ty_m), 'XTickLabel', ...
                cellfun(@(ty) pert_type_display.(ty), ty_m, 'UniformOutput', false), ...
                'TickDir', 'out', 'FontSize', 9);
            grid(ax, 'on');
            if k == 1, subtitle(ax, sub, 'FontWeight', 'bold'); end
            if b == 1, ylabel(ax, tested(k).label); end
        end
    end
    title(tl, sprintf(['%s — perturbation types compared (%s, %s, whole array)\n' ...
        'one point per realisation; stars = FDR-adjusted p'], MOD, method_label(method), ...
        ori_titles.(headline_orientation)), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['types_compared_' method], fig_resolution);
end


% =========================================================================
% E. GAIN VS RDM, per realisation
% =========================================================================

fprintf('[E] Gain vs RDM\n');
mk = {'o', 's', '^'};
for m = 1:numel(all_methods)
    method = all_methods{m};
    fig = figure('Color', 'w', 'Position', [40 40 380*numel(oris), 420]);
    tl  = tiledlayout(1, numel(oris), 'TileSpacing', 'compact');
    for o = 1:numel(oris)
        ax = nexttile(tl); hold(ax, 'on');
        for t = 1:numel(types)
            ty = types{t}; P = PM.pert.(ty);
            if ~isfield(P.M, method), continue; end
            g = pt_shift_summary(P.M.(method), 'gain', c_arr, o);
            d = pt_shift_summary(P.M.(method), 'rdm',  c_arr, o);
            for b = 1:nb
                sel = P.bundle_idx == b;
                scatter(ax, g(sel), d(sel), 30 + 25*b, pert_type_colors.(ty), mk{b}, ...
                    'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.75, ...
                    'DisplayName', sprintf('%s, %s', pert_type_display.(ty), bundle_names{b}));
            end
        end
        xline(ax, 0, ':k', 'HandleVisibility', 'off');
        grid(ax, 'on'); set(ax, 'TickDir', 'out');
        title(ax, ori_titles.(oris{o}));
        xlabel(ax, 'Gain (%)');
        if o == 1, ylabel(ax, 'RDM (topography change)'); end
        if o == numel(oris), legend(ax, 'Location', 'bestoutside', 'FontSize', 7, 'Box', 'off'); end
    end
    title(tl, sprintf(['%s — amplitude vs topography per realisation (%s, whole array)\n' ...
        'near the horizontal axis = pure rescaling'], MOD, method_label(method)), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['gain_vs_rdm_' method], fig_resolution);
end


% =========================================================================
% F. FORWARD MODEL TYPES COMPARED (paired by realisation)
% =========================================================================

Tmeth = table();
if numel(all_methods) >= 2
    fprintf('[F] Forward model types compared\n');
    mpairs = nchoosek(1:numel(all_methods), 2);
    tabs = {};
    for t = 1:numel(types)
        ty = types{t}; P = PM.pert.(ty);
        for q = 1:size(mpairs, 1)
            ma = all_methods{mpairs(q, 1)}; mb = all_methods{mpairs(q, 2)};
            if ~isfield(P.M, ma) || ~isfield(P.M, mb), continue; end
            for c = 1:numel(cfgs)
                for o = 1:numel(oris)
                    for k = 1:numel(tested)
                        va = pt_shift_summary(P.M.(ma), tested(k).name, c, o);
                        vb = pt_shift_summary(P.M.(mb), tested(k).name, c, o);
                        T = pt_stats_between(va, vb, P.bundle_idx, P.bundle_idx, true, ...
                            sopts, sprintf('%s vs %s', method_label(ma), method_label(mb)));
                        T.perturbation = repmat({pert_type_display.(ty)}, height(T), 1);
                        T.orientation  = repmat(oris(o), height(T), 1);
                        T.axis_set     = repmat(cfgs(c), height(T), 1);
                        T.metric       = repmat({tested(k).txt}, height(T), 1);
                        tabs{end+1} = T; %#ok<SAGROW>
                    end
                end
            end
        end
    end
    if ~isempty(tabs)
        Tmeth = vertcat(tabs{:});
        Tmeth = pt_fdr_by_group(Tmeth, {'metric'}, stats_alpha);
        Tmeth = movevars(Tmeth, {'perturbation', 'metric', 'orientation', 'axis_set'}, 'Before', 1);
        Tmeth.group = strcat(Tmeth.comparison, {' — '}, Tmeth.metric);
        pt_write_table(Tmeth, fullfile(out_dir, 'within_methods'), struct( ...
            'cols', {{ 'perturbation', 'Perturbation', 'Perturbation', '%s'; ...
                       'subset', 'Bundle', 'Bundle', '%s'; ...
                       'n_a', 'n', '$n$', '%d'; ...
                       'median_a', 'median A', 'Median A', '%.3g'; ...
                       'median_b', 'median B', 'Median B', '%.3g'; ...
                       'median_diff', 'median A-B', 'Median A$-$B', '%+.3g'; ...
                       'effect', 'rank-biserial', 'Rank-biserial', '%+.2f'; ...
                       'p_fdr', 'p (FDR)', '$p_\mathrm{FDR}$', 'p' }}, ...
            'title', sprintf('STAGE %d — %s: FORWARD MODEL TYPES COMPARED (noise-free)', stage_no, MOD), ...
            'notes', {{'Paired by realisation: the same perturbed geometry solved with both models.', ...
                       'Positive effect = model A shows the larger value of the metric.'}}, ...
            'caption', sprintf(['Sensitivity of each forward model type to the same ' ...
                       'perturbations (%s). Sign-flip permutation tests paired by ' ...
                       'realisation; $p_\\mathrm{FDR}$ adjusted within each metric.'], MOD), ...
            'label', sprintf('tab:within-%s-methods', PM.modality), 'group_by', 'group', ...
            'tex_rows', strcmp(Tmeth.orientation, headline_orientation) & ...
                        strcmp(Tmeth.axis_set, 'Array')));

        % Figure: paired lines per realisation, one column per type
        ma = all_methods{1}; mb = all_methods{2};
        ty_m = types(cellfun(@(ty) isfield(PM.pert.(ty).M, ma) && ...
                               isfield(PM.pert.(ty).M, mb), types));
        fig = figure('Color', 'w', 'Position', [40 40 360*numel(ty_m) + 100, 1050]);
        tl  = tiledlayout(numel(tested), numel(ty_m), 'TileSpacing', 'compact');
        for k = 1:numel(tested)
            for t = 1:numel(ty_m)
                ty = ty_m{t}; P = PM.pert.(ty);
                va = pt_shift_summary(P.M.(ma), tested(k).name, c_arr, o_head);
                vb = pt_shift_summary(P.M.(mb), tested(k).name, c_arr, o_head);
                ax = nexttile(tl); hold(ax, 'on');
                for i = 1:numel(va)
                    plot(ax, [1 2], [va(i) vb(i)], '-o', 'Color', bcol.(ty)(P.bundle_idx(i), :), ...
                        'MarkerFaceColor', bcol.(ty)(P.bundle_idx(i), :), 'MarkerSize', 4);
                end
                set(ax, 'XTick', [1 2], 'XTickLabel', {method_label(ma), method_label(mb)}, ...
                    'XLim', [0.6 2.4], 'TickDir', 'out', 'FontSize', 9);
                grid(ax, 'on');
                row = Tmeth(strcmp(Tmeth.perturbation, pert_type_display.(ty)) & ...
                    strcmp(Tmeth.subset, 'pooled') & strcmp(Tmeth.metric, tested(k).txt) & ...
                    strcmp(Tmeth.orientation, headline_orientation) & strcmp(Tmeth.axis_set, 'Array'), :);
                if ~isempty(row)
                    title(ax, sprintf('r_{rb} = %+.2f  %s', row.effect(1), pt_sig_label(row.p_fdr(1))), ...
                        'FontSize', 9, 'FontWeight', 'normal');
                end
                if k == 1, subtitle(ax, pert_type_display.(ty), 'FontWeight', 'bold'); end
                if t == 1, ylabel(ax, tested(k).label); end
            end
        end
        title(tl, sprintf('%s — forward model types under identical perturbations (%s, whole array)', ...
            MOD, ori_titles.(headline_orientation)), 'FontWeight', 'bold');
        pt_save_figure(fig, out_dir, 'methods_compared', fig_resolution);
    end
end


% =========================================================================
% G. SOURCE VS SENSOR ALONG THE CORD (per-position paired test)
% =========================================================================

ps_rows = {};
if all(isfield(PM.pert, {'source', 'sensor'})) && source_sensor_paired
    fprintf('[G] Source vs sensor along the cord\n');
    Psrc = PM.pert.source; Psen = PM.pert.sensor;
    meth_g = intersect(Psrc.methods, Psen.methods, 'stable');
    for m = 1:numel(meth_g)
        method = meth_g{m};
        fig = figure('Color', 'w', 'Position', [40 40 1500 1050]);
        tl  = tiledlayout(numel(tested), nb + 1, 'TileSpacing', 'compact');
        for k = 1:numel(tested)
            [~, Xa] = pt_shift_summary(Psrc.M.(method), tested(k).name, c_arr, o_head);
            [~, Xb] = pt_shift_summary(Psen.M.(method), tested(k).name, c_arr, o_head);
            Dd = Xa - Xb;
            for b = [1:nb, 0]
                if b == 0, sel = true(size(Dd, 1), 1); sub = 'pooled';
                else,      sel = Psrc.bundle_idx == b; sub = bundle_names{b}; end
                p = nan(1, numel(dist));
                for s = 1:numel(dist)
                    p(s) = st_signflip_test(Dd(sel, s), stats_n_perm, 'both');
                end
                [padj, sig] = st_bh_fdr(p, stats_alpha);
                for s = 1:numel(dist)
                    ps_rows(end+1, :) = {method_label(method), tested(k).txt, sub, ...
                        dist(s), median(Dd(sel, s), 'omitnan'), p(s), padj(s), sig(s)}; %#ok<SAGROW>
                end
                ax = nexttile(tl); hold(ax, 'on');
                pt_plot_band(ax, dist, Dd(sel, :), [0.3 0.3 0.3]);
                yline(ax, 0, ':k');
                yl = ylim(ax);
                plot(ax, dist(sig), repmat(yl(2) - 0.04*diff(yl), 1, sum(sig)), 's', ...
                    'MarkerSize', 3, 'MarkerFaceColor', [0.8 0.1 0.1], 'MarkerEdgeColor', 'none');
                grid(ax, 'on'); set(ax, 'TickDir', 'out', 'FontSize', 9);
                if k == 1, title(ax, sprintf('%s (n = %d)', sub, sum(sel))); end
                if b == 1, ylabel(ax, ['\Delta ' tested(k).label]); end
                if k == numel(tested), xlabel(ax, 'Distance along cord (mm)'); end
            end
        end
        title(tl, sprintf(['%s — source minus sensor shift, along the cord (%s, %s, whole array)\n' ...
            'line = median paired difference, band = IQR, red = FDR-significant position'], ...
            MOD, method_label(method), ori_titles.(headline_orientation)), 'FontWeight', 'bold');
        pt_save_figure(fig, out_dir, ['source_minus_sensor_along_cord_' method], fig_resolution);
    end
end
if ~isempty(ps_rows)
    writetable(cell2table(ps_rows, 'VariableNames', {'method', 'metric', 'bundle', ...
        'distance_mm', 'median_source_minus_sensor', 'p', 'p_fdr', 'sig'}), ...
        fullfile(out_dir, 'within_per_source.csv'));
end

writetable(long_rows, fullfile(out_dir, sprintf('stage%d_long.csv', stage_no)));
fprintf('\nStage %d complete: %s\n', stage_no, out_dir);

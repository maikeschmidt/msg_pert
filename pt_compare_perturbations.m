% pt_compare_perturbations - Cross-perturbation and cross-modality r² comparison
%
% Compares r² distributions across:
%   1. Perturbation types within the same modality — source vs sensor shifts
%      at matched bundle magnitudes (same mm range).
%   2. MSG vs ESG within the same perturbation type and forward model method
%      (source, sensor, or cond).
%
% Two statistical approaches are computed for BOTH comparisons above:
%
%   A) Wilcoxon (signed-rank for source-vs-sensor, rank-sum for MSG-vs-ESG),
%      Bonferroni-corrected within each comparison group. Kept for
%      reference/backward compatibility → comparison_stats.tsv. No longer
%      drives any figure.
%
%   B) Permutation test (label-shuffling), computed both as a per-combo
%      aggregate AND at every individual cord position:
%        - Aggregate: one p-value per method/orientation/axis/bundle combo,
%          no correction needed (a single test) → permutation_stats.tsv.
%          This is what now drives the heatmaps, dumbbell plots, and the
%          significance-rate summary, for BOTH comparisons.
%        - Cord-position: one p-value per cord position within each combo,
%          CLUSTER-CORRECTED across the ~100+ positions tested (see
%          PERMUTATION-BASED ANALYSIS below for the method) → 
%          cord_position_stats.tsv. Drives the cord-position figures.
%
% Bonferroni correction (approach A) is applied within each comparison group
% (not globally). Cluster correction (approach B, cord-position only) is
% applied within each method/orientation/axis/bundle combo, across cord
% position — i.e. it answers "given I tested ~100+ positions along the cord
% for this one combo, how many false positives would that produce by chance,
% and does my finding survive that?"
%
% OUTPUTS (saved to <msg_results_path>/perturbation_analysis/comparison/):
%   comparison_stats.tsv        — Wilcoxon/Bonferroni results, reference only
%   permutation_stats.tsv       — permutation aggregate results (both comparisons)
%   cord_position_stats.tsv     — permutation per-cord-position results,
%                                 uncorrected AND cluster-corrected (both comparisons)
%   source_vs_sensor_<modality>.png/fig — grouped box plots, descriptive only
%   msg_vs_esg_<pert_type>.png/fig      — grouped box plots, descriptive only
%   dumbbell_source_vs_sensor_<modality>_<bundle>.png/fig
%   dumbbell_msg_vs_esg_<pert_type>_<bundle>.png/fig
%   source_vs_sensor_heatmap_<modality>.png/fig
%   msg_vs_esg_heatmap_<pert_type>.png/fig
%   cord_position_sourcesensor_<modality>_<method>_bundle<N>_<name>.png/fig
%   cord_position_msgesg_<pert_type>_<method>_bundle<N>_<name>.png/fig
%   significance_summary.png/fig
%
% DEPENDENCIES:
%   pert_source_rsq.mat, pert_sensor_rsq.mat, pert_cond_rsq.mat
%   (from both MSG and ESG analysis directories when run_msg_vs_esg = true)
%   Statistics and Machine Learning Toolbox (signrank, ranksum, iqr, prctile)
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

fprintf('pt_compare_perturbations\n\n');

%% =========================================================================
%% USER CONFIGURATION
%% =========================================================================

% Path containing pert_source_rsq.mat / pert_sensor_rsq.mat / pert_cond_rsq.mat
% (forward_fields_base from config_pert for the MSG dataset)
msg_results_path = 'D:\Simulations\Pertubations\fields\mag';   % SET THIS

% Path to ESG analysis results (leave '' to skip MSG vs ESG comparisons)
esg_results_path = 'D:\Simulations\Pertubations\fields\elec';   % SET THIS

% Which comparisons to run
run_source_vs_sensor = true;   % source vs sensor (matched bundle, mm-based)
run_msg_vs_esg       = true;   % MSG vs ESG per perturbation type

% Significance level for the Wilcoxon/Bonferroni comparisons (reference only)
alpha = 0.05;

% Significance level for the permutation-based tests (aggregate AND
% cord-position, before and after cluster correction)
alpha_perm = 0.05;

% Number of label-shuffles used to build the permutation null distribution
% (and, for cord-position tests, the cluster-mass null distribution too).
% 100 is the bare minimum for a stable p-value; 2000 gives much finer
% resolution since nothing here is Bonferroni-inflating the threshold.
n_perm = 2000;

% Number of bootstrap resamples for the uncertainty ribbon around the
% observed cord-position difference (resampling within each group, with
% replacement).
n_boot = 1000;

% Extra comparison figures: dumbbell plots, effect-size heatmaps,
% cord-position significance figures, significance-rate summary.
make_stat_figures = true;

% Set false to keep the box plots purely descriptive (box + whiskers only,
% no significance brackets).
show_significance_on_boxplots = false;

% Output directory
save_dir = fullfile(msg_results_path, 'perturbation_analysis', 'comparison');

% Bundle labels (must match the three-bundle structure in config_pert)
bundle_display = {'Small', 'Medium', 'Large'};

% Cord position (mm) for each of the n_positions source points, used as the
% x-axis for the cord-position significance figures. UPDATE THIS if it
% doesn't match your actual source grid, or point it at the real field in
% your .mat files if one exists (e.g. msg.source.source_positions_mm).
cord_position_mm = 10:5:560;   % SET THIS to match your source grid

% Figure colours
col_source = [0.90 0.55 0.10];   % orange — source-position perturbation
col_sensor = [0.25 0.50 0.80];   % blue   — sensor-array perturbation
col_cond   = [0.35 0.70 0.45];
col_MSG    = [0.25 0.50 0.80];   % used for MSG-vs-ESG box plots + dumbbells
col_ESG    = [0.80 0.50 0.20];
% MSG-vs-ESG heatmaps and cord-position figures use a SEPARATE purple/teal
% pair rather than col_MSG/col_ESG, because col_MSG happens to equal
% col_sensor (both blue) — using the same blue for "sensor" in one figure
% and "MSG" in another would be confusing when the two figure families sit
% side by side. Box plots and dumbbells keep col_MSG/col_ESG unchanged.
col_MSG_spatial = [0.31 0.29 0.65];   % purple
col_ESG_spatial = [0.06 0.44 0.35];   % teal

% Publication figure defaults
pub_lw = 1.2;
pub_ms = 6;

%% =========================================================================
%% LOAD RESULTS
%% =========================================================================

fprintf('Loading MSG results from:\n  %s\n', msg_results_path);
msg = load_pert_results(msg_results_path);
fprintf('  Source: %s\n', yesno(msg.have_source));
fprintf('  Sensor: %s\n', yesno(msg.have_sensor));
fprintf('  Cond:   %s\n', yesno(msg.have_cond));

have_esg = ~isempty(esg_results_path) && isfolder(esg_results_path);
if have_esg && run_msg_vs_esg
    fprintf('\nLoading ESG results from:\n  %s\n', esg_results_path);
    esg = load_pert_results(esg_results_path);
    fprintf('  Source: %s\n', yesno(esg.have_source));
    fprintf('  Sensor: %s\n', yesno(esg.have_sensor));
    fprintf('  Cond:   %s\n', yesno(esg.have_cond));
else
    esg = [];
    if run_msg_vs_esg
        fprintf('\nNo ESG path set — skipping MSG vs ESG comparisons.\n');
    end
    have_esg = false;
end

if ~exist(save_dir, 'dir'); mkdir(save_dir); end

%% =========================================================================
%% COLLECT RESULTS TABLE AND RAW PLOT DATA (Wilcoxon / Bonferroni — reference)
%% =========================================================================

all_rows  = {};    % cell of row tables, concatenated at end
plot_data = struct();   % stores raw per-position means for box-plot figures

%% -------------------------------------------------------------------------
%% COMPARISON 1 — Source vs Sensor within each modality (Wilcoxon, Bonferroni)
%% -------------------------------------------------------------------------

if run_source_vs_sensor
    fprintf('\n--- Source vs Sensor comparison (Wilcoxon, Bonferroni) ---\n');

    for mod_idx = 1:2
        if mod_idx == 1
            dat = msg;  mod_label = 'MSG';
        else
            if ~have_esg; continue; end
            dat = esg;  mod_label = 'ESG';
        end

        if ~dat.have_source || ~dat.have_sensor
            fprintf('  %s: source or sensor r² missing — skipping.\n', mod_label);
            continue
        end

        src  = dat.source;
        sens = dat.sensor;
        src_bund  = src.valid_source_bundle_idx;
        sens_bund = sens.valid_bundle_idx;

        methods_both = intersect(src.loaded_methods, sens.loaded_methods);
        if isempty(methods_both)
            fprintf('  %s: no shared methods — skipping.\n', mod_label);
            continue
        end

        ori_labels = fieldnames(src.rsq_by_method.(methods_both{1}));
        n_ori     = numel(ori_labels);
        n_bundles = numel(bundle_display);
        n_axes    = min(src.n_axes, sens.n_axes);

        n_tests = numel(methods_both) * n_bundles * n_ori * n_axes;
        fprintf('  %s: %d methods × %d bundles × %d ori × %d axes = %d tests (Bonferroni)\n', ...
            mod_label, numel(methods_both), n_bundles, n_ori, n_axes, n_tests);

        pd_key = sprintf('sv_s_%s', mod_label);
        plot_data.(pd_key) = struct();

        for m = 1:numel(methods_both)
            method = methods_both{m};
            for b = 1:n_bundles
                sm = src_bund  == b;
                ss = sens_bund == b;
                for o = 1:n_ori
                    ori = ori_labels{o};
                    for ax = 1:n_axes
                        src_rsq_block  = src.rsq_by_method.(method).(ori)(sm, :, ax);
                        sens_rsq_block = sens.rsq_by_method.(method).(ori)(ss, :, ax);

                        src_pos_mean  = mean(src_rsq_block,  1, 'omitnan');
                        sens_pos_mean = mean(sens_rsq_block, 1, 'omitnan');

                        ok = ~isnan(src_pos_mean) & ~isnan(sens_pos_mean);
                        n_ok = sum(ok);
                        if n_ok < 3; continue; end

                        [p_raw, ~] = signrank(src_pos_mean(ok), sens_pos_mean(ok));
                        p_corr     = min(1, p_raw * n_tests);

                        row = build_row('source_vs_sensor', mod_label, 'source', 'sensor', ...
                            method, b, bundle_display{b}, ori, ax, n_ok, ...
                            src_pos_mean(ok), sens_pos_mean(ok), p_raw, p_corr, alpha);
                        all_rows{end+1} = row; %#ok<AGROW>

                        fk = sprintf('%s_%s_%s_ax%d', method, ori, mod_label, ax);
                        if ~isfield(plot_data.(pd_key), fk)
                            plot_data.(pd_key).(fk) = struct('y', [], 'g_bundle', {{}}, 'g_cond', {{}});
                        end
                        pd = plot_data.(pd_key).(fk);
                        for b2 = 1:n_bundles
                            sm2  = src_bund  == b2;
                            ss2  = sens_bund == b2;
                            sv2  = mean(src.rsq_by_method.(method).(ori)(sm2, :, ax), 1, 'omitnan');
                            se2  = mean(sens.rsq_by_method.(method).(ori)(ss2, :, ax), 1, 'omitnan');
                            ok2  = ~isnan(sv2) & ~isnan(se2);
                            pd.y        = [pd.y; sv2(ok2)'; se2(ok2)'];
                            pd.g_bundle = [pd.g_bundle; repmat({bundle_display{b2}}, sum(ok2)*2, 1)];
                            pd.g_cond   = [pd.g_cond;
                                           repmat({'Source'}, sum(ok2), 1);
                                           repmat({'Sensor'}, sum(ok2), 1)];
                        end
                        plot_data.(pd_key).(fk) = pd;
                    end
                end
            end
        end
        fprintf('  %s: done\n', mod_label);
    end
end

%% -------------------------------------------------------------------------
%% COMPARISON 2 — MSG vs ESG per perturbation type (Wilcoxon, Bonferroni)
%% -------------------------------------------------------------------------

if run_msg_vs_esg && have_esg
    fprintf('\n--- MSG vs ESG comparison (Wilcoxon, Bonferroni) ---\n');

    pert_types = {'source', 'sensor', 'cond'};

    for pt = 1:numel(pert_types)
        pert = pert_types{pt};

        switch pert
            case 'source'
                have_both = msg.have_source && esg.have_source;
                if ~have_both; fprintf('  source: missing — skipping.\n'); continue; end
                dat_msg = msg.source; dat_esg = esg.source;
                bund_msg = dat_msg.valid_source_bundle_idx;
                bund_esg = dat_esg.valid_source_bundle_idx;
            case 'sensor'
                have_both = msg.have_sensor && esg.have_sensor;
                if ~have_both; fprintf('  sensor: missing — skipping.\n'); continue; end
                dat_msg = msg.sensor; dat_esg = esg.sensor;
                bund_msg = dat_msg.valid_bundle_idx;
                bund_esg = dat_esg.valid_bundle_idx;
            case 'cond'
                have_both = msg.have_cond && esg.have_cond;
                if ~have_both; fprintf('  cond: missing — skipping.\n'); continue; end
                dat_msg = msg.cond; dat_esg = esg.cond;
                bund_msg = dat_msg.valid_cond_bundle_idx;
                bund_esg = dat_esg.valid_cond_bundle_idx;
        end

        methods_both = intersect(dat_msg.loaded_methods, dat_esg.loaded_methods);
        if isempty(methods_both)
            fprintf('  %s: no shared methods — skipping.\n', pert);
            continue
        end

        ori_labels = fieldnames(dat_msg.rsq_by_method.(methods_both{1}));
        n_ori     = numel(ori_labels);
        n_bundles = numel(bundle_display);
        n_axes    = min(dat_msg.n_axes, dat_esg.n_axes);

        n_tests = numel(methods_both) * n_bundles * n_ori * n_axes;
        fprintf('  %s: %d methods × %d bundles × %d ori × %d axes = %d tests\n', ...
            pert, numel(methods_both), n_bundles, n_ori, n_axes, n_tests);

        pd_key = sprintf('mve_%s', pert);
        plot_data.(pd_key) = struct();

        for m = 1:numel(methods_both)
            method = methods_both{m};
            for b = 1:n_bundles
                mm = bund_msg == b;
                me = bund_esg == b;
                for o = 1:n_ori
                    ori = ori_labels{o};
                    for ax = 1:n_axes
                        rsq_msg_block = dat_msg.rsq_by_method.(method).(ori)(mm, :, ax);
                        rsq_esg_block = dat_esg.rsq_by_method.(method).(ori)(me, :, ax);

                        vals_msg = rsq_msg_block(:);
                        vals_esg = rsq_esg_block(:);
                        ok_msg = ~isnan(vals_msg);
                        ok_esg = ~isnan(vals_esg);
                        if sum(ok_msg) < 3 || sum(ok_esg) < 3; continue; end

                        [p_raw, ~] = ranksum(vals_msg(ok_msg), vals_esg(ok_esg));
                        p_corr     = min(1, p_raw * n_tests);

                        row = build_row('msg_vs_esg', pert, 'MSG', 'ESG', ...
                            method, b, bundle_display{b}, ori, ax, ...
                            min(sum(ok_msg), sum(ok_esg)), ...
                            vals_msg(ok_msg), vals_esg(ok_esg), p_raw, p_corr, alpha);
                        all_rows{end+1} = row; %#ok<AGROW>

                        fk = sprintf('%s_%s_ax%d', method, ori, ax);
                        if ~isfield(plot_data.(pd_key), fk)
                            plot_data.(pd_key).(fk) = struct('y', [], 'g_bundle', {{}}, 'g_cond', {{}});
                        end
                        pd = plot_data.(pd_key).(fk);
                        for b2 = 1:n_bundles
                            mm2 = bund_msg == b2;
                            me2 = bund_esg == b2;
                            vm2 = mean(dat_msg.rsq_by_method.(method).(ori)(mm2, :, ax), 1, 'omitnan');
                            ve2 = mean(dat_esg.rsq_by_method.(method).(ori)(me2, :, ax), 1, 'omitnan');
                            vm2 = vm2(~isnan(vm2));
                            ve2 = ve2(~isnan(ve2));
                            pd.y        = [pd.y; vm2(:); ve2(:)];
                            pd.g_bundle = [pd.g_bundle;
                                           repmat({bundle_display{b2}}, numel(vm2)+numel(ve2), 1)];
                            pd.g_cond   = [pd.g_cond;
                                           repmat({'MSG'}, numel(vm2), 1);
                                           repmat({'ESG'}, numel(ve2), 1)];
                        end
                        plot_data.(pd_key).(fk) = pd;
                    end
                end
            end
        end
        fprintf('  %s: done\n', pert);
    end
end

%% =========================================================================
%% SAVE WILCOXON/BONFERRONI TABLE (reference only — figures use permutation
%% results below)
%% =========================================================================

if isempty(all_rows)
    fprintf('\nNo comparisons completed — check input paths and flags.\n');
    return
end

results = vertcat(all_rows{:});
outfile = fullfile(save_dir, 'comparison_stats.tsv');
writetable(results, outfile, 'Delimiter', '\t', 'FileType', 'text');
fprintf('\nSaved: %s\n', outfile);

sig = results(results.significant == 1, :);
fprintf('\n=== Significant comparisons (Wilcoxon/Bonferroni): %d of %d ===\n', height(sig), height(results));
for r = 1:height(sig)
    dir_str = sig.effect_direction{r};
    fprintf('  [%s | %s | %s | B%d | %s | ax%d]\n', ...
        sig.comparison_type{r}, sig.modality_or_pert{r}, sig.method{r}, ...
        sig.bundle(r), sig.orientation{r}, sig.sensor_axis(r));
    fprintf('    %s median=%.4f  vs  %s median=%.4f  p_corr=%.4f  → %s\n', ...
        sig.label_A{r}, sig.median_A(r), sig.label_B{r}, sig.median_B(r), ...
        sig.p_corrected(r), dir_str);
end

%% =========================================================================
%% PERMUTATION-BASED ANALYSIS (cord-position + aggregate)
%% Source vs Sensor, AND MSG vs ESG
%% =========================================================================
%
% For every combination, and for every cord position within it, this runs a
% two-sample permutation test:
%
%   1. Observed effect = median(group A) - median(group B), using the
%      (typically 8) shift realisations in that bundle.
%   2. Null distribution: pool group A + group B, shuffle the labels
%      n_perm times, recompute the statistic each time.
%   3. p_perm = fraction of null values at least as extreme as the observed
%      effect — for the AGGREGATE test (one number per combo), this is used
%      directly, no correction needed.
%   4. For the CORD-POSITION test (one number per position within a combo),
%      a cluster correction is applied on top:
%        a. Threshold every position's null distribution at alpha_perm to
%           get a per-position exceedance cutoff.
%        b. In the observed map, group contiguous same-sign exceedances
%           into clusters; cluster "mass" = sum(|effect|) across the cluster.
%        c. In every one of the n_perm shuffled maps, apply the same
%           cutoffs, find clusters the same way, and keep only the single
%           largest cluster mass from that shuffle — building a null
%           distribution of "biggest cluster you'd see anywhere along the
%           cord by chance."
%        d. An observed cluster survives correction if its mass beats
%           1-alpha_perm of that null distribution.
%   5. A bootstrap resampling (n_boot draws, with replacement, within each
%      true group) gives a 90% interval around the observed effect for the
%      uncertainty ribbon on the cord-position figures.
%
% Group A / group B:
%   source_vs_sensor : A = source, B = sensor (within MSG, and within ESG)
%   msg_vs_esg       : A = MSG,    B = ESG    (within source, sensor, cond)
%
% Two outputs:
%   - cord_stats  (long format, one row per cord position) → cord_position_stats.tsv
%   - perm_stats  (one row per combo, all positions pooled) → permutation_stats.tsv
% perm_stats now drives ALL heatmaps, dumbbell plots, and the significance
% summary (both comparisons). comparison_stats.tsv above is reference only.
%
% NOTE on MSG vs ESG cord-position tests: if MSG and ESG have a different
% number of valid cord positions for a given combo, position-by-position
% comparison isn't meaningful (position p in one may not be the same
% physical location as position p in the other) — that combo's aggregate
% test still runs, but its cord-position test is skipped with a warning.

cord_rows = {};
perm_rows = {};
cord_plot_data = struct();

fprintf('\n--- Permutation-based analysis (n_perm=%d, n_boot=%d, alpha=%.3g) ---\n', ...
    n_perm, n_boot, alpha_perm);

%% --- Source vs Sensor -------------------------------------------------

if run_source_vs_sensor
    fprintf('\nSource vs Sensor:\n');

    for mod_idx = 1:2
        if mod_idx == 1
            dat = msg;  mod_label = 'MSG';
        else
            if ~have_esg; continue; end
            dat = esg;  mod_label = 'ESG';
        end
        if ~dat.have_source || ~dat.have_sensor; continue; end

        src  = dat.source;
        sens = dat.sensor;
        src_bund  = src.valid_source_bundle_idx;
        sens_bund = sens.valid_bundle_idx;

        methods_both = intersect(src.loaded_methods, sens.loaded_methods);
        if isempty(methods_both); continue; end

        ori_labels = fieldnames(src.rsq_by_method.(methods_both{1}));
        n_ori     = numel(ori_labels);
        n_bundles = numel(bundle_display);
        n_axes    = min(src.n_axes, sens.n_axes);

        pd_key = sprintf('sourcesensor_%s', mod_label);
        cord_plot_data.(pd_key) = struct();

        for m = 1:numel(methods_both)
            method = methods_both{m};
            for b = 1:n_bundles
                sm = src_bund  == b;
                ss = sens_bund == b;
                for o = 1:n_ori
                    ori = ori_labels{o};
                    for ax = 1:n_axes
                        X_A = src.rsq_by_method.(method).(ori)(sm, :, ax);
                        X_B = sens.rsq_by_method.(method).(ori)(ss, :, ax);

                        [ok_run, cord_rows, perm_rows, cord_plot_data.(pd_key)] = ...
                            run_permutation_combo(X_A, X_B, cord_position_mm, n_perm, ...
                            n_boot, alpha_perm, 'source_vs_sensor', mod_label, ...
                            'source', 'sensor', method, b, bundle_display{b}, ori, ax, ...
                            cord_rows, perm_rows, cord_plot_data.(pd_key));
                        if ~ok_run
                            fprintf(['  skipped %s/%s/%s/B%d/ax%d (all-NaN block)\n'], ...
                                mod_label, method, ori, b, ax); %#ok<*NBRAK2>
                        end
                    end
                end
            end
        end
        fprintf('  %s: done\n', mod_label);
    end
end

%% --- MSG vs ESG ---------------------------------------------------------

if run_msg_vs_esg && have_esg
    fprintf('\nMSG vs ESG:\n');

    pert_types = {'source', 'sensor', 'cond'};

    for pt = 1:numel(pert_types)
        pert = pert_types{pt};

        switch pert
            case 'source'
                have_both = msg.have_source && esg.have_source;
                if ~have_both; continue; end
                dat_msg = msg.source; dat_esg = esg.source;
                bund_msg = dat_msg.valid_source_bundle_idx;
                bund_esg = dat_esg.valid_source_bundle_idx;
            case 'sensor'
                have_both = msg.have_sensor && esg.have_sensor;
                if ~have_both; continue; end
                dat_msg = msg.sensor; dat_esg = esg.sensor;
                bund_msg = dat_msg.valid_bundle_idx;
                bund_esg = dat_esg.valid_bundle_idx;
            case 'cond'
                have_both = msg.have_cond && esg.have_cond;
                if ~have_both; continue; end
                dat_msg = msg.cond; dat_esg = esg.cond;
                bund_msg = dat_msg.valid_cond_bundle_idx;
                bund_esg = dat_esg.valid_cond_bundle_idx;
        end

        methods_both = intersect(dat_msg.loaded_methods, dat_esg.loaded_methods);
        if isempty(methods_both); continue; end

        ori_labels = fieldnames(dat_msg.rsq_by_method.(methods_both{1}));
        n_ori     = numel(ori_labels);
        n_bundles = numel(bundle_display);
        n_axes    = min(dat_msg.n_axes, dat_esg.n_axes);

        pd_key = sprintf('msgesg_%s', pert);
        cord_plot_data.(pd_key) = struct();

        for m = 1:numel(methods_both)
            method = methods_both{m};
            for b = 1:n_bundles
                mm = bund_msg == b;
                me = bund_esg == b;
                for o = 1:n_ori
                    ori = ori_labels{o};
                    for ax = 1:n_axes
                        X_A = dat_msg.rsq_by_method.(method).(ori)(mm, :, ax);
                        X_B = dat_esg.rsq_by_method.(method).(ori)(me, :, ax);

                        if size(X_A, 2) ~= size(X_B, 2)
                            warning(['MSG (%d positions) and ESG (%d positions) don''t ' ...
                                'match for %s/%s/%s/B%d/ax%d — skipping cord-position test, ' ...
                                'running aggregate only.'], size(X_A,2), size(X_B,2), ...
                                pert, method, ori, b, ax);
                            if all(isnan(X_A(:))) || all(isnan(X_B(:))); continue; end
                            [agg_diff, agg_p, agg_ci_lo, agg_ci_hi] = permutation_test_diff( ...
                                X_A(:), X_B(:), n_perm, n_boot);
                            perm_rows = append_perm_row(perm_rows, 'msg_vs_esg', pert, ...
                                'MSG', 'ESG', method, b, bundle_display{b}, ori, ax, ...
                                X_A, X_B, agg_diff, agg_p, agg_ci_lo, agg_ci_hi, alpha_perm);
                            continue
                        end

                        [ok_run, cord_rows, perm_rows, cord_plot_data.(pd_key)] = ...
                            run_permutation_combo(X_A, X_B, cord_position_mm, n_perm, ...
                            n_boot, alpha_perm, 'msg_vs_esg', pert, 'MSG', 'ESG', ...
                            method, b, bundle_display{b}, ori, ax, ...
                            cord_rows, perm_rows, cord_plot_data.(pd_key));
                        if ~ok_run
                            fprintf('  skipped %s/%s/%s/B%d/ax%d (all-NaN block)\n', ...
                                pert, method, ori, b, ax);
                        end
                    end
                end
            end
        end
        fprintf('  %s: done\n', pert);
    end
end

if ~isempty(cord_rows)
    cord_stats = vertcat(cord_rows{:});
    cord_outfile = fullfile(save_dir, 'cord_position_stats.tsv');
    writetable(cord_stats, cord_outfile, 'Delimiter', '\t', 'FileType', 'text');
    fprintf('\nSaved: %s\n', cord_outfile);
else
    cord_stats = table();
end

if ~isempty(perm_rows)
    perm_stats = vertcat(perm_rows{:});
    perm_outfile = fullfile(save_dir, 'permutation_stats.tsv');
    writetable(perm_stats, perm_outfile, 'Delimiter', '\t', 'FileType', 'text');
    fprintf('Saved: %s\n', perm_outfile);
else
    perm_stats = table();
end

%% =========================================================================
%% FIGURES
%% =========================================================================

fprintf('\nGenerating figures...\n');

%% Figure set 1: source vs sensor within MSG / ESG — descriptive box plots

for mod_idx = 1:2
    if mod_idx == 1; mod_label = 'MSG'; else; mod_label = 'ESG'; end
    pd_key = sprintf('sv_s_%s', mod_label);
    if ~isfield(plot_data, pd_key); continue; end

    fk_list = fieldnames(plot_data.(pd_key));
    if isempty(fk_list); continue; end

    methods_fig = {};  oris_fig = {};
    for fi = 1:numel(fk_list)
        parts = strsplit(fk_list{fi}, '_');
        methods_fig{end+1} = parts{1}; %#ok<AGROW>
        oris_fig{end+1}    = parts{2}; %#ok<AGROW>
    end
    methods_fig = unique(methods_fig);
    oris_fig    = unique(oris_fig);
    n_m = numel(methods_fig);
    n_o = numel(oris_fig);

    res_sub = results(strcmp(results.comparison_type, 'source_vs_sensor') & ...
                      strcmp(results.modality_or_pert, mod_label), :);

    fig = figure('Color', 'w', 'Position', [50 50 max(400, 380*n_m) max(360, 300*n_o)]);
    tl  = tiledlayout(n_o, n_m, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(tl, sprintf('%s — Source vs Sensor r² (descriptive)', mod_label), ...
        'FontSize', 12, 'FontWeight', 'bold');

    for o = 1:n_o
        ori = oris_fig{o};
        for m = 1:n_m
            method = methods_fig{m};
            ax_tile = nexttile(tl);
            hold(ax_tile, 'on');

            y_all = [];  gb_all = {};  gc_all = {};
            for ax = 1:10
                fk = sprintf('%s_%s_%s_ax%d', method, ori, mod_label, ax);
                if isfield(plot_data.(pd_key), fk)
                    pd = plot_data.(pd_key).(fk);
                    y_all  = [y_all;  pd.y]; %#ok<AGROW>
                    gb_all = [gb_all; pd.g_bundle]; %#ok<AGROW>
                    gc_all = [gc_all; pd.g_cond]; %#ok<AGROW>
                end
            end
            if isempty(y_all); continue; end

            bc = boxchart(ax_tile, categorical(gb_all, bundle_display), y_all, ...
                'GroupByColor', categorical(gc_all, {'Source', 'Sensor'}), ...
                'BoxFaceAlpha', 0.7, 'LineWidth', pub_lw, 'MarkerSize', 3);
            bc(1).BoxFaceColor = col_source;
            if numel(bc) > 1; bc(2).BoxFaceColor = col_sensor; end

            if show_significance_on_boxplots
                for b = 1:numel(bundle_display)
                    rr = res_sub(strcmp(res_sub.method, method) & ...
                                 res_sub.bundle == b & ...
                                 strcmp(res_sub.orientation, ori), :);
                    if isempty(rr); continue; end
                    p_c = min(rr.p_corrected);
                    if p_c < alpha
                        y_top = max(y_all(strcmp(gb_all, bundle_display{b}))) + 0.02;
                        plot(ax_tile, [b-0.3 b+0.3], [y_top y_top], 'k-', 'LineWidth', 0.8);
                        if p_c < 0.001; txt = '***'; elseif p_c < 0.01; txt = '**'; else; txt = '*'; end
                        text(ax_tile, b, y_top + 0.01, txt, 'HorizontalAlignment', 'center', 'FontSize', 9);
                    end
                end
            end

            ylim(ax_tile, [max(0, min(y_all)-0.05) 1]);
            title(ax_tile, sprintf('%s | %s', upper(method), ori), 'FontSize', 10);
            ylabel(ax_tile, 'r²', 'FontSize', 9);
            if o == 1 && m == 1; legend(ax_tile, 'Location', 'southwest', 'FontSize', 8); else; legend(ax_tile, 'off'); end
            box(ax_tile, 'on');
        end
    end

    fname = sprintf('source_vs_sensor_%s', lower(mod_label));
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

%% Figure set 2: MSG vs ESG per perturbation type — descriptive box plots

pert_info = {
    'source', 'Source shifts (mm)',  {col_MSG, col_ESG};
    'sensor', 'Sensor shifts (mm)',  {col_MSG, col_ESG};
    'cond',   'Conductivity (%%)',   {col_MSG, col_ESG}
};

for pt = 1:size(pert_info, 1)
    pert     = pert_info{pt, 1};
    pert_dsp = pert_info{pt, 2};
    cols     = pert_info{pt, 3};
    pd_key   = sprintf('mve_%s', pert);
    if ~isfield(plot_data, pd_key); continue; end

    fk_list = fieldnames(plot_data.(pd_key));
    if isempty(fk_list); continue; end

    methods_fig = {};  oris_fig = {};
    for fi = 1:numel(fk_list)
        parts = strsplit(fk_list{fi}, '_');
        methods_fig{end+1} = parts{1}; %#ok<AGROW>
        oris_fig{end+1}    = parts{2}; %#ok<AGROW>
    end
    methods_fig = unique(methods_fig);
    oris_fig    = unique(oris_fig);
    n_m = numel(methods_fig);
    n_o = numel(oris_fig);

    res_sub = results(strcmp(results.comparison_type, 'msg_vs_esg') & ...
                      strcmp(results.modality_or_pert, pert), :);

    fig = figure('Color', 'w', 'Position', [50 50 max(400, 380*n_m) max(360, 300*n_o)]);
    tl  = tiledlayout(n_o, n_m, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(tl, sprintf('MSG vs ESG — %s r² (descriptive) | %s', pert, pert_dsp), ...
        'FontSize', 12, 'FontWeight', 'bold');

    for o = 1:n_o
        ori = oris_fig{o};
        for m = 1:n_m
            method = methods_fig{m};
            ax_tile = nexttile(tl);
            hold(ax_tile, 'on');

            y_all = [];  gb_all = {};  gc_all = {};
            for ax = 1:10
                fk = sprintf('%s_%s_ax%d', method, ori, ax);
                if isfield(plot_data.(pd_key), fk)
                    pd = plot_data.(pd_key).(fk);
                    y_all  = [y_all;  pd.y]; %#ok<AGROW>
                    gb_all = [gb_all; pd.g_bundle]; %#ok<AGROW>
                    gc_all = [gc_all; pd.g_cond]; %#ok<AGROW>
                end
            end
            if isempty(y_all); continue; end

            bc = boxchart(ax_tile, categorical(gb_all, bundle_display), y_all, ...
                'GroupByColor', categorical(gc_all, {'MSG', 'ESG'}), ...
                'BoxFaceAlpha', 0.7, 'LineWidth', pub_lw, 'MarkerSize', 3);
            bc(1).BoxFaceColor = cols{1};
            if numel(bc) > 1; bc(2).BoxFaceColor = cols{2}; end

            if show_significance_on_boxplots
                for b = 1:numel(bundle_display)
                    rr = res_sub(strcmp(res_sub.method, method) & ...
                                 res_sub.bundle == b & ...
                                 strcmp(res_sub.orientation, ori), :);
                    if isempty(rr); continue; end
                    p_c = min(rr.p_corrected);
                    if p_c < alpha
                        y_top = max(y_all(strcmp(gb_all, bundle_display{b}))) + 0.02;
                        plot(ax_tile, [b-0.3 b+0.3], [y_top y_top], 'k-', 'LineWidth', 0.8);
                        if p_c < 0.001; txt = '***'; elseif p_c < 0.01; txt = '**'; else; txt = '*'; end
                        text(ax_tile, b, y_top + 0.01, txt, 'HorizontalAlignment', 'center', 'FontSize', 9);
                    end
                end
            end

            ylim(ax_tile, [max(0, min(y_all)-0.05) 1]);
            title(ax_tile, sprintf('%s | %s', upper(method), ori), 'FontSize', 10);
            ylabel(ax_tile, 'r²', 'FontSize', 9);
            if o == 1 && m == 1; legend(ax_tile, 'Location', 'southwest', 'FontSize', 8); else; legend(ax_tile, 'off'); end
            box(ax_tile, 'on');
        end
    end

    fname = sprintf('msg_vs_esg_%s', pert);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

if make_stat_figures

    %% Figure set 3: dumbbell plots — both comparisons, from perm_stats

    fprintf('\nGenerating dumbbell plots...\n');

    if run_source_vs_sensor && ~isempty(perm_stats)
        for mod_idx = 1:2
            if mod_idx == 1; mod_label = 'MSG'; else; mod_label = 'ESG'; end
            if mod_idx == 2 && ~have_esg; continue; end
            plot_dumbbell_group(perm_stats, 'source_vs_sensor', mod_label, save_dir, ...
                'source', 'sensor', col_source, col_sensor, 'source_vs_sensor');
        end
    end
    if run_msg_vs_esg && have_esg && ~isempty(perm_stats)
        for pt = {'source', 'sensor', 'cond'}
            plot_dumbbell_group(perm_stats, 'msg_vs_esg', pt{1}, save_dir, ...
                'MSG', 'ESG', col_MSG_spatial, col_ESG_spatial, 'msg_vs_esg');
        end
    end

    %% Figure set 4: effect-size heatmaps — both comparisons, from perm_stats

    fprintf('\nGenerating effect-size heatmaps...\n');

    cmap_sourcesensor = orangeblue_cmap(col_source, col_sensor);
    cmap_msgesg        = orangeblue_cmap(col_MSG_spatial, col_ESG_spatial);

    if run_source_vs_sensor && ~isempty(perm_stats)
        plot_effect_heatmap_group(perm_stats, 'source_vs_sensor', 'MSG', save_dir, ...
            cmap_sourcesensor, 'source_vs_sensor', bundle_display);
        if have_esg
            plot_effect_heatmap_group(perm_stats, 'source_vs_sensor', 'ESG', save_dir, ...
                cmap_sourcesensor, 'source_vs_sensor', bundle_display);
        end
    end
    if run_msg_vs_esg && have_esg && ~isempty(perm_stats)
        for pt = {'source', 'sensor', 'cond'}
            plot_effect_heatmap_group(perm_stats, 'msg_vs_esg', pt{1}, save_dir, ...
                cmap_msgesg, 'msg_vs_esg', bundle_display);
        end
    end

    %% Figure set 5: overall significance-rate summary (permutation-based,
    %% both comparisons — like-for-like now, both uncorrected-aggregate)

    fprintf('\nGenerating significance summary...\n');

    if ~isempty(perm_stats)
        plot_significance_summary(perm_stats, save_dir, alpha_perm);
    end

    %% Figure set 6: cord-position significance figures — both comparisons

    fprintf('\nGenerating cord-position significance figures...\n');

    if run_source_vs_sensor
        for mod_idx = 1:2
            if mod_idx == 1; mod_label = 'MSG'; else; mod_label = 'ESG'; end
            if mod_idx == 2 && ~have_esg; continue; end
            pd_key = sprintf('sourcesensor_%s', mod_label);
            plot_cord_position_significance(cord_plot_data, pd_key, ...
                sprintf('%s source vs sensor', mod_label), bundle_display, save_dir, ...
                col_source, col_sensor, sprintf('sourcesensor_%s', lower(mod_label)));
        end
    end
    if run_msg_vs_esg && have_esg
        for pt = {'source', 'sensor', 'cond'}
            pd_key = sprintf('msgesg_%s', pt{1});
            plot_cord_position_significance(cord_plot_data, pd_key, ...
                sprintf('MSG vs ESG (%s)', pt{1}), bundle_display, save_dir, ...
                col_MSG_spatial, col_ESG_spatial, sprintf('msgesg_%s', pt{1}));
        end
    end

end

fprintf('\npt_compare_perturbations complete.\n');
fprintf('Outputs: %s\n', save_dir);


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function dat = load_pert_results(base_path)
    dat.have_source = false;
    dat.have_sensor = false;
    dat.have_cond   = false;

    f = fullfile(base_path, 'pert_source_rsq.mat');
    if isfile(f); dat.source = load(f); dat.have_source = true; end

    f = fullfile(base_path, 'pert_sensor_rsq.mat');
    if isfile(f); dat.sensor = load(f); dat.have_sensor = true; end

    f = fullfile(base_path, 'pert_cond_rsq.mat');
    if isfile(f); dat.cond = load(f); dat.have_cond = true; end
end

function row = build_row(comp_type, mod_or_pert, label_A, label_B, ...
    method, b, b_name, ori, ax, n_pos, vals_A, vals_B, p_raw, p_corr, alpha)

    med_A  = median(vals_A, 'omitnan');
    med_B  = median(vals_B, 'omitnan');
    iqr_A  = iqr(vals_A);
    iqr_B  = iqr(vals_B);
    delta  = med_A - med_B;
    sig    = double(p_corr < alpha);
    if sig
        if delta < 0; dir = 'A_lower_rsq'; else; dir = 'B_lower_rsq'; end
    else
        dir = 'no_significant_difference';
    end

    row = table( ...
        {comp_type}, {mod_or_pert}, {label_A}, {label_B}, {method}, ...
        b, {b_name}, {ori}, ax, n_pos, ...
        med_A, med_B, delta, iqr_A, iqr_B, ...
        p_raw, p_corr, sig, {dir}, ...
        'VariableNames', { ...
            'comparison_type', 'modality_or_pert', 'label_A', 'label_B', 'method', ...
            'bundle', 'bundle_name', 'orientation', 'sensor_axis', 'n_positions', ...
            'median_A', 'median_B', 'delta_median_AminusB', 'iqr_A', 'iqr_B', ...
            'p_raw', 'p_corrected', 'significant', 'effect_direction' ...
        });
end

function s = yesno(v)
    if v; s = 'found'; else; s = 'not found'; end
end

function [obs_diff, p_perm, ci_lo, ci_hi] = permutation_test_diff(X_A, X_B, n_perm, n_boot)
    % Two-sample permutation test + bootstrap CI, per column of X_A / X_B.
    % Rows = samples within a group (e.g. shift realisations); columns =
    % independent things to test (pass column vectors for one aggregate
    % test). X_A and X_B may have different row counts.
    [nA, n_col] = size(X_A);
    nB = size(X_B, 1);
    n_tot = nA + nB;

    obs_diff = median(X_A, 1, 'omitnan') - median(X_B, 1, 'omitnan');
    pooled = [X_A; X_B];

    null_diff = nan(n_perm, n_col);
    for it = 1:n_perm
        [~, idx]  = sort(rand(n_tot, n_col), 1);
        lin_idx   = idx + (0:n_col-1) * n_tot;
        perm_pool = pooled(lin_idx);
        null_diff(it, :) = median(perm_pool(1:nA, :), 1, 'omitnan') - ...
                            median(perm_pool(nA+1:end, :), 1, 'omitnan');
    end
    p_perm = (sum(abs(null_diff) >= abs(obs_diff), 1) + 1) / (n_perm + 1);

    boot_diff = nan(n_boot, n_col);
    for it = 1:n_boot
        idxA = randi(nA, nA, 1);
        idxB = randi(nB, nB, 1);
        boot_diff(it, :) = median(X_A(idxA, :), 1, 'omitnan') - median(X_B(idxB, :), 1, 'omitnan');
    end
    ci_lo = prctile(boot_diff, 5, 1);
    ci_hi = prctile(boot_diff, 95, 1);
end

function [obs_diff, p_perm, ci_lo, ci_hi, sig_uncorrected, sig_cluster] = ...
    permutation_test_diff_cluster(X_A, X_B, n_perm, n_boot, alpha_perm)
    % As permutation_test_diff, but additionally returns cluster-corrected
    % significance across columns (intended for column = cord position).
    %
    % Cluster correction (Maris & Oostenveld style, adapted to 1-D):
    %   1. Threshold every column's null distribution at alpha_perm
    %      (two-tailed) → a per-position exceedance cutoff.
    %   2. In the observed map, group contiguous same-sign exceedances into
    %      clusters; cluster mass = sum(|effect|) over the cluster.
    %   3. In every shuffled map (each row of the null distribution), apply
    %      the SAME cutoffs, find clusters the same way, and keep the
    %      single largest cluster mass from that shuffle.
    %   4. An observed cluster survives correction if its mass beats
    %      1-alpha_perm of the "largest shuffled cluster" null distribution.
    %
    % sig_uncorrected : per-position significance, no spatial correction
    % sig_cluster      : per-position significance after cluster correction
    %                    (a stricter subset of sig_uncorrected)

    [nA, n_col] = size(X_A);
    nB = size(X_B, 1);
    n_tot = nA + nB;

    obs_diff = median(X_A, 1, 'omitnan') - median(X_B, 1, 'omitnan');
    pooled = [X_A; X_B];

    null_diff = nan(n_perm, n_col);
    for it = 1:n_perm
        [~, idx]  = sort(rand(n_tot, n_col), 1);
        lin_idx   = idx + (0:n_col-1) * n_tot;
        perm_pool = pooled(lin_idx);
        null_diff(it, :) = median(perm_pool(1:nA, :), 1, 'omitnan') - ...
                            median(perm_pool(nA+1:end, :), 1, 'omitnan');
    end
    p_perm = (sum(abs(null_diff) >= abs(obs_diff), 1) + 1) / (n_perm + 1);

    boot_diff = nan(n_boot, n_col);
    for it = 1:n_boot
        idxA = randi(nA, nA, 1);
        idxB = randi(nB, nB, 1);
        boot_diff(it, :) = median(X_A(idxA, :), 1, 'omitnan') - median(X_B(idxB, :), 1, 'omitnan');
    end
    ci_lo = prctile(boot_diff, 5, 1);
    ci_hi = prctile(boot_diff, 95, 1);

    thresh_lo = prctile(null_diff, 100*alpha_perm/2, 1);
    thresh_hi = prctile(null_diff, 100*(1-alpha_perm/2), 1);
    sig_uncorrected = p_perm < alpha_perm;

    max_cluster_null = zeros(n_perm, 1);
    for it = 1:n_perm
        row = null_diff(it, :);
        mask = row < thresh_lo | row > thresh_hi;
        runs = find_sig_runs(mask, sign(row));
        m = 0;
        for r = 1:size(runs, 1)
            m = max(m, sum(abs(row(runs(r,1):runs(r,2)))));
        end
        max_cluster_null(it) = m;
    end

    sig_cluster = false(1, n_col);
    obs_runs = find_sig_runs(sig_uncorrected, sign(obs_diff));
    for r = 1:size(obs_runs, 1)
        i0 = obs_runs(r,1); i1 = obs_runs(r,2);
        cmass = sum(abs(obs_diff(i0:i1)));
        p_cluster = (sum(max_cluster_null >= cmass) + 1) / (n_perm + 1);
        if p_cluster < alpha_perm
            sig_cluster(i0:i1) = true;
        end
    end
end

function perm_rows = append_perm_row(perm_rows, comp_type, mop, label_A, label_B, ...
    method, b, b_name, ori, ax, X_A, X_B, agg_diff, agg_p, agg_ci_lo, agg_ci_hi, alpha_perm)
    % Appends one aggregate-level permutation result row (used for heatmaps
    % / dumbbell plots).
    agg_sig = agg_p < alpha_perm;
    if agg_sig
        if agg_diff < 0; agg_dir = 'A_lower_rsq'; else; agg_dir = 'B_lower_rsq'; end
    else
        agg_dir = 'no_significant_difference';
    end
    perm_rows{end+1} = table( ...
        {comp_type}, {mop}, {label_A}, {label_B}, {method}, ...
        b, {b_name}, {ori}, ax, numel(X_A), ...
        median(X_A(:), 'omitnan'), median(X_B(:), 'omitnan'), agg_diff, ...
        iqr(X_A(:)), iqr(X_B(:)), agg_p, double(agg_sig), {agg_dir}, ...
        agg_ci_lo, agg_ci_hi, ...
        'VariableNames', {'comparison_type', 'modality_or_pert', 'label_A', ...
        'label_B', 'method', 'bundle', 'bundle_name', 'orientation', ...
        'sensor_axis', 'n_samples', 'median_A', 'median_B', ...
        'delta_median_AminusB', 'iqr_A', 'iqr_B', 'p_perm', 'significant', ...
        'effect_direction', 'boot_ci_lo', 'boot_ci_hi'}); %#ok<AGROW>
end

function [ok, cord_rows, perm_rows, cord_pd] = run_permutation_combo(X_A, X_B, ...
    cord_position_mm, n_perm, n_boot, alpha_perm, comp_type, mop, label_A, label_B, ...
    method, b, b_name, ori, ax, cord_rows, perm_rows, cord_pd)
    % Runs both the cord-position (cluster-corrected) and aggregate
    % permutation tests for one combo, appending results into the running
    % cord_rows / perm_rows cell arrays and the cord_pd plotting struct.
    % Returns ok=false (no-op) if the block is entirely NaN (e.g. a
    % degenerate axis/orientation combination for a given method).

    if all(isnan(X_A(:))) || all(isnan(X_B(:)))
        ok = false;
        return
    end
    ok = true;

    n_pos = size(X_A, 2);
    if n_pos == numel(cord_position_mm)
        pos_mm = cord_position_mm;
    else
        warning(['cord_position_mm length (%d) does not match n_pos (%d) for ' ...
            '%s/%s/%s/%s/B%d/ax%d — using position index instead. Update ' ...
            'cord_position_mm in the config section.'], numel(cord_position_mm), ...
            n_pos, comp_type, mop, method, ori, b, ax);
        pos_mm = 1:n_pos;
    end

    % --- cord position (cluster-corrected) ---
    [obs_diff, p_perm, ci_lo, ci_hi, sig_unc, sig_clu] = permutation_test_diff_cluster( ...
        X_A, X_B, n_perm, n_boot, alpha_perm);

    pos_dir = repmat({'not_significant'}, 1, n_pos);
    pos_dir(sig_clu & obs_diff < 0) = {'A_lower_rsq'};
    pos_dir(sig_clu & obs_diff > 0) = {'B_lower_rsq'};

    fk = sprintf('%s_%s_ax%d_B%d', method, ori, ax, b);
    cord_pd.(fk) = struct('position_mm', pos_mm, 'obs_diff', obs_diff, 'p_perm', p_perm, ...
        'ci_lo', ci_lo, 'ci_hi', ci_hi, 'sig_uncorrected', sig_unc, 'sig_cluster', sig_clu, ...
        'method', method, 'orientation', ori, 'sensor_axis', ax, 'bundle', b);

    for p = 1:n_pos
        cord_rows{end+1} = table( ...
            {comp_type}, {mop}, {label_A}, {label_B}, {method}, b, {b_name}, {ori}, ax, ...
            pos_mm(p), obs_diff(p), p_perm(p), ci_lo(p), ci_hi(p), ...
            double(sig_unc(p)), double(sig_clu(p)), pos_dir(p), ...
            'VariableNames', {'comparison_type', 'modality_or_pert', 'label_A', 'label_B', ...
            'method', 'bundle', 'bundle_name', 'orientation', 'sensor_axis', 'position_mm', ...
            'obs_diff_AminusB', 'p_perm', 'ci_lo', 'ci_hi', 'significant_uncorrected', ...
            'significant_cluster_corrected', 'effect_direction'}); %#ok<AGROW>
    end

    % --- aggregate (all shifts x all positions pooled) ---
    [agg_diff, agg_p, agg_ci_lo, agg_ci_hi] = permutation_test_diff(X_A(:), X_B(:), n_perm, n_boot);
    perm_rows = append_perm_row(perm_rows, comp_type, mop, label_A, label_B, method, b, ...
        b_name, ori, ax, X_A, X_B, agg_diff, agg_p, agg_ci_lo, agg_ci_hi, alpha_perm);
end

function plot_dumbbell_group(stats_table, comp_type, mop, save_dir, label_A, label_B, ...
    col_A, col_B, fname_prefix)
    % Dumbbell chart: for each (orientation, sensor axis) combination within
    % a bundle, shows median_A and median_B as connected dots. Filled =
    % significant (permutation, uncorrected aggregate test), hollow = not.
    % One figure per bundle size, tiled by method.

    sub = stats_table(strcmp(stats_table.comparison_type, comp_type) & ...
                       strcmp(stats_table.modality_or_pert, mop), :);
    if isempty(sub); return; end

    methods_fig  = unique(sub.method, 'stable');
    bundle_names = unique(sub.bundle_name, 'stable');
    n_m = numel(methods_fig);

    for bi = 1:numel(bundle_names)
        bname = bundle_names{bi};
        rows  = sub(strcmp(sub.bundle_name, bname), :);
        if isempty(rows); continue; end
        [~, sort_idx] = sort(abs(rows.delta_median_AminusB), 'ascend');
        rows = rows(sort_idx, :);

        n_r_max = 0;
        for m = 1:n_m
            n_r_max = max(n_r_max, sum(strcmp(rows.method, methods_fig{m})));
        end

        fig = figure('Color', 'w', 'Position', ...
            [50 50 max(420, 380*n_m) max(360, 24*n_r_max + 100)]);
        tl  = tiledlayout(1, n_m, 'TileSpacing', 'compact', 'Padding', 'normal');
        title(tl, sprintf('%s (%s) — %s bundle: %s vs %s median r²', ...
            strrep(comp_type, '_', ' '), mop, bname, label_A, label_B), ...
            'FontSize', 12, 'FontWeight', 'bold');

        for m = 1:n_m
            method = methods_fig{m};
            mrows  = rows(strcmp(rows.method, method), :);
            n_r    = height(mrows);
            ax_t   = nexttile(tl);
            hold(ax_t, 'on');

            labels = cell(n_r, 1);
            for r = 1:n_r
                labels{r} = sprintf('%s | ax%d', mrows.orientation{r}, mrows.sensor_axis(r));

                y = r;
                is_sig = mrows.significant(r) == 1;
                plot(ax_t, [mrows.median_A(r) mrows.median_B(r)], [y y], ...
                    'Color', [0.7 0.7 0.7], 'LineWidth', 1);
                if is_sig; faceA = col_A; faceB = col_B; else; faceA = 'none'; faceB = 'none'; end
                scatter(ax_t, mrows.median_A(r), y, 55, 'MarkerFaceColor', faceA, ...
                    'MarkerEdgeColor', col_A, 'LineWidth', 1.2);
                scatter(ax_t, mrows.median_B(r), y, 55, 'MarkerFaceColor', faceB, ...
                    'MarkerEdgeColor', col_B, 'LineWidth', 1.2);
            end

            set(ax_t, 'YTick', 1:n_r, 'YTickLabel', labels, 'FontSize', 8);
            ylim(ax_t, [0.5 n_r+0.5]);
            xlabel(ax_t, 'median r²', 'FontSize', 9);
            title(ax_t, upper(method), 'FontSize', 10);
            box(ax_t, 'on');
        end

        fname = sprintf('dumbbell_%s_%s_%s', fname_prefix, lower(mop), lower(bname));
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end

function plot_effect_heatmap_group(stats_table, comp_type, mop, save_dir, cmap, ...
    fname_prefix, bundle_display)
    % Orientation x bundle-size heatmap of effect size, tiled by method
    % (columns) x sensor axis (rows). Significant cells marked '*'.

    sub = stats_table(strcmp(stats_table.comparison_type, comp_type) & ...
                      strcmp(stats_table.modality_or_pert, mop), :);
    if isempty(sub); return; end

    methods_fig = unique(sub.method, 'stable');
    n_m     = numel(methods_fig);
    n_axes  = max(sub.sensor_axis);
    ori_order = unique(sub.orientation, 'stable');
    n_ori   = numel(ori_order);
    n_bund  = numel(bundle_display);

    clim_max = max(abs(sub.delta_median_AminusB), [], 'omitnan');
    if isnan(clim_max) || clim_max == 0; clim_max = 0.01; end

    fig = figure('Color', 'w', 'Position', [50 50 max(420, 300*n_m) max(300, 220*n_axes)]);
    tl  = tiledlayout(n_axes, n_m, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(tl, sprintf('%s (%s) — effect size heatmap (* = significant)', ...
        strrep(comp_type, '_', ' '), mop), 'FontSize', 12, 'FontWeight', 'bold');

    last_tile = [];
    for ax_i = 1:n_axes
        for m = 1:n_m
            method = methods_fig{m};
            tile = nexttile(tl);
            last_tile = tile;

            grid_delta = nan(n_ori, n_bund);
            grid_sig   = false(n_ori, n_bund);

            for o = 1:n_ori
                for b = 1:n_bund
                    rr = sub(strcmp(sub.method, method) & sub.sensor_axis == ax_i & ...
                              strcmp(sub.orientation, ori_order{o}) & sub.bundle == b, :);
                    if isempty(rr); continue; end
                    grid_delta(o, b) = rr.delta_median_AminusB(1);
                    grid_sig(o, b)   = rr.significant(1) == 1;
                end
            end

            imagesc(tile, grid_delta, [-clim_max clim_max]);
            colormap(tile, cmap);
            set(tile, 'XTick', 1:n_bund, 'XTickLabel', bundle_display, ...
                      'YTick', 1:n_ori, 'YTickLabel', ori_order, 'YDir', 'normal');
            title(tile, sprintf('%s | axis %d', upper(method), ax_i), 'FontSize', 10);

            for o = 1:n_ori
                for b = 1:n_bund
                    if isnan(grid_delta(o,b)); continue; end
                    txt = sprintf('%.3f', grid_delta(o,b));
                    if grid_sig(o,b); txt = [txt '*']; end %#ok<AGROW>
                    text(tile, b, o, txt, 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0 0 0]);
                end
            end
        end
    end

    if ~isempty(last_tile)
        cb = colorbar(last_tile);
        cb.Layout.Tile = 'east';
        ylabel(cb, 'median_A - median_B', 'FontSize', 9);
    end

    fname = sprintf('%s_heatmap_%s', fname_prefix, lower(mop));
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

function plot_significance_summary(perm_stats, save_dir, alpha_perm)
    % Proportion of aggregate permutation tests significant per orientation
    % and bundle size, one panel per comparison type (source_vs_sensor,
    % msg_vs_esg) — both computed the same way now, so directly comparable.

    comp_types = unique(perm_stats.comparison_type, 'stable');
    n_c = numel(comp_types);
    ori_order    = unique(perm_stats.orientation, 'stable');
    bundle_order = unique(perm_stats.bundle_name, 'stable');

    fig = figure('Color', 'w', 'Position', [50 50 max(420, 380*n_c) 380]);
    tl = tiledlayout(1, n_c, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(tl, sprintf('Proportion of aggregate tests significant (permutation, \x03b1=%.3g, uncorrected)', ...
        alpha_perm), 'FontSize', 12, 'FontWeight', 'bold');

    for c = 1:n_c
        sub  = perm_stats(strcmp(perm_stats.comparison_type, comp_types{c}), :);
        ax_t = nexttile(tl);
        hold(ax_t, 'on');

        prop = nan(numel(ori_order), numel(bundle_order));
        for o = 1:numel(ori_order)
            for b = 1:numel(bundle_order)
                rr = sub(strcmp(sub.orientation, ori_order{o}) & ...
                          strcmp(sub.bundle_name, bundle_order{b}), :);
                if isempty(rr); continue; end
                prop(o,b) = mean(rr.significant == 1);
            end
        end

        bar(ax_t, prop, 'grouped');
        set(ax_t, 'XTick', 1:numel(ori_order), 'XTickLabel', ori_order);
        ylim(ax_t, [0 1]);
        ylabel(ax_t, 'Proportion significant', 'FontSize', 9);
        title(ax_t, strrep(comp_types{c}, '_', ' '), 'FontSize', 10);
        if c == n_c; legend(ax_t, bundle_order, 'Location', 'eastoutside', 'FontSize', 8); end
        box(ax_t, 'on');
    end

    fname = 'significance_summary';
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

function plot_cord_position_significance(cord_plot_data, pd_key, title_label, ...
    bundle_display, save_dir, col_A, col_B, fname_prefix)
    % One figure per (method, bundle), tiled by orientation x sensor axis.
    % Each tile: observed effect (A - B) along the cord, a bootstrap
    % uncertainty ribbon, and two layers of shading —
    %   faint  = significant before cluster correction, but not after
    %            ("would have been flagged without correcting for the
    %            number of cord positions tested")
    %   solid  = still significant after cluster correction
    % col_A shading = A (source / MSG) lower there; col_B shading = B
    % (sensor / ESG) lower there.

    if ~isfield(cord_plot_data, pd_key); return; end
    pd = cord_plot_data.(pd_key);
    fk_list = fieldnames(pd);
    if isempty(fk_list); return; end

    methods_fig = {};  oris_fig = {};  n_axes = 0;
    for fi = 1:numel(fk_list)
        d = pd.(fk_list{fi});
        methods_fig{end+1} = d.method; %#ok<AGROW>
        oris_fig{end+1}    = d.orientation; %#ok<AGROW>
        n_axes = max(n_axes, d.sensor_axis);
    end
    methods_fig = unique(methods_fig, 'stable');
    oris_fig    = unique(oris_fig, 'stable');
    n_o = numel(oris_fig);

    for m = 1:numel(methods_fig)
        method = methods_fig{m};
        for b = 1:numel(bundle_display)
            fig = figure('Color', 'w', 'Position', ...
                [50 50 max(500, 320*n_axes) max(400, 260*n_o)]);
            tl = tiledlayout(n_o, n_axes, 'TileSpacing', 'compact', 'Padding', 'normal');
            title(tl, sprintf(['%s — %s %s bundle: effect along cord ' ...
                '(permutation test, cluster-corrected)'], title_label, upper(method), ...
                bundle_display{b}), 'FontSize', 12, 'FontWeight', 'bold');

            any_tile = false;
            for o = 1:n_o
                ori = oris_fig{o};
                for ax = 1:n_axes
                    fk = sprintf('%s_%s_ax%d_B%d', method, ori, ax, b);
                    ax_t = nexttile(tl);
                    hold(ax_t, 'on');
                    if ~isfield(pd, fk)
                        axis(ax_t, 'off');
                        continue
                    end
                    any_tile = true;
                    d = pd.(fk);

                    fill(ax_t, [d.position_mm fliplr(d.position_mm)], ...
                        [d.ci_lo fliplr(d.ci_hi)], [0.5 0.5 0.5], ...
                        'FaceAlpha', 0.25, 'EdgeColor', 'none');
                    plot(ax_t, d.position_mm, d.obs_diff, 'k-', 'LineWidth', 1.2);
                    yline(ax_t, 0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':');

                    ylim(ax_t, [min([d.ci_lo -0.01]) max([d.ci_hi 0.01])]);
                    xlim(ax_t, [min(d.position_mm) max(d.position_mm)]);
                    yl = ylim(ax_t);

                    % Faint: uncorrected-only (didn't survive cluster correction)
                    unc_only = d.sig_uncorrected & ~d.sig_cluster;
                    runs = find_sig_runs(unc_only, sign(d.obs_diff));
                    for r = 1:size(runs, 1)
                        i0 = runs(r,1); i1 = runs(r,2); s = runs(r,3);
                        if s < 0; c = col_A; else; c = col_B; end
                        xs = [d.position_mm(i0) d.position_mm(i1) d.position_mm(i1) d.position_mm(i0)];
                        ph = patch(ax_t, xs, [yl(1) yl(1) yl(2) yl(2)], c, ...
                            'FaceAlpha', 0.06, 'EdgeColor', 'none');
                        uistack(ph, 'bottom');
                    end

                    % Solid: cluster-corrected significant
                    runs = find_sig_runs(d.sig_cluster, sign(d.obs_diff));
                    for r = 1:size(runs, 1)
                        i0 = runs(r,1); i1 = runs(r,2); s = runs(r,3);
                        if s < 0; c = col_A; else; c = col_B; end
                        xs = [d.position_mm(i0) d.position_mm(i1) d.position_mm(i1) d.position_mm(i0)];
                        ph = patch(ax_t, xs, [yl(1) yl(1) yl(2) yl(2)], c, ...
                            'FaceAlpha', 0.18, 'EdgeColor', 'none');
                        uistack(ph, 'bottom');
                    end

                    title(ax_t, sprintf('%s | axis %d', ori, ax), 'FontSize', 10);
                    xlabel(ax_t, 'Cord position (mm)', 'FontSize', 8);
                    ylabel(ax_t, '\x0394 r²', 'FontSize', 8);
                    box(ax_t, 'on');
                end
            end

            if ~any_tile; close(fig); continue; end

            fname = sprintf('cord_position_%s_%s_bundle%d_%s', ...
                fname_prefix, lower(method), b, lower(bundle_display{b}));
            exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
            saveas(fig, fullfile(save_dir, [fname '.fig']));
            close(fig);
            fprintf('  Saved: %s\n', fname);
        end
    end
end

function runs = find_sig_runs(sig, s)
    % Returns an [N x 3] matrix of [start_idx end_idx sign] for each
    % contiguous run where sig is true and sign(effect) is constant.
    runs = [];
    n = numel(sig);
    i = 1;
    while i <= n
        if ~sig(i); i = i + 1; continue; end
        j = i;
        while j < n && sig(j+1) && s(j+1) == s(i)
            j = j + 1;
        end
        runs = [runs; i j s(i)]; %#ok<AGROW>
        i = j + 1;
    end
end

function cmap = orangeblue_cmap(col_A, col_B)
    % Diverging colormap: negative effect (A lower) shades toward col_A;
    % positive effect (B lower) shades toward col_B. White at zero.
    n = 128;
    lo = [linspace(col_A(1),1,n)' linspace(col_A(2),1,n)' linspace(col_A(3),1,n)'];
    hi = [linspace(1,col_B(1),n)' linspace(1,col_B(2),n)' linspace(1,col_B(3),n)'];
    cmap = [lo; hi];
end
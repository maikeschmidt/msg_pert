% pt_compare_perturbations - Cross-perturbation and cross-modality r² comparison
%
% Compares r² distributions across:
%   1. Perturbation types within the same modality — source vs sensor shifts
%      at matched bundle magnitudes (same mm range). Paired Wilcoxon signed-rank
%      test, pairing by source position (mean over shifts per bundle).
%
%   2. MSG vs ESG within the same perturbation type and forward model method.
%      Unpaired Wilcoxon rank-sum test (different sensor geometry).
%      Includes source, sensor, and conductivity (cond only compared MSG/ESG,
%      not across perturbation types since cond uses % units, not mm).
%
% Bonferroni correction is applied within each comparison group (not globally).
%
% OUTPUTS (saved to <msg_results_path>/perturbation_analysis/comparison/):
%   comparison_stats.tsv                — full results table
%   source_vs_sensor_<modality>.png/fig — grouped box plots per modality
%   msg_vs_esg_<pert_type>.png/fig      — grouped box plots per pert type
%
% DEPENDENCIES:
%   pert_source_rsq.mat, pert_sensor_rsq.mat, pert_cond_rsq.mat
%   (from both MSG and ESG analysis directories when run_msg_vs_esg = true)
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
msg_results_path = 'D:\Simulations\Pertubations\fields';   % SET THIS

% Path to ESG analysis results (leave '' to skip MSG vs ESG comparisons)
esg_results_path = '';   % SET THIS

% Which comparisons to run
run_source_vs_sensor = true;   % source vs sensor (matched bundle, mm-based)
run_msg_vs_esg       = true;   % MSG vs ESG per perturbation type

% Significance level (after Bonferroni correction within each group)
alpha = 0.05;

% Output directory
save_dir = fullfile(msg_results_path, 'perturbation_analysis', 'comparison');

% Bundle labels (must match the three-bundle structure in config_pert)
bundle_display = {'Small', 'Medium', 'Large'};

% Figure colours: source, sensor, cond / MSG, ESG
col_source = [0.25 0.50 0.80];
col_sensor = [0.85 0.38 0.25];
col_cond   = [0.35 0.70 0.45];
col_MSG    = [0.25 0.50 0.80];
col_ESG    = [0.80 0.50 0.20];

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
%% COLLECT RESULTS TABLE AND RAW PLOT DATA
%% =========================================================================

all_rows  = {};    % cell of row tables, concatenated at end
plot_data = struct();   % stores raw per-position means for figures

%% -------------------------------------------------------------------------
%% COMPARISON 1 — Source vs Sensor within each modality
%% -------------------------------------------------------------------------

if run_source_vs_sensor
    fprintf('\n--- Source vs Sensor comparison ---\n');

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

        src  = dat.source;   % pert_source_rsq fields
        sens = dat.sensor;   % pert_sensor_rsq fields

        % source uses 'valid_source_bundle_idx'; sensor uses 'valid_bundle_idx'
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

        pd_key = sprintf('sv_s_%s', mod_label);   % plot_data key
        plot_data.(pd_key) = struct();

        for m = 1:numel(methods_both)
            method = methods_both{m};
            for b = 1:n_bundles
                sm = src_bund  == b;
                ss = sens_bund == b;
                for o = 1:n_ori
                    ori = ori_labels{o};
                    for ax = 1:n_axes
                        % Mean across shifts per source position → [n_src] vectors
                        src_rsq_block  = src.rsq_by_method.(method).(ori)(sm, :, ax);
                        sens_rsq_block = sens.rsq_by_method.(method).(ori)(ss, :, ax);

                        src_pos_mean  = mean(src_rsq_block,  1, 'omitnan');  % [1×n_src]
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

                        % Store raw data for figures
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
%% COMPARISON 2 — MSG vs ESG per perturbation type
%% -------------------------------------------------------------------------

if run_msg_vs_esg && have_esg
    fprintf('\n--- MSG vs ESG comparison ---\n');

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
                        % Pool all shift × position values (MSG and ESG may differ in n_src)
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

                        % Store raw for figures (per-position means for cleaner display)
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
%% SAVE TABLE
%% =========================================================================

if isempty(all_rows)
    fprintf('\nNo comparisons completed — check input paths and flags.\n');
    return
end

results = vertcat(all_rows{:});
outfile = fullfile(save_dir, 'comparison_stats.tsv');
writetable(results, outfile, 'Delimiter', '\t', 'FileType', 'text');
fprintf('\nSaved: %s\n', outfile);

% Print significant results
sig = results(results.significant == 1, :);
fprintf('\n=== Significant comparisons: %d of %d ===\n', height(sig), height(results));
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
%% FIGURES
%% =========================================================================

fprintf('\nGenerating figures...\n');

%% Figure set 1: source vs sensor within MSG / ESG

for mod_idx = 1:2
    if mod_idx == 1; mod_label = 'MSG'; else; mod_label = 'ESG'; end
    pd_key = sprintf('sv_s_%s', mod_label);
    if ~isfield(plot_data, pd_key); continue; end

    fk_list = fieldnames(plot_data.(pd_key));
    if isempty(fk_list); continue; end

    % Parse keys to find unique methods and oris
    methods_fig = {};  oris_fig = {};
    for fi = 1:numel(fk_list)
        parts = strsplit(fk_list{fi}, '_');
        % key format: <method>_<ori>_<modality>_ax<N>
        methods_fig{end+1} = parts{1}; %#ok<AGROW>
        oris_fig{end+1}    = parts{2}; %#ok<AGROW>
    end
    methods_fig = unique(methods_fig);
    oris_fig    = unique(oris_fig);
    n_m = numel(methods_fig);
    n_o = numel(oris_fig);

    % Get significance lookup from results table
    res_sub = results(strcmp(results.comparison_type, 'source_vs_sensor') & ...
                      strcmp(results.modality_or_pert, mod_label), :);

    fig = figure('Color', 'w', 'Position', [50 50 max(400, 380*n_m) max(360, 300*n_o)]);
    tl  = tiledlayout(n_o, n_m, 'TileSpacing', 'compact', 'Padding', 'normal');
    title(tl, sprintf('%s — Source vs Sensor r²  (paired signed-rank, Bonferroni)', mod_label), ...
        'FontSize', 12, 'FontWeight', 'bold');

    for o = 1:n_o
        ori = oris_fig{o};
        for m = 1:n_m
            method = methods_fig{m};
            ax_tile = nexttile(tl);
            hold(ax_tile, 'on');

            % Concatenate across sensor axes
            y_all = [];  gb_all = {};  gc_all = {};
            for ax = 1:10   % generous upper bound
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

            % Significance annotations
            for b = 1:numel(bundle_display)
                rr = res_sub(strcmp(res_sub.method, method) & ...
                             res_sub.bundle == b & ...
                             strcmp(res_sub.orientation, ori), :);
                if isempty(rr); continue; end
                p_c = min(rr.p_corrected);
                if p_c < alpha
                    y_top = max(y_all(strcmp(gb_all, bundle_display{b}))) + 0.02;
                    plot(ax_tile, [b-0.3 b+0.3], [y_top y_top], 'k-', 'LineWidth', 0.8);
                    if p_c < 0.001; txt = '***';
                    elseif p_c < 0.01; txt = '**';
                    else; txt = '*';
                    end
                    text(ax_tile, b, y_top + 0.01, txt, ...
                        'HorizontalAlignment', 'center', 'FontSize', 9);
                end
            end

            ylim(ax_tile, [max(0, min(y_all)-0.05) 1]);
            title(ax_tile, sprintf('%s | %s', upper(method), ori), 'FontSize', 10);
            ylabel(ax_tile, 'r²', 'FontSize', 9);
            if o == 1 && m == 1
                legend(ax_tile, 'Location', 'southwest', 'FontSize', 8);
            else
                legend(ax_tile, 'off');
            end
            box(ax_tile, 'on');
        end
    end

    fname = sprintf('source_vs_sensor_%s', lower(mod_label));
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

%% Figure set 2: MSG vs ESG per perturbation type

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
    title(tl, sprintf('MSG vs ESG — %s r²  (rank-sum, Bonferroni) | %s', pert, pert_dsp), ...
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

            for b = 1:numel(bundle_display)
                rr = res_sub(strcmp(res_sub.method, method) & ...
                             res_sub.bundle == b & ...
                             strcmp(res_sub.orientation, ori), :);
                if isempty(rr); continue; end
                p_c = min(rr.p_corrected);
                if p_c < alpha
                    y_top = max(y_all(strcmp(gb_all, bundle_display{b}))) + 0.02;
                    plot(ax_tile, [b-0.3 b+0.3], [y_top y_top], 'k-', 'LineWidth', 0.8);
                    if p_c < 0.001; txt = '***';
                    elseif p_c < 0.01; txt = '**';
                    else; txt = '*';
                    end
                    text(ax_tile, b, y_top + 0.01, txt, ...
                        'HorizontalAlignment', 'center', 'FontSize', 9);
                end
            end

            ylim(ax_tile, [max(0, min(y_all)-0.05) 1]);
            title(ax_tile, sprintf('%s | %s', upper(method), ori), 'FontSize', 10);
            ylabel(ax_tile, 'r²', 'FontSize', 9);
            if o == 1 && m == 1
                legend(ax_tile, 'Location', 'southwest', 'FontSize', 8);
            else
                legend(ax_tile, 'off');
            end
            box(ax_tile, 'on');
        end
    end

    fname = sprintf('msg_vs_esg_%s', pert);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
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
    if isfile(f)
        dat.source = load(f);
        dat.have_source = true;
    end

    f = fullfile(base_path, 'pert_sensor_rsq.mat');
    if isfile(f)
        dat.sensor = load(f);
        dat.have_sensor = true;
    end

    f = fullfile(base_path, 'pert_cond_rsq.mat');
    if isfile(f)
        dat.cond = load(f);
        dat.have_cond = true;
    end
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

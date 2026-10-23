% pt_summary_table - Stage 6: every factor on one scale
%
% Pulls the long-format results of stages 1, 2, 3 and 5 into one ranking per
% modality: the forward-model-type differences (Biot-Savart and FEM against
% the MRI-derived BEM), every perturbation bundle for every model, the FEM
% re-meshing floor, sensor noise alone, and perturbation plus noise at the
% reference noise level. All rows use the same four metrics, so the table answers "which
% source of error matters most" directly — the msg_pert counterpart of
% msg_fwd's compute_hierarchy_table. Nothing is recomputed, so this table
% can never disagree with the stage outputs.
%
% USAGE:
%   pt_summary_table            % run last; missing stages are skipped
%
% INPUT (from <stage_results_dir>):
%   1_baseline/stage1_long.csv, 2_within_msg/stage2_long.csv,
%   3_within_esg/stage3_long.csv, 5_noise/stage5_long.csv
%
% OUTPUT (to <stage_results_dir>/6_summary/):
%   summary_all.csv                     every row of every stage, one file
%   summary_<modality>.csv/.txt/.tex    ranking at the headline orientation,
%                                       whole array, sorted by RE; every row a
%                                       cord median, as in stages 1-4
%   summary_<modality>_axes.csv         the same for every axis set
%   summary_<modality>_c7.*             ranking at the evoked C7 source, from
%                                       stage 5: each perturbation noise-free,
%                                       noise alone, and perturbation + noise
%   summary_forest.png/.fig             RE with CI per factor, both modalities
%   summary_forest_c7.png/.fig          the same at the C7 source
%
% READING IT
%   The unit behind each CI differs by row and is stated in the unit column:
%   'source' rows (model type) resample cord positions, 'realisation' rows
%   resample perturbation draws, 'noise realisation' rows resample noise
%   draws. Compare magnitudes across rows; do not treat overlapping CIs from
%   different units as a test.
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

fprintf('Stage 6 — summary across stages\n\n');

out_dir = fullfile(stage_results_dir, '6_summary');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

files = {fullfile('1_baseline', 'stage1_long.csv'), ...
         fullfile('2_within_msg', 'stage2_long.csv'), ...
         fullfile('3_within_esg', 'stage3_long.csv'), ...
         fullfile('5_noise', 'stage5_long.csv')};
parts = {};
for f = files
    p = fullfile(stage_results_dir, f{1});
    if isfile(p)
        T = readtable(p, 'TextType', 'char', 'Delimiter', ',');
        parts{end+1} = T; %#ok<SAGROW>
        fprintf('  read %s (%d rows)\n', f{1}, height(T));
    else
        fprintf('  missing %s — skipped\n', f{1});
    end
end
if isempty(parts), error('No stage outputs found under %s.', stage_results_dir); end
L = vertcat(parts{:});
for v = {'stage', 'modality', 'system', 'method', 'factor', 'level', ...
         'orientation', 'axis_set', 'site', 'unit'}
    L.(v{1}) = cellstr(string(L.(v{1})));
end
writetable(L, fullfile(out_dir, 'summary_all.csv'));

ref_tag  = sprintf('@ %gx', noise_reference_factor);
zero_tag = '@ 0x';
meth     = modality_compare_method;
is5      = strcmp(L.stage, '5_noise');

% ---- Cord ranking: model type (stage 1), every perturbation bundle (stages
% 2/3, every method), noise alone and perturbation + noise at the reference
% level, all as cord medians
keep = strcmp(L.factor, 'model_type') | startsWith(L.stage, {'2_', '3_'}) | ...
       (is5 & strcmp(L.site, 'cord') & endsWith(L.system, ref_tag));
R = label_rows(L(keep, :), ref_tag, pert_type_display);

% ---- C7 ranking: stage 5 at the evoked source. Noise-free rows (0x) give
% the perturbation alone at C7; they are identical for systems that share a
% forward model, so one is kept per modality, method, type and bundle.
keep7 = is5 & strcmp(L.site, 'C7') & ...
        (endsWith(L.system, ref_tag) | ...
         (endsWith(L.system, zero_tag) & ~strcmp(L.factor, 'noise_only')));
R7 = L(keep7, :);
if height(R7) > 0
    nf = endsWith(R7.system, zero_tag);
    key = strcat(R7.modality, '|', R7.method, '|', R7.factor, '|', R7.level, '|', ...
                 R7.orientation, '|', R7.axis_set);
    key(~nf) = strcat(key(~nf), '|', R7.system(~nf));     % noisy rows all kept
    [~, first] = unique(key, 'stable');
    R7 = R7(sort(first), :);
end
R7 = label_rows(R7, ref_tag, pert_type_display);

write_ranking(R, '', 'cord median', out_dir, headline_orientation, ori_titles, ...
    noise_reference_factor);
write_ranking(R7, '_c7', sprintf('%s source', evoked_src_name), out_dir, ...
    headline_orientation, ori_titles, noise_reference_factor);

plot_forest(R, meth, '', 'cord median', out_dir, headline_orientation, ori_titles, ...
    noise_reference_factor, pert_type_colors, fig_resolution);
plot_forest(R7, meth, '_c7', sprintf('%s source', evoked_src_name), out_dir, ...
    headline_orientation, ori_titles, noise_reference_factor, pert_type_colors, ...
    fig_resolution);

fprintf('\nStage 6 complete: %s\n', out_dir);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function R = label_rows(R, ref_tag, tdisp)
% Family and row label for every kept row
    R.row_label = cell(height(R), 1);
    R.family    = cell(height(R), 1);
    for i = 1:height(R)
        sys = strtrim(regexprep(R.system{i}, '@ \S+$', ''));
        switch R.factor{i}
            case 'model_type'
                R.family{i}    = 'Forward model type';
                R.row_label{i} = R.level{i};
            case 'noise_only'
                R.family{i}    = 'Sensor noise';
                R.row_label{i} = sprintf('%s noise alone', sys);
            case 'remesh'
                R.family{i}    = 'FEM re-meshing';
                R.row_label{i} = sprintf('FEM re-meshing (%s)', R.level{i});
            otherwise
                fam = tdisp.(R.factor{i});
                if ~strcmp(R.stage{i}, '5_noise')
                    R.family{i}    = sprintf('Perturbation (%s)', pt_method_label(R.method{i}));
                    R.row_label{i} = sprintf('%s, %s', fam, R.level{i});
                elseif endsWith(R.system{i}, ref_tag)
                    R.family{i}    = 'Perturbation + noise';
                    R.row_label{i} = sprintf('%s, %s + %s noise', fam, R.level{i}, sys);
                else
                    R.family{i}    = sprintf('Perturbation (%s)', pt_method_label(R.method{i}));
                    R.row_label{i} = sprintf('%s, %s (%s)', fam, R.level{i}, upper(R.modality{i}));
                end
        end
    end
end

function write_ranking(R, suffix, site_txt, out_dir, ori_head, ori_titles, ref_factor)
    mods = unique(R.modality, 'stable');
    for m = 1:numel(mods)
        md = mods{m};
        Rm = R(strcmp(R.modality, md), :);
        writetable(Rm, fullfile(out_dir, sprintf('summary_%s%s_axes.csv', md, suffix)));

        H = Rm(strcmp(Rm.orientation, ori_head) & strcmp(Rm.axis_set, 'Array'), :);
        H = sortrows(H, 're_median', 'descend');
        H.re = H.re_median;   % pt_write_table's 'ci' format reads <var>, <var>_ci_lo, <var>_ci_hi
        pt_write_table(H, fullfile(out_dir, sprintf('summary_%s%s', md, suffix)), struct( ...
            'cols', {{ 'row_label', 'Source of error', 'Source of error', '%s'; ...
                       'family', 'Family', 'Family', '%s'; ...
                       're', 'RE % [95% CI]', 'RE (\%) [95\% CI]', 'ci'; ...
                       're_p95', 'RE p95 [95% CI]', 'RE P95 (\%) [95\% CI]', 'ci'; ...
                       'rsq_median', 'r2', '$r^2$', '%.4f'; ...
                       'rdm_median', 'RDM', 'RDM', '%.4f'; ...
                       'gain_median', 'gain %', 'Gain (\%)', '%+.2f'; ...
                       'absgain_median', '|gain| %', '$|$Gain$|$ (\%)', '%.2f'; ...
                       'unit', 'CI over', 'CI over', '%s' }}, ...
            'title', sprintf('STAGE 6 — EVERY SOURCE OF ERROR ON ONE SCALE: %s, %s (%s, whole array)', ...
                upper(md), site_txt, ori_titles.(ori_head)), ...
            'notes', {{'Sorted by median RE. Noise rows are at the reference noise level.', ...
                       'CI resampling unit differs by row (column "CI over"): compare magnitudes, not CI overlap.'}}, ...
            'caption', sprintf(['All sources of forward-field error for %s on one scale ' ...
                       '(%s): each perturbation bundle, sensor noise alone and ' ...
                       'perturbation plus noise at %g$\\times$ the noise floor. Sorted ' ...
                       'by median relative error.'], upper(md), site_txt, ref_factor), ...
            'label', sprintf('tab:summary-%s%s', md, strrep(suffix, '_', '-'))));
    end
end

function plot_forest(R, meth, suffix, site_txt, out_dir, ori_head, ori_titles, ...
        ref_factor, type_cols, res)
% RE with CI per factor, one panel per modality, headline orientation, whole
% array, compare method (+ noise rows)
    mods = unique(R.modality, 'stable');
    if isempty(mods), return; end
    fig = figure('Color', 'w', 'Position', [40 40 700*numel(mods), 900]);
    tl  = tiledlayout(1, numel(mods), 'TileSpacing', 'compact');
    fam_cols = containers.Map( ...
        {'Forward model type', 'Sensor noise', 'Perturbation + noise', 'FEM re-meshing'}, ...
        {[0.3 0.3 0.3], [0.6 0.6 0.6], [0.55 0.25 0.60], [0.75 0.75 0.75]});
    for m = 1:numel(mods)
        H = R(strcmp(R.modality, mods{m}) & strcmp(R.orientation, ori_head) & ...
              strcmp(R.axis_set, 'Array') & ...
              (~startsWith(R.family, 'Perturbation (') | strcmp(R.method, meth)), :);
        H = sortrows(H, 're_median', 'ascend');
        ax = nexttile(tl); hold(ax, 'on');
        for i = 1:height(H)
            if isKey(fam_cols, H.family{i})
                col = fam_cols(H.family{i});
            else
                col = type_cols.(H.factor{i});
            end
            lo = H.re_ci_lo(i); hi = H.re_ci_hi(i);
            if ~isnan(lo), plot(ax, [lo hi], [i i], '-', 'Color', col, 'LineWidth', 2); end
            plot(ax, H.re_median(i), i, 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k');
        end
        set(ax, 'YTick', 1:height(H), 'YTickLabel', H.row_label, 'XScale', 'log', ...
            'TickLabelInterpreter', 'none', 'FontSize', 8, 'TickDir', 'out');
        grid(ax, 'on');
        xlabel(ax, 'RE (%), log scale');
        title(ax, upper(mods{m}));
    end
    title(tl, sprintf('Every source of error on one scale — %s (%s, whole array; noise at %gx floor)', ...
        site_txt, ori_titles.(ori_head), ref_factor), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['summary_forest' suffix], res);
end

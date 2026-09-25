% run_perturbation_analysis - Master script for the msg_pert analysis (Phase 2)
%
% Runs every analysis stage in order, for both modalities, from the
% leadfields computed in msg_fwd. See RUN_ORDER.md for what each stage needs
% and costs, and INTERPRETATION.md for how to read the outputs.
%
% STAGES
%   0  Load and compute     per modality: pt_load_leadfields, pt_compute_metrics
%   1  Baseline             pt_baseline            noise-free fields: topoplots,
%                                                  amplitude, model-type metrics
%   2  Within MSG           pt_within_modality     perturbation types and sizes
%   3  Within ESG           pt_within_modality     the same for ESG
%   4  MSG vs ESG           pt_compare_modalities  the two modalities compared,
%                                                  with statistics
%   5  Sensor noise         pt_noise_simulate,     stages 2-4 repeated through
%                           pt_noise_analyse       realistic sensor noise
%   6  Summary              pt_summary_table       every factor on one scale
%
%   Optional, off by default: the earlier r2-only figures and tables
%   (pt_compute_rsq, pt_plot_curves, pt_plot_heatmaps, pt_plot_displacement,
%   pt_plot_slope_vs_position, pt_compute_table, pt_compare_perturbations).
%
% USAGE:
%   run_perturbation_analysis
%
%   Set the flags below to run a subset. Each stage reads only files written
%   by earlier stages, so any stage can be re-run on its own once those exist:
%     pt_modality('set', 'esg');  pt_within_modality
%
% CONFIGURATION:
%   config_pert.m — paths, methods, perturbation definitions, statistics and
%   noise settings. Nothing here needs editing apart from the flags.
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

% ---- Stage flags --------------------------------------------------------
run_stage0 = true;    % load leadfields + compute metrics (per modality)
run_stage1 = true;    % noise-free baseline
run_stage2 = true;    % within MSG
run_stage3 = true;    % within ESG
run_stage4 = true;    % MSG vs ESG
run_stage5 = true;    % sensor noise (simulation is the slow part)
run_stage6 = true;    % summary across stages
run_legacy = false;   % earlier r2-only figures and tables, per modality

fprintf('  MSG Perturbation Analysis Pipeline\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n\n');

pt_modality('clear');
config_pert;
pt_add_functions;
mods_to_run = pert_modalities;
t0 = tic;
status = {};

% =========================================================================
% STAGE 0 — load leadfields and compute all metrics, per modality
% =========================================================================
if run_stage0
    for mi = 1:numel(mods_to_run)
        pt_modality('set', mods_to_run{mi});
        banner(sprintf('STAGE 0 — %s: load and compute metrics', upper(mods_to_run{mi})));
        ok = pt_run_step('pt_load_leadfields.m');
        if ok, ok = pt_run_step('pt_compute_metrics.m'); end
        status(end+1, :) = {sprintf('0 %s', mods_to_run{mi}), ok}; %#ok<SAGROW>
    end
end

% =========================================================================
% STAGE 1 — baseline
% =========================================================================
if run_stage1
    pt_modality('clear');
    banner('STAGE 1 — noise-free baseline');
    status(end+1, :) = {'1 baseline', pt_run_step('pt_baseline.m')};
end

% =========================================================================
% STAGES 2 and 3 — within each modality
% =========================================================================
stage_flag = struct('msg', run_stage2, 'esg', run_stage3);
stage_no   = struct('msg', 2, 'esg', 3);
for mi = 1:numel(mods_to_run)
    md = mods_to_run{mi};
    if ~stage_flag.(md), continue; end
    pt_modality('set', md);
    banner(sprintf('STAGE %d — within %s', stage_no.(md), upper(md)));
    status(end+1, :) = {sprintf('%d within %s', stage_no.(md), md), ...
        pt_run_step('pt_within_modality.m')}; %#ok<SAGROW>
end

% =========================================================================
% STAGE 4 — MSG vs ESG
% =========================================================================
if run_stage4 && all(ismember({'msg', 'esg'}, mods_to_run))
    pt_modality('clear');
    banner('STAGE 4 — MSG vs ESG');
    status(end+1, :) = {'4 msg vs esg', pt_run_step('pt_compare_modalities.m')};
end

% =========================================================================
% STAGE 5 — sensor noise
% =========================================================================
if run_stage5
    pt_modality('clear');
    banner('STAGE 5 — sensor noise');
    ok = pt_run_step('pt_noise_simulate.m');
    if ok, ok = pt_run_step('pt_noise_analyse.m'); end
    status(end+1, :) = {'5 noise', ok};
end

% =========================================================================
% STAGE 6 — summary
% =========================================================================
if run_stage6
    pt_modality('clear');
    banner('STAGE 6 — summary across stages');
    status(end+1, :) = {'6 summary', pt_run_step('pt_summary_table.m')};
end

% =========================================================================
% OPTIONAL — earlier r2-only outputs
% =========================================================================
if run_legacy
    legacy = {'pt_compute_rsq.m', 'pt_plot_curves.m', 'pt_plot_heatmaps.m', ...
              'pt_plot_displacement.m', 'pt_plot_slope_vs_position.m', ...
              'pt_compute_table.m'};
    for mi = 1:numel(mods_to_run)
        pt_modality('set', mods_to_run{mi});
        banner(sprintf('LEGACY r2 outputs — %s', upper(mods_to_run{mi})));
        for k = 1:numel(legacy)
            pt_run_step(legacy{k});
        end
    end
    pt_modality('clear');
    pt_run_step('pt_compare_perturbations.m');
end

pt_modality('clear');
fprintf('\n=========================================\n');
fprintf('  Pipeline finished in %.1f min\n', toc(t0) / 60);
for k = 1:size(status, 1)
    if status{k, 2}, s = 'ok'; else, s = 'FAILED'; end
    fprintf('  stage %-16s %s\n', status{k, 1}, s);
end
fprintf('  Results: %s\n', stage_results_dir);
fprintf('=========================================\n');


function banner(txt)
    fprintf('\n==================================================\n');
    fprintf('  %s\n', txt);
    fprintf('==================================================\n\n');
end

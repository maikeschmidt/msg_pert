function pt_run_one_modality()
% pt_run_one_modality - Run Phase 2 analysis steps for the ACTIVE modality
%
% Runs the per-modality perturbation analysis (load → r² → curves → heatmaps →
% displacement → slope → tables) for whichever modality pt_modality currently
% selects. run_perturbation_analysis calls this once per modality.
%
% Why a function, not inlined in the master loop: every sub-script starts with
% `clearvars`, and run() executes them in the caller's workspace — so a plain
% loop in the master would have its loop counter wiped by the first sub-script.
% Wrapping the step sequence in a function isolates those clearvars to THIS
% function's workspace, leaving the master's loop intact.
%
% The active modality's paths (forward_fields_base, save_base_dir, ...) come from
% config_pert, which reads pt_modality. Each sub-script re-runs config_pert, so
% the config variables are refreshed in this workspace after every run() call.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

config_pert;   % establishes forward_fields_base / save_base_dir for this modality
fprintf('  Results → %s\n\n', save_base_dir);

% -------------------------------------------------------------------------
% STEP 1: Load and organise perturbed leadfields
% -------------------------------------------------------------------------
fprintf('[1/7] Loading and organising perturbed leadfields...\n');
try
    run('pt_load_leadfields.m');
    fprintf('[1/7] Complete.\n\n');
catch err
    fprintf('ERROR: pt_load_leadfields failed:\n  %s\n', err.message);
    fprintf('Skipping the rest of this modality.\n');
    return;
end

source_rsq_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
sensor_rsq_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
cond_rsq_file   = fullfile(forward_fields_base, 'pert_cond_rsq.mat');

% -------------------------------------------------------------------------
% STEP 2: Compute r²
% -------------------------------------------------------------------------
fprintf('[2/7] Computing perturbation r² (source, sensor, conductivity)...\n');
try
    run('pt_compute_rsq.m');
    fprintf('[2/7] Complete.\n\n');
catch err
    fprintf('ERROR: pt_compute_rsq failed:\n  %s\n', err.message);
    return;
end

have_source_rsq = isfile(source_rsq_file);
have_sensor_rsq = isfile(sensor_rsq_file);
have_cond_rsq   = isfile(cond_rsq_file);
have_any_rsq    = have_source_rsq || have_sensor_rsq || have_cond_rsq;

% -------------------------------------------------------------------------
% STEP 3: r² curve figures
% -------------------------------------------------------------------------
if have_any_rsq
    fprintf('[3/7] Plotting perturbation curves...\n');
    try
        run('pt_plot_curves.m');
        fprintf('[3/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_curves failed:\n  %s\nContinuing...\n\n', err.message);
    end
else
    fprintf('[3/7] Skipping pt_plot_curves — no r² files found.\n\n');
end

% -------------------------------------------------------------------------
% STEP 4: Heatmap summaries
% -------------------------------------------------------------------------
if have_source_rsq || have_sensor_rsq || have_cond_rsq
    fprintf('[4/7] Plotting heatmaps...\n');
    try
        run('pt_plot_heatmaps.m');
        fprintf('[4/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_heatmaps failed:\n  %s\nContinuing...\n\n', err.message);
    end
else
    fprintf('[4/7] Skipping pt_plot_heatmaps — no r² files found.\n\n');
end

% -------------------------------------------------------------------------
% STEP 5: Displacement / perturbation vs r²
% -------------------------------------------------------------------------
if have_any_rsq
    fprintf('[5/7] Plotting displacement vs r²...\n');
    fprintf('       (individual figures: cervical region; combined + table: full cord)\n');
    try
        run('pt_plot_displacement.m');
        fprintf('[5/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_displacement failed:\n  %s\nContinuing...\n\n', err.message);
    end
else
    fprintf('[5/7] Skipping pt_plot_displacement — no r² files found.\n\n');
end

% -------------------------------------------------------------------------
% STEP 6: Slope of r² vs cord position
% -------------------------------------------------------------------------
sensor_tbl = fullfile(save_base_dir, 'perturbation_analysis', 'sensor', ...
    'sensor_disp_trend_table.tsv');
source_tbl = fullfile(save_base_dir, 'perturbation_analysis', 'source', ...
    'source_disp_trend_table.tsv');
cond_tbl   = fullfile(save_base_dir, 'perturbation_analysis', 'cond', ...
    'cond_disp_trend_table.tsv');
have_any_tbl = isfile(sensor_tbl) || isfile(source_tbl) || isfile(cond_tbl);

if have_any_tbl
    fprintf('[6/7] Plotting slope vs cord position (full cord)...\n');
    try
        run('pt_plot_slope_vs_position.m');
        fprintf('[6/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_slope_vs_position failed:\n  %s\nContinuing...\n\n', err.message);
    end
else
    fprintf('[6/7] Skipping pt_plot_slope_vs_position — no trend tables found.\n\n');
end

% -------------------------------------------------------------------------
% STEP 7: Summary tables
% -------------------------------------------------------------------------
if have_any_rsq
    fprintf('[7/7] Computing summary tables...\n');
    try
        run('pt_compute_table.m');
        fprintf('[7/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_compute_table failed:\n  %s\nContinuing...\n\n', err.message);
    end
else
    fprintf('[7/7] Skipping pt_compute_table — no r² files found.\n\n');
end
end

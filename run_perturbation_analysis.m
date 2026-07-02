% run_perturbation_analysis - Master script for MSG perturbation analysis pipeline
%
% Runs the post-forward-model analysis pipeline for msg_pert. All leadfields
% for shifted geometries must have been computed in msg_fwd before running
% this script.
%
% This script runs the analysis in the correct order. Each step loads
% config_pert and its required .mat files independently, so individual
% steps can also be run standalone.
%
% WORKFLOW (two-phase):
%
%   Phase 1 — Perturbation generation (run ONCE per study, in msg_fwd):
%     1. pt_generate_source_shifts         24 source-shift geometry files
%                                          (3 bundles x 8 random shifts)
%     2. pt_generate_sensor_shifts         24 sensor-shift geometry files
%                                          (3 bundles x 8 random shifts)
%     3. run_bem_leadfields                BEM leadfields for source + sensor
%                                          shift geometries
%     4. run_conductivity_perturbation     BEM leadfields with perturbed tissue
%                                          conductivities (3 bundles x 8 shifts)
%     5. [Run other forward models as needed: BS law, sphere, FEM]
%
%   Phase 2 — Analysis (this script, 7 steps):
%     1. pt_load_leadfields    Load + organise all leadfields into one struct
%     2. pt_compute_rsq        Compute r² for source, sensor, and cond shifts
%     3. pt_plot_curves        r² vs cord distance (detail, summary, cross-model)
%     4. pt_plot_heatmaps      Pairwise RE and r² heatmaps
%     5. pt_plot_displacement  Displacement / perturbation vs r² figures
%                              (individual figures: cervical region only;
%                               combined + trend table: full cord)
%     6. pt_plot_slope_vs_position  Slope of r² change vs cord position
%     7. pt_compute_table      Summary tables (.txt and .csv)
%
% USAGE:
%   run_perturbation_analysis
%
% CONFIGURATION:
%   - Update config_pert.m paths and parameters before running
%   - Configure which methods are available (BEM/FEM/BS/sphere/BEM-cond)
%     in pt_load_leadfields.m
%   - Set n_cond_compartments in config_pert.m to match your BEM geometry
%
% NOTES:
%   - Steps 3-7 auto-skip if their required .mat file does not exist
%   - pt_plot_slope_vs_position requires pt_plot_displacement to have run first
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
% Date:   June 2026
%
% This file is part of the MSG Perturbation Toolbox (msg_pert).
% Used in conjunction with msg_coreg and msg_fwd:
%   https://github.com/maikeschmidt/msg_coreg
%   https://github.com/maikeschmidt/msg_fwd

clearvars
close all
clc

fprintf('  MSG Perturbation Analysis Pipeline\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n\n');

config_pert;

% =========================================================================
% STEP 1: Load and organise perturbed leadfields
% =========================================================================

fprintf('[1/7] Loading and organising perturbed leadfields...\n');
try
    run('pt_load_leadfields.m');
    fprintf('[1/7] Complete.\n\n');
catch err
    fprintf('ERROR: pt_load_leadfields failed:\n  %s\n', err.message);
    fprintf('Check that forward model output paths and method flags are\n');
    fprintf('configured in pt_load_leadfields.m before running.\n');
    return;
end

source_rsq_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
sensor_rsq_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
cond_rsq_file   = fullfile(forward_fields_base, 'pert_cond_rsq.mat');

% =========================================================================
% STEP 2: Compute r²
% =========================================================================

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

% =========================================================================
% STEP 3: r² curve figures
% =========================================================================

if have_any_rsq
    fprintf('[3/7] Plotting perturbation curves...\n');
    try
        run('pt_plot_curves.m');
        fprintf('[3/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_curves failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[3/7] Skipping pt_plot_curves — no r² files found.\n\n');
end

% =========================================================================
% STEP 4: Heatmap summaries
% =========================================================================

if have_source_rsq || have_sensor_rsq
    fprintf('[4/7] Plotting heatmaps...\n');
    try
        run('pt_plot_heatmaps.m');
        fprintf('[4/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_heatmaps failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[4/7] Skipping pt_plot_heatmaps — no source/sensor r² files found.\n\n');
end

% =========================================================================
% STEP 5: Displacement / perturbation vs r²
% =========================================================================

if have_any_rsq
    fprintf('[5/7] Plotting displacement vs r²...\n');
    fprintf('       (individual figures: cervical region; combined + table: full cord)\n');
    try
        run('pt_plot_displacement.m');
        fprintf('[5/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_displacement failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[5/7] Skipping pt_plot_displacement — no r² files found.\n\n');
end

% =========================================================================
% STEP 6: Slope of r² vs cord position
% =========================================================================

% Check that at least one trend table exists (produced by pt_plot_displacement)
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
        fprintf('WARNING: pt_plot_slope_vs_position failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[6/7] Skipping pt_plot_slope_vs_position — no trend tables found.\n');
    fprintf('       Run pt_plot_displacement first.\n\n');
end

% =========================================================================
% STEP 7: Summary tables
% =========================================================================

if have_any_rsq
    fprintf('[7/7] Computing summary tables...\n');
    try
        run('pt_compute_table.m');
        fprintf('[7/7] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_compute_table failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[7/7] Skipping pt_compute_table — no r² files found.\n\n');
end

fprintf('=========================================\n');
fprintf('  Perturbation analysis pipeline complete.\n');
fprintf('  Figures saved to: %s\n', fullfile(save_base_dir, 'perturbation_analysis'));
fprintf('=========================================\n');

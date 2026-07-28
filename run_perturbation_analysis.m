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

config_pert;                 % defines pert_modalities + combined_results_dir
mods_to_run = pert_modalities;

% =========================================================================
% PER-MODALITY ANALYSIS (loops MSG and ESG — no manual switching)
% =========================================================================
% pt_run_one_modality is a FUNCTION so its clearvars-heavy sub-scripts cannot
% wipe this loop. Each modality writes to its own forward_fields_base and
% save_base_dir (set per modality in config_pert).

for mi = 1:numel(mods_to_run)
    md = mods_to_run{mi};
    pt_modality('set', md);
    fprintf('\n==================================================\n');
    fprintf('  MODALITY %d/%d: %s\n', mi, numel(mods_to_run), upper(md));
    fprintf('==================================================\n\n');
    pt_run_one_modality();
end

% =========================================================================
% COMBINED: cross-perturbation and MSG-vs-ESG comparison
% =========================================================================
% Reads each modality's pert_*_rsq.mat (from its forward_fields_base) and writes
% comparison figures to combined_results_dir. Both are configured in config_pert.

pt_modality('clear');
fprintf('\n==================================================\n');
fprintf('  COMBINED: cross-perturbation / MSG vs ESG\n');
fprintf('==================================================\n\n');
pt_run_compare();   % function file: isolates its clearvars, callable interactively too

fprintf('=========================================\n');
fprintf('  Perturbation analysis pipeline complete.\n');
fprintf('  Ran modalities: %s\n', strjoin(mods_to_run, ', '));
fprintf('=========================================\n');

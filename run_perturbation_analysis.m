% run_perturbation_analysis - Master script for MSG perturbation analysis pipeline
%
% Runs the post-forward-model analysis pipeline for msg_pert. All leadfields
% for shifted geometries must have been computed in msg_fwd before running
% this script.
%
% This script runs the analysis in the correct order. Each analysis script
% loads config_pert and the relevant .mat files independently, so individual
% steps can also be run standalone.
%
% WORKFLOW (two-phase):
%
%   Phase 1 — Geometry generation (run ONCE per study):
%     1. pt_generate_source_shifts   generate 18 source-shift geometry files
%     2. pt_generate_sensor_shifts   generate 24 sensor-shift geometry files
%     3. [Copy printed filenames into msg_fwd and run forward models]
%
%   Phase 2 — Analysis (this script):
%     1. load_and_organise_leadfields  (run in msg_fwd to produce organised .mat)
%     2. pt_compute_rsq                compute r² for source and sensor shifts
%     3. pt_plot_curves                r² vs cord distance figures
%     4. pt_plot_displacement          displacement vs r² figures (sensor mode)
%     5. pt_compute_table              summary tables (.txt and .csv)
%
% USAGE:
%   run_perturbation_analysis
%
% NOTES:
%   - Run load_and_organise_leadfields.m in msg_fwd first to produce
%     leadfields_organised.mat from the shifted leadfield files
%   - All paths are configured in config_pert.m — update that file first
%   - Steps 3-5 are auto-skipped if their required .mat file does not yet exist
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

source_rsq_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
sensor_rsq_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');

% STEP 1: Compute r²
% Loads leadfields_organised.mat — ensure this exists (run load_and_organise_leadfields
% in msg_fwd on the perturbed leadfield files first).

fprintf('[1/4] Computing perturbation r²...\n');
try
    run('pt_compute_rsq.m');
    fprintf('[1/4] Complete.\n\n');
catch err
    fprintf('ERROR: pt_compute_rsq failed:\n  %s\n', err.message);
    fprintf('Ensure leadfields_organised.mat exists in forward_fields_base.\n');
    fprintf('Run load_and_organise_leadfields.m in msg_fwd first.\n');
    return;
end

have_source_rsq = isfile(source_rsq_file);
have_sensor_rsq = isfile(sensor_rsq_file);

% STEP 2: Sensitivity curve figures

if have_source_rsq || have_sensor_rsq
    fprintf('[2/4] Plotting perturbation curves...\n');
    try
        run('pt_plot_curves.m');
        fprintf('[2/4] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_curves failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[2/4] Skipping pt_plot_curves — no r² files found.\n\n');
end

% STEP 3: Displacement vs r² (sensor mode only)

if have_sensor_rsq
    fprintf('[3/4] Plotting displacement vs r²...\n');
    try
        run('pt_plot_displacement.m');
        fprintf('[3/4] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_plot_displacement failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[3/4] Skipping pt_plot_displacement — no sensor r² file found.\n\n');
end

% STEP 4: Summary tables

if have_source_rsq || have_sensor_rsq
    fprintf('[4/4] Computing summary tables...\n');
    try
        run('pt_compute_table.m');
        fprintf('[4/4] Complete.\n\n');
    catch err
        fprintf('WARNING: pt_compute_table failed:\n  %s\n', err.message);
        fprintf('Continuing...\n\n');
    end
else
    fprintf('[4/4] Skipping pt_compute_table — no r² files found.\n\n');
end

fprintf('  Perturbation analysis pipeline complete.\n');
fprintf('  Figures saved to: %s\n', save_base_dir);

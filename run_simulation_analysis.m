% run_simulation_analysis - Master script for the realistic-measurement pipeline
%
% The pt_* perturbation pipeline answers: how wrong does the forward field go
% when the MODEL is wrong? This pipeline answers a complementary question: given
% a PERFECT forward model, how much does real SENSOR NOISE degrade what you can
% actually measure, and does that differ between SQUID MSG, OP-MSG, and ESG?
%
% Motivation: the perturbation results showed MSG (not ESG) was generally the
% more sensitive modality, with almost no MSG/ESG difference for conductivity.
% The candidate explanation is that BEM MSG produces highly individualised,
% spatially sharp fields, whereas Biot-Savart MSG and BEM ESG produce smoother,
% more diffuse fields that tolerate error better. Step 1 tests that directly by
% putting the three field maps side by side.
%
% WORKFLOW (4 steps):
%   1. sim_plot_topoplots       Perfect forward fields for Biot-Savart MSG,
%                               BEM MSG, and BEM ESG at one source, for both
%                               arrays and every sensor axis and orientation.
%   2. sim_simulate_noise       Simulate a 1 nA*m 90 Hz source at every point on
%                               the cord, add sensor noise for each system across
%                               a sweep of noise levels, and score the result
%                               against the noise-free field with r-squared.
%   3. sim_plot_noise_curves    r^2 vs noise level per system — averaged across
%                               the cord, at one chosen source, and as a function
%                               of cord position.
%   4. sim_plot_noise_topoplot  What one system actually measures at a chosen
%                               source and noise level, next to the perfect field.
%
% USAGE:
%   run_simulation_analysis
%
% CONFIGURATION:
%   All paths, models, systems, noise floors, and the signal definition live in
%   config_sim.m. Set those before running.
%
% A NOTE ON COMPARING MSG WITH ESG:
%   MSG leadfields are fT/nAm and ESG leadfields are uV/nAm. Absolute noise
%   levels are NOT comparable across the two. Everything here is expressed in
%   multiples of each system's own baseline noise floor, which IS comparable:
%   it asks how far each system sits from its own real-world operating point.
%
%   The ESG noise floor (~1 uV/sqrt(Hz)) is an amplifier-noise estimate. In a
%   real recording, cardiac artefact (~29 uV cervical, ~657 uV lumbar) dwarfs
%   it. The ESG curves here are therefore an OPTIMISTIC bound on ESG, not a
%   realistic operating point — worth stating explicitly in any write-up.
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

fprintf('  MSG Realistic-Measurement Simulation Pipeline\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n\n');

config_sim;

rsq_file = fullfile(sim_out_dir, 'sim_noise_rsq.mat');


% =========================================================================
% STEP 1: Perfect forward field topoplots
% =========================================================================

fprintf('[1/4] Plotting perfect forward fields (Biot-Savart MSG, BEM MSG, BEM ESG)...\n');
try
    run('sim_plot_topoplots.m');
    fprintf('[1/4] Complete.\n\n');
catch err
    fprintf('WARNING: sim_plot_topoplots failed:\n  %s\n', err.message);
    fprintf('Check the model paths in config_sim.m. Continuing...\n\n');
end


% =========================================================================
% STEP 2: Simulate data + sensor noise, compute r^2
% =========================================================================

fprintf('[2/4] Simulating source data with sensor noise...\n');
try
    run('sim_simulate_noise.m');
    fprintf('[2/4] Complete.\n\n');
catch err
    fprintf('ERROR: sim_simulate_noise failed:\n  %s\n', err.message);
    fprintf('Steps 3-4 depend on this — stopping.\n');
    return;
end

if ~isfile(rsq_file)
    fprintf('ERROR: %s was not produced. Stopping.\n', rsq_file);
    return;
end


% =========================================================================
% STEP 3: r^2 vs noise level curves
% =========================================================================

fprintf('[3/4] Plotting r^2 vs noise level...\n');
try
    run('sim_plot_noise_curves.m');
    fprintf('[3/4] Complete.\n\n');
catch err
    fprintf('WARNING: sim_plot_noise_curves failed:\n  %s\n', err.message);
    fprintf('Continuing...\n\n');
end


% =========================================================================
% STEP 4: Topoplot at a chosen source and noise level
% =========================================================================

fprintf('[4/4] Plotting measured (noisy) topoplots...\n');
try
    run('sim_plot_noise_topoplot.m');
    fprintf('[4/4] Complete.\n\n');
catch err
    fprintf('WARNING: sim_plot_noise_topoplot failed:\n  %s\n', err.message);
    fprintf('Continuing...\n\n');
end

fprintf('=========================================\n');
fprintf('  Simulation pipeline complete.\n');
fprintf('  Figures saved to: %s\n', sim_save_dir);
fprintf('=========================================\n');

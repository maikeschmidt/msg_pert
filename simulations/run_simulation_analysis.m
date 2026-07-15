% run_simulation_analysis - Master script for the realistic-measurement pipeline
%
% Self-contained simulation package (msg_pert/simulations). The pt_* perturbation
% pipeline answers how wrong the forward field goes when the MODEL is wrong; this
% package answers a complementary question: how much does real SENSOR NOISE
% degrade a measured evoked response, and how does that differ between systems
% (SQUID MSG, OP-MSG, ESG) AND between geometry variants (the unperturbed model
% vs representative source / sensor / conductivity perturbations)?
%
% WORKFLOW (4 steps):
%   1. sim_plot_topoplots     Perfect forward fields (Biot-Savart MSG, BEM MSG,
%                             BEM ESG) for the chosen geometry variant, both
%                             arrays, every sensor axis and orientation.
%   2. sim_run_geometries     For EVERY geometry variant in sim_geometries
%                             (original + one representative shift per source /
%                             sensor / conductivity bundle) and every applicable
%                             system: simulate the evoked response, add
%                             trial-averaged sensor noise across the sweep, and
%                             score r^2 vs the noise-free field. One output
%                             subfolder per variant.
%   3. sim_plot_comparison    Overlay the r^2-vs-noise curves across variants,
%                             one figure per system: rows = perturbation family,
%                             lines = baseline / small / medium / large.
%   4. sim_plot_noise_topoplot  What one system actually measures at a chosen
%                             source and noise level, next to the noise-free field.
%
% To compare noise levels for a SINGLE geometry only, put just that one entry in
% sim_geometries — the loop handles the one-variant case identically.
%
% USAGE:
%   run_simulation_analysis
%
% CONFIGURATION:
%   Everything (paths, models, systems, noise floors, signal, and the geometry
%   variant list) lives in config_sim.m. Set those before running. In particular:
%     - the leadfield roots and geometry stems
%     - which representative shift per bundle to use (src_reps/sen_reps/cond_reps)
%     - bem_patched / the unit scales (a wrong scale flattens every curve)
%
% A NOTE ON COMPARING MSG WITH ESG:
%   MSG leadfields are fT/nAm and ESG leadfields are uV/nAm. Absolute noise
%   levels are NOT comparable across the two. Everything here is expressed in
%   multiples of each system's own baseline noise floor, which IS comparable.
%   The ESG floor (~1 uV/sqrt(Hz)) is an amplifier-noise estimate; real ESG is
%   dominated by much larger cardiac artefact, so ESG curves are an OPTIMISTIC
%   bound — state that in any write-up.
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

Metadata;

fprintf('  MSG Realistic-Measurement Simulation Pipeline\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n\n');

config_sim;

idxfile = fullfile(sim_out_dir, 'sim_geometry_index.mat');


% =========================================================================
% STEP 1: Perfect forward-field topoplots (chosen geometry variant)
% =========================================================================

fprintf('[1/4] Plotting perfect forward fields...\n');
try
    run('sim_plot_topoplots.m');
    fprintf('[1/4] Complete.\n\n');
catch err
    fprintf('WARNING: sim_plot_topoplots failed:\n  %s\n', err.message);
    fprintf('Check the model paths in config_sim.m. Continuing...\n\n');
end


% =========================================================================
% STEP 2: Simulate evoked + noise for every geometry variant
% =========================================================================

fprintf('[2/4] Simulating evoked response + noise across geometry variants...\n');
try
    run('sim_run_geometries.m');
    fprintf('[2/4] Complete.\n\n');
catch err
    fprintf('ERROR: sim_run_geometries failed:\n  %s\n', err.message);
    fprintf('Steps 3-4 depend on this — stopping.\n');
    return;
end

if ~isfile(idxfile)
    fprintf('ERROR: %s was not produced. Stopping.\n', idxfile);
    return;
end


% =========================================================================
% STEP 3: Comparison curves across variants
% =========================================================================

fprintf('[3/4] Plotting comparison curves across geometry variants...\n');
try
    run('sim_plot_comparison.m');
    fprintf('[3/4] Complete.\n\n');
catch err
    fprintf('WARNING: sim_plot_comparison failed:\n  %s\n', err.message);
    fprintf('Continuing...\n\n');
end


% =========================================================================
% STEP 4: Noisy measured topoplot at a chosen source and noise level
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

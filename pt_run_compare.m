function pt_run_compare()
% pt_run_compare - Run the combined MSG-vs-ESG comparison in an isolated workspace
%
% Thin wrapper around pt_compare_perturbations. It exists as its own function
% file (not a local function of run_perturbation_analysis) for two reasons:
%   1. Isolation — pt_compare_perturbations starts with clearvars; running it
%      through this function means that clearvars lands in this wrapper's
%      workspace, not the caller's, so run_perturbation_analysis keeps its loop
%      state (mods_to_run) for the final summary.
%   2. Interactive use — a function file is callable from the command window and
%      from run-section execution, whereas a script-local function is not.
%
% USAGE:
%   pt_run_compare
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

try
    run('pt_compare_perturbations.m');
    fprintf('Combined comparison complete.\n\n');
catch me
    fprintf('WARNING: pt_compare_perturbations failed:\n  %s\nContinuing...\n\n', ...
        me.message);
end
end

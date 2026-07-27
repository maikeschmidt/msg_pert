function pt_run_one_modality()
% pt_run_one_modality - Run Phase 2 analysis steps for the ACTIVE modality
%
% Runs the per-modality perturbation analysis (load → r² → curves → heatmaps →
% displacement → slope → tables) for whichever modality pt_modality currently
% selects. run_perturbation_analysis calls this once per modality.
%
% WORKSPACE ISOLATION
%   Every sub-script starts with `clearvars`, and run() executes a script in the
%   CALLER's workspace — so calling them directly here would let each sub-script
%   wipe this function's own variables (the *_rsq_file paths, the have_* flags)
%   between steps. To prevent that, every sub-script is launched through the
%   local run_script() wrapper: the clearvars then lands in the WRAPPER's
%   workspace, and this function's variables survive across all the steps.
%
%   Because of that isolation, config values are NOT left here by the sub-scripts;
%   this function calls config_pert once at the top to establish
%   forward_fields_base / save_base_dir for its own file checks. The active
%   modality (set by the master via pt_modality) is read by config_pert.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

config_pert;   % forward_fields_base / save_base_dir for the active modality
fprintf('  Results → %s\n\n', save_base_dir);

source_rsq_file = fullfile(forward_fields_base, 'pert_source_rsq.mat');
sensor_rsq_file = fullfile(forward_fields_base, 'pert_sensor_rsq.mat');
cond_rsq_file   = fullfile(forward_fields_base, 'pert_cond_rsq.mat');

% -------------------------------------------------------------------------
% STEP 1: Load and organise perturbed leadfields
% -------------------------------------------------------------------------
fprintf('[1/7] Loading and organising perturbed leadfields...\n');
if ~run_script('pt_load_leadfields.m')
    fprintf('ERROR: pt_load_leadfields failed — skipping the rest of this modality.\n');
    return;
end
fprintf('[1/7] Complete.\n\n');

% -------------------------------------------------------------------------
% STEP 2: Compute r²
% -------------------------------------------------------------------------
fprintf('[2/7] Computing perturbation r² (source, sensor, conductivity)...\n');
if ~run_script('pt_compute_rsq.m')
    fprintf('ERROR: pt_compute_rsq failed — skipping the rest of this modality.\n');
    return;
end
fprintf('[2/7] Complete.\n\n');

have_source_rsq = isfile(source_rsq_file);
have_sensor_rsq = isfile(sensor_rsq_file);
have_cond_rsq   = isfile(cond_rsq_file);
have_any_rsq    = have_source_rsq || have_sensor_rsq || have_cond_rsq;

% -------------------------------------------------------------------------
% STEP 3: r² curve figures
% -------------------------------------------------------------------------
if have_any_rsq
    fprintf('[3/7] Plotting perturbation curves...\n');
    run_script('pt_plot_curves.m');
    fprintf('[3/7] Done.\n\n');
else
    fprintf('[3/7] Skipping pt_plot_curves — no r² files found.\n\n');
end

% -------------------------------------------------------------------------
% STEP 4: Heatmap summaries
% -------------------------------------------------------------------------
if have_any_rsq
    fprintf('[4/7] Plotting heatmaps...\n');
    run_script('pt_plot_heatmaps.m');
    fprintf('[4/7] Done.\n\n');
else
    fprintf('[4/7] Skipping pt_plot_heatmaps — no r² files found.\n\n');
end

% -------------------------------------------------------------------------
% STEP 5: Displacement / perturbation vs r²
% -------------------------------------------------------------------------
if have_any_rsq
    fprintf('[5/7] Plotting displacement vs r²...\n');
    fprintf('       (individual figures: cervical region; combined + table: full cord)\n');
    run_script('pt_plot_displacement.m');
    fprintf('[5/7] Done.\n\n');
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

if isfile(sensor_tbl) || isfile(source_tbl) || isfile(cond_tbl)
    fprintf('[6/7] Plotting slope vs cord position (full cord)...\n');
    run_script('pt_plot_slope_vs_position.m');
    fprintf('[6/7] Done.\n\n');
else
    fprintf('[6/7] Skipping pt_plot_slope_vs_position — no trend tables found.\n\n');
end

% -------------------------------------------------------------------------
% STEP 7: Summary tables
% -------------------------------------------------------------------------
if have_any_rsq
    fprintf('[7/7] Computing summary tables...\n');
    run_script('pt_compute_table.m');
    fprintf('[7/7] Done.\n\n');
else
    fprintf('[7/7] Skipping pt_compute_table — no r² files found.\n\n');
end
end


% =========================================================================
% Local function: run a sub-script in an ISOLATED workspace
% =========================================================================
% The sub-script's clearvars lands in THIS wrapper's workspace, not in
% pt_run_one_modality, so the caller's variables survive across all steps.
%
% Note the ordering: `ok` is assigned AFTER run() returns, because the
% sub-script's clearvars would otherwise wipe it. In the catch, `me` is bound by
% the catch statement itself (after any clearvars), so it is always valid.
function ok = run_script(script_name)
    try
        run(script_name);
        ok = true;    % set after run() — clearvars inside has already happened
    catch me
        ok = false;
        fprintf('WARNING: a pipeline sub-script failed:\n  %s\n', me.message);
    end
end

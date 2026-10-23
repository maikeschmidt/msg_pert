% run_all_perturbations - Phase 1: perturbed geometries and their lead fields
%
% Generates every perturbed geometry and computes every lead field the
% analysis needs, for each modality in pert_modalities, using msg_fwd's
% solvers. When it finishes, run_all_pert_analysis runs the analysis.
%
% STEPS, per modality
%   1  Source shifts       pt_generate_source_shifts   25 geometries (original
%                                                      + 3 bundles x 8)
%   2  Sensor shifts       pt_generate_sensor_shifts   25 geometries
%   3  BEM                 pt_compute_leadfields       every geometry above
%   4  Biot-Savart (MSG)   pt_compute_leadfields       every geometry above
%   5  Conductivity (BEM)  pt_compute_leadfields       the original geometry,
%                                                      3 bundles x 8 draws
%   6  FEM (MSG)           pt_compute_leadfields       every geometry above
%   7  Conductivity (FEM)  pt_compute_leadfields       the same 24 draws as 5
%   8  FEM re-meshing      pt_compute_leadfields       the original geometry,
%                                                      one mesh per
%                                                      fem_remesh_seeds
%
%   pt_compute_leadfields follows msg_fwd's run_bem_leadfields,
%   run_biot_savart_leadfields, run_conductivity_perturbation and
%   run_fem_leadfields and writes msg_fwd's file layout, into the folders of
%   config_pert_paths. Only the array the analysis uses (config_pert:
%   analysis_array) is computed. FEM is MSG only: the DUNEuro solve msg_fwd
%   uses has no EEG route.
%
% USAGE:
%   run_all_perturbations
%
%   Set the flags below to run a subset. Every step skips files that already
%   exist, so an interrupted run resumes where it stopped, and re-running
%   after adding a method computes only what is missing.
%
% BEFORE RUNNING
%   - config_pert_paths: data_root, and the original geometry from msg_coreg
%     saved as <msg_geoms_dir>/<base_geom_name>.mat (and the ESG one in
%     esg_geoms_dir)
%   - config_pert: base_geom_name, pert_modalities, analysis_array, the
%     have_* flags per modality (they decide which methods are computed)
%   - SPM (FieldTrip), msg_coreg with HBF, and msg_fwd as sibling folders;
%     pt_add_functions checks all three
%   - for FEM: Windows, the DUNEuro binary (config_pert_paths:
%     duneuro_binpath), ISO2Mesh and fem_calc_fwds on the path, as for
%     msg_fwd's run_fem_leadfields. The FEM steps stop with a message if
%     any is missing; the BEM and Biot-Savart steps still run.
%
% COST
%   BEM: one head model per geometry (50 per modality) and one per
%   conductivity draw (24 per modality); hours. FEM: 25 + 8 meshes (source
%   geometries and re-meshes; sensor shifts and conductivity draws reuse the
%   unshifted mesh) and 50 + 24 + 8 DUNEuro solves; expect days, so run it
%   with do_fem = false first if you want the BEM results sooner.
%   Biot-Savart takes seconds per geometry.
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

% ---- Step flags ---------------------------------------------------------
do_geometries   = true;    % steps 1-2: generate shifted geometries
do_bem          = true;    % step 3
do_bslaw        = true;    % step 4 (MSG only)
do_conductivity = true;    % steps 5 and 7 (BEM and FEM draws)
do_fem          = true;    % steps 6-8 (MSG only)
overwrite       = false;   % true recomputes lead fields that already exist

fprintf('  MSG Perturbation Pipeline — Phase 1: geometries and lead fields\n');
fprintf('  University College London\n');
fprintf('  Department of Imaging Neuroscience\n\n');

pt_modality('clear');
config_pert;
pt_add_functions;
mods_to_run = pert_modalities;
t0 = tic;
status = cell(0, 2);

for mi = 1:numel(mods_to_run)
    md = mods_to_run{mi};
    pt_modality('set', md);
    config_pert;                        % this modality's paths and flags
    banner(sprintf('%s — geometries in %s', upper(md), perturbed_geoms_path));

    src_keys = [{sensitivity_ref_key}, sensitivity_keys];
    sen_keys = [{sensor_sensitivity_ref_key}, sensor_sensitivity_keys];

    % ---- Steps 1-2: geometries ------------------------------------------
    if do_geometries
        ok = try_step(@() pt_generate_source_shifts(), 'source-shift geometries');
        status(end+1, :) = {sprintf('%s source geometries', md), ok}; %#ok<SAGROW>
        ok = try_step(@() pt_generate_sensor_shifts(), 'sensor-shift geometries');
        status(end+1, :) = {sprintf('%s sensor geometries', md), ok}; %#ok<SAGROW>
        config_pert;                    % re-read: picks up the saved shift vectors
    end

    lopts = struct('arrays', {{analysis_array}}, 'overwrite', overwrite);
    all_keys = [src_keys, sen_keys];

    % ---- Step 3: BEM ----------------------------------------------------
    if do_bem && have_bem
        banner(sprintf('%s — BEM lead fields (%d geometries)', upper(md), numel(all_keys)));
        ok = try_step(@() report(pt_compute_leadfields('bem', perturbed_geoms_path, ...
            all_keys, bem_path, lopts)), 'BEM');
        status(end+1, :) = {sprintf('%s BEM', md), ok}; %#ok<SAGROW>
    end

    % ---- Step 4: Biot-Savart ---------------------------------------------
    if do_bslaw && have_bslaw && sensor_is_meg
        banner(sprintf('%s — Biot-Savart lead fields', upper(md)));
        ok = try_step(@() report(pt_compute_leadfields('bslaw', perturbed_geoms_path, ...
            all_keys, bslaw_path, lopts)), 'Biot-Savart');
        status(end+1, :) = {sprintf('%s Biot-Savart', md), ok}; %#ok<SAGROW>
    end

    % ---- Step 5: conductivity --------------------------------------------
    if do_conductivity && have_bem_cond
        banner(sprintf('%s — BEM conductivity perturbations', upper(md)));
        copts = lopts;
        copts.cond = struct('seed', cond_seed, 'bundle_pct', cond_bundle_pct, ...
            'n_shifts', n_cond_shifts, 'n_comp', n_cond_compartments);
        ok = try_step(@() report(pt_compute_leadfields('bem_cond', perturbed_geoms_path, ...
            cond_sensitivity_ref_key, bem_cond_path, copts)), 'conductivity');
        status(end+1, :) = {sprintf('%s conductivity', md), ok}; %#ok<SAGROW>
    end

    % ---- Steps 6-8: FEM ---------------------------------------------------
    if do_fem && have_fem && sensor_is_meg
        fopts = lopts;
        fopts.fem = struct('binpath', duneuro_binpath, 'maxvol_mm3', fem_maxvol_mm3, ...
            'reduce_torso', fem_reduce_torso, 'torso_keep', fem_torso_keep, ...
            'seed', fem_mesh_seed, 'remesh_seeds', fem_remesh_seeds, ...
            'cond', fem_cond_nominal);
        fopts.cond = struct('seed', cond_seed, 'bundle_pct', cond_bundle_pct, ...
            'n_shifts', n_cond_shifts, 'n_comp', n_cond_compartments);

        banner(sprintf('%s — FEM lead fields (%d geometries)', upper(md), numel(all_keys)));
        ok = try_step(@() report(pt_compute_leadfields('fem', perturbed_geoms_path, ...
            all_keys, fem_path, fopts)), 'FEM');
        status(end+1, :) = {sprintf('%s FEM', md), ok}; %#ok<SAGROW>

        if do_conductivity && have_fem_cond
            banner(sprintf('%s — FEM conductivity perturbations', upper(md)));
            ok = try_step(@() report(pt_compute_leadfields('fem_cond', perturbed_geoms_path, ...
                cond_sensitivity_ref_key, fem_path, fopts)), 'FEM conductivity');
            status(end+1, :) = {sprintf('%s FEM conductivity', md), ok}; %#ok<SAGROW>
        end

        if have_fem_remesh
            banner(sprintf('%s — FEM re-meshing control (%d seeds)', upper(md), ...
                numel(fem_remesh_seeds)));
            ok = try_step(@() report(pt_compute_leadfields('fem_remesh', perturbed_geoms_path, ...
                sensitivity_ref_key, fem_path, fopts)), 'FEM re-meshing');
            status(end+1, :) = {sprintf('%s FEM re-meshing', md), ok}; %#ok<SAGROW>
        end
    elseif do_fem && have_fem && ~sensor_is_meg
        fprintf('  %s: FEM skipped — the DUNEuro solve msg_fwd uses is MEG only.\n', upper(md));
    end
end

pt_modality('clear');
fprintf('\n=========================================\n');
fprintf('  Phase 1 finished in %.1f min\n', toc(t0) / 60);
for k = 1:size(status, 1)
    if status{k, 2}, s = 'ok'; else, s = 'FAILED'; end
    fprintf('  %-28s %s\n', status{k, 1}, s);
end
fprintf('  Next: run_all_pert_analysis\n');
fprintf('=========================================\n');


function ok = try_step(fn, label)
% Run one step; an error is reported and the next step still runs
    try
        fn();
        ok = true;
    catch me
        ok = false;
        fprintf(2, 'WARNING: %s failed:\n  %s\n', label, me.message);
        if ~isempty(me.stack)
            fprintf(2, '  at %s, line %d\n', me.stack(1).name, me.stack(1).line);
        end
    end
end

function report(R)
    fprintf('  written %d, already present %d, geometries missing %d\n', ...
        R.written, R.skipped, R.missing);
    if R.missing > 0
        error('%d geometry file(s) could not be computed (see warnings above).', ...
            R.missing);
    end
end

function banner(txt)
    fprintf('\n==================================================\n');
    fprintf('  %s\n', txt);
    fprintf('==================================================\n\n');
end

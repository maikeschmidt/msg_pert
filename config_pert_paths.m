% config_pert_paths - Every filesystem location used by msg_pert
%
% All paths live here. No other msg_pert script contains a hard-coded
% directory: config_pert (the analysis) and simulations/config_sim (the
% standalone noise package) both call this file and read the variables
% below. Change a drive letter or move the data once, here, and every script
% follows.
%
% The split matches msg_fwd's config_paths / config_models:
%   config_pert_paths   WHERE the data is
%   config_pert         WHAT is perturbed, analysed and how it is drawn
%
% The file is not called config_paths because msg_pert puts msg_fwd on the
% MATLAB path, and msg_fwd's config_paths would then shadow it.
%
% USAGE
%   Called automatically by config_pert and config_sim. Call it directly only
%   to see the paths without the rest of the configuration.
%
% LAYOUT UNDER data_root (the defaults; change any line below)
%   geometries/              MSG geometries: <base>.mat from msg_coreg, plus
%                            every geometries_<base>_source_* / _sensor_*
%                            file run_all_perturbations writes
%   geoms_elec/              the same for ESG
%   fields/mag/              MSG lead fields; also leadfields_organised.mat
%     bem/                     and pert_metrics.mat (stage 0)
%     bs_law/
%     bem_cond_msg/
%     fem/                     FEM: shifted geometries, conductivity draws and
%                              re-meshing control, plus meshes/ (cached
%                              volume meshes) and work/ (DUNEuro scratch)
%   fields/elec/             ESG lead fields
%     bem_elec/
%     bem_cond_esg/
%   fields/simulation/       simulations/ package working files
%   results/                 every figure and table (results_root)
%     staged/                  stages 1-6
%     msg/, esg/, combined/    earlier r2-only outputs (run_legacy)
%     simulation/              simulations/ package
%
% VARIABLES DEFINED
%   Roots:        data_root, results_root
%   Geometries:   msg_geoms_dir, esg_geoms_dir
%   Lead fields:  msg_fields_dir, msg_bem_dir, msg_bslaw_dir, msg_fem_dir,
%                 msg_sphere_dir, msg_bem_cond_dir, and esg_* likewise
%                 (esg_bslaw_dir and esg_fem_dir are always '': Biot-Savart
%                 is magnetic, and the DUNEuro solve msg_fwd uses is MEG only)
%   Binaries:     duneuro_binpath
%   Results:      stage_results_dir, msg_results_dir, esg_results_dir,
%                 combined_results_dir, sim_results_dir, sim_fields_dir
%   Helpers:      pert_geom_file(dir, key)   full path of geometries_<key>.mat
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk


% ROOTS — in most cases these are the only two lines that need changing

data_root    = 'D:\Thesis\chapter_6';      % SET THIS
results_root = fullfile(data_root);     % SET THIS (or leave under data_root)

% DUNEuro, for the FEM lead fields. The FOLDER holding
% bst_duneuro_meeg_win64.exe (not the executable), somewhere group policy
% lets executables run from. Same value as msg_fwd's config_paths.
duneuro_binpath = 'C:\wtcnapps\duneuro';           % SET THIS


% GEOMETRIES
% The original geometry from msg_coreg (<base_geom_name>.mat, see config_pert)
% and every perturbed copy live in one folder per modality.

msg_geoms_dir = fullfile(data_root, 'geometries');
esg_geoms_dir = fullfile(data_root, 'geoms_elec');


% LEAD FIELDS
% The folder layout inside each method folder is msg_fwd's: BEM and FEM keep
% one geometries_<key> subfolder per geometry, Biot-Savart and sphere are flat.
% Leave a method '' if it was not computed (and set its have_* flag in
% config_pert to false).

msg_fields_dir   = fullfile(data_root, 'fields', 'mag');
msg_bem_dir      = fullfile(msg_fields_dir, 'bem');
msg_bslaw_dir    = fullfile(msg_fields_dir, 'bs_law');
msg_bem_cond_dir = fullfile(msg_fields_dir, 'bem_cond_msg');
msg_fem_dir      = fullfile(msg_fields_dir, 'fem');   % '' if FEM is not computed
msg_sphere_dir   = '';

esg_fields_dir   = fullfile(data_root, 'fields', 'elec');
esg_bem_dir      = fullfile(esg_fields_dir, 'bem_elec');
esg_bem_cond_dir = fullfile(esg_fields_dir, 'bem_cond_esg');
esg_bslaw_dir    = '';                       % Biot-Savart is magnetic: none for ESG
esg_fem_dir      = '';                       % DUNEuro solve is MEG only: none for ESG
esg_sphere_dir   = '';


% RESULTS
% Every analysis writes to its own subfolder of results_root.

stage_results_dir    = fullfile(results_root, 'staged');
msg_results_dir      = fullfile(results_root, 'msg');
esg_results_dir      = fullfile(results_root, 'esg');
combined_results_dir = fullfile(results_root, 'combined');

% simulations/ package
sim_results_dir = fullfile(results_root, 'simulation');
sim_fields_dir  = fullfile(data_root, 'fields', 'simulation');


% PATH HELPERS
% Anonymous functions rather than files, so they are available to any script
% that has run config_pert_paths without adding anything to the MATLAB path.

pert_geom_file = @(dir, key) fullfile(dir, sprintf('geometries_%s.mat', key));

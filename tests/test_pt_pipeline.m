% test_pt_pipeline - End-to-end run of every stage on synthetic leadfields
%
% Builds a small synthetic data set in a temporary folder — a straight cord,
% a triaxial MSG array and a two-set ESG electrode array, with dipole fields
% in an unbounded medium plus a mirror-image term standing in for the volume
% conductor — then runs run_all_pert_analysis on it with every stage
% enabled. The files deliberately use the same mixture of unit conventions
% as real msg_fwd output (raw BEM per A*m, Biot-Savart already in fT/nAm,
% conductivity files with an extra x1e15), plus a decoy front-array file, so
% the loader's array filter and per-file scaling are exercised too.
%
% Checks:
%   - every stage finishes and writes its main outputs
%   - only back-array files are loaded, and every perturbed leadfield is
%     scaled to within 10x of its reference
%   - stage 5 at zero noise reproduces the stage 2 metrics exactly, for the
%     cord median and at the C7 source
%   - larger source shifts give larger RE (dose response), with an interval
%     that excludes zero
%   - the Biot-Savart noise systems exist and skip conductivity
%   - the anatomy-contribution table and the C7 ranking are written
%   - FEM lead fields, FEM conductivity draws and the FEM re-meshing control
%     are loaded and scored, and the FEM noise systems are simulated
%   - stage 7 writes its divergence, detectability and perturbation tables,
%     and the FEM recording differs from the BEM one
%
% USAGE:
%   test_pt_pipeline            % needs no toolboxes beyond MATLAB
%
% The temporary folder is printed and kept, so the figures can be inspected.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

clearvars
close all

root = fileparts(fileparts(mfilename('fullpath')));
tmp  = tempname;
mkdir(tmp);
fprintf('Synthetic test folder: %s\n', tmp);

% ---- Config pointing at the temporary folder, with fast settings --------
% config_pert_paths is copied with data_root moved to the temporary folder;
% config_pert is copied with fast settings appended. Both shadow the
% repository copies because cfg_dir goes to the top of the path.
ptxt = fileread(fullfile(root, 'config_pert_paths.m'));
[i0, i1] = regexp(ptxt, '^data_root\s*=[^\n]*', 'once', 'lineanchors');
assert(~isempty(i0), 'data_root line not found in config_pert_paths.m');
ptxt = [ptxt(1:i0-1), sprintf('data_root    = ''%s'';', tmp), ptxt(i1+1:end)];
txt = fileread(fullfile(root, 'config_pert.m'));
txt = [txt, sprintf(['\n\n%% ---- test overrides ----\n' ...
    'stats_n_boot = 200; noise_n_real = 4;\n' ...
    'noise_factors = [0.5, 1, 4]; fig_resolution = 50;\n' ...
    'baseline_topo_src_mm = 100; evoked_src_mm = 120;\n' ...
    'cord_fs = 5000; cord_n13_fwhm_mm = [0, 10]; cord_n13_head_fwhm = 10;\n'])];
cfg_dir = fullfile(tmp, 'cfg'); mkdir(cfg_dir);
fid = fopen(fullfile(cfg_dir, 'config_pert.m'), 'w'); fwrite(fid, txt); fclose(fid);
fid = fopen(fullfile(cfg_dir, 'config_pert_paths.m'), 'w'); fwrite(fid, ptxt); fclose(fid);

addpath(root);
addpath(cfg_dir);                     % shadows the repository config
pt_modality('clear');
config_pert;
pt_add_functions;

% ---- Synthetic geometry -------------------------------------------------
n_src  = 60;
src0   = [zeros(n_src, 1), (0:n_src-1)' * src_spacing_mm, zeros(n_src, 1)];   % mm
[gx, gy] = meshgrid(linspace(-60, 60, 5), linspace(-20, 320, 12));
msg_pos0 = [gx(:), gy(:), -60 * ones(numel(gx), 1)];
esg_pos0 = {[gx(:), gy(:), -50 * ones(numel(gx), 1)], ...
            [gx(:) + 6, gy(:) + 3, -55 * ones(numel(gx), 1)]};

src_vec = vertcat(source_shift_vectors{:});
sen_vec = vertcat(sensor_shift_vectors{:});
rng(cond_seed);
cond_d = cell(1, n_cond_bundles);
for b = 1:n_cond_bundles
    cond_d{b} = rand(n_cond_shifts, n_cond_compartments) * cond_bundle_pct(b);
end

geom_variants = {'original_source_original', zeros(1, 3), zeros(1, 3)
                 'original_sensor_original', zeros(1, 3), zeros(1, 3)};
for i = 1:size(src_vec, 1)
    b = ceil(i / n_source_shifts); s = i - (b-1)*n_source_shifts;
    geom_variants(end+1, :) = {sprintf('original_source_bundle%d_shift%d', b, s), src_vec(i, :), zeros(1, 3)}; %#ok<SAGROW>
    geom_variants(end+1, :) = {sprintf('original_sensor_bundle%d_shift%d', b, s), zeros(1, 3), sen_vec(i, :)}; %#ok<SAGROW>
end

for md = {'msg', 'esg'}
    C = mods_cfg.(md{1});
    for g = 1:size(geom_variants, 1)
        short = geom_variants{g, 1};
        ssh = geom_variants{g, 2}; nsh = geom_variants{g, 3};
        if strcmp(md{1}, 'esg'), nsh(3) = 0; end
        src = src0 + ssh;
        bem = synth(md{1}, src, msg_pos0, esg_pos0, nsh, 0.4, 1e9);     % raw, per A*m
        save_lf(fullfile(C.bem_path, ['geometries_' short]), ...
            sprintf('leadfield_%s_bem_back.mat', short), 'leadfield_cord', bem);
        if strcmp(md{1}, 'msg')
            bs = synth('msg', src, msg_pos0, esg_pos0, nsh, 0, 1e15);   % fT/nAm
            save_lf(C.bslaw_path, sprintf('leadfield_geometries_%s_bslaw_back.mat', short), ...
                'leadfield_bs', bs);
            % FEM: a slightly different volume-conductor weight, in fT/nAm
            fem = synth('msg', src, msg_pos0, esg_pos0, nsh, 0.45, 1e15);
            save_lf(fullfile(C.fem_path, ['geometries_' short]), ...
                sprintf('cord_leadfield_%s_fem_back.mat', short), 'leadfield_ft', fem);
        end
        if g == 1
            % Decoy front-array file: must never be loaded
            save_lf(fullfile(C.bem_path, ['geometries_' short]), ...
                sprintf('leadfield_%s_bem_front.mat', short), 'leadfield_cord', scale_lf(bem, 1e6));
            geom = struct();
            if strcmp(md{1}, 'msg')
                geom.back_coils_3axis.chanpos = repmat(msg_pos0, 3, 1);
            else
                geom.back_coils_2axis.elecpos = vertcat(esg_pos0{:});
            end
            if ~exist(C.geoms_path, 'dir'), mkdir(C.geoms_path); end
            save(fullfile(C.geoms_path, ['geometries_' short '.mat']), '-struct', 'geom');
        end
    end
    % Conductivity: changes the volume-conductor weight and the overall gain
    for b = 1:n_cond_bundles
        for s = 1:n_cond_shifts
            d = cond_d{b}(s, :);
            lf = synth(md{1}, src0, msg_pos0, esg_pos0, zeros(1, 3), ...
                 0.4 * (1 + d(1) - d(2)), (1 - 0.3 * mean(d)) * 1e9 * 1e15);   % baked x1e15
            save_lf(fullfile(C.bem_cond_path, 'geometries_original_source_original'), ...
                sprintf('leadfield_original_source_original_bem_cond_bundle%d_shift%d_back.mat', b, s), ...
                'leadfield_cord', lf);
            if strcmp(md{1}, 'msg')
                lf = synth('msg', src0, msg_pos0, esg_pos0, zeros(1, 3), ...
                     0.45 * (1 + d(1) - d(2)), (1 - 0.3 * mean(d)) * 1e15);
                save_lf(fullfile(C.fem_path, 'geometries_original_source_original'), ...
                    sprintf('cord_leadfield_original_source_original_fem_cond_bundle%d_shift%d_back.mat', b, s), ...
                    'leadfield_ft', lf);
            end
        end
    end
    % FEM re-meshing: the unshifted FEM with a small random perturbation
    if strcmp(md{1}, 'msg')
        fem0 = synth('msg', src0, msg_pos0, esg_pos0, zeros(1, 3), 0.45, 1e15);
        for k = 1:numel(fem_remesh_seeds)
            lf = fem0;
            lf.leadfield = cellfun(@(x) x .* (1 + 0.01 * randn(size(x))), lf.leadfield, ...
                'UniformOutput', false);
            save_lf(fullfile(C.fem_path, 'geometries_original_source_original'), ...
                sprintf('cord_leadfield_original_source_original_fem_remesh%d_back.mat', k), ...
                'leadfield_ft', lf);
        end
    end
end
fprintf('Synthetic leadfields written.\n\n');

% ---- Run every stage ----------------------------------------------------
run_all_pert_analysis;

% ---- Checks -------------------------------------------------------------
% run_all_pert_analysis starts with clearvars, so re-derive what is needed
config_pert;
tmp = fileparts(fileparts(stage_results_dir));
n_fail = 0;
expect = {fullfile('1_baseline', 'baseline_model_type.csv'), ...
          fullfile('1_baseline', 'amplitude_profiles.png'), ...
          fullfile('2_within_msg', 'within_stats.csv'), ...
          fullfile('2_within_msg', 'methods_compared.png'), ...
          fullfile('3_within_esg', 'within_descriptive.tex'), ...
          fullfile('4_msg_vs_esg', 'modality_stats.csv'), ...
          fullfile('4_msg_vs_esg', 'modality_interaction.csv'), ...
          fullfile('2_within_msg', 'within_anatomy.csv'), ...
          fullfile('2_within_msg', 'anatomy_contribution_bem.png'), ...
          fullfile('2_within_msg', 'anatomy_contribution_fem.png'), ...
          fullfile('2_within_msg', 'within_fem_remesh.csv'), ...
          fullfile('2_within_msg', 'fem_remesh_floor.png'), ...
          fullfile('5_noise', 'data', 'noise_squid_msg_fem.mat'), ...
          fullfile('5_noise', 'c7', 'noise_combined.csv'), ...
          fullfile('5_noise', 'c7', 'noise_system_stats.csv'), ...
          fullfile('5_noise', 'cord', 'noise_combined.csv'), ...
          fullfile('5_noise', 'stage5_long.csv'), ...
          fullfile('5_noise', 'data', 'noise_squid_msg_bs.mat'), ...
          fullfile('6_summary', 'summary_msg.tex'), ...
          fullfile('6_summary', 'summary_msg_c7.tex'), ...
          fullfile('7_cord_activity', 'cord_divergence.csv'), ...
          fullfile('7_cord_activity', 'cord_critical_levels.csv'), ...
          fullfile('7_cord_activity', 'cord_bem_vs_fem_noise.csv'), ...
          fullfile('7_cord_activity', 'cord_perturbation_bem_vs_fem.csv'), ...
          fullfile('7_cord_activity', 'source_model.png')};
for k = 1:numel(expect)
    n_fail = n_fail + report(['exists: ' expect{k}], ...
        double(isfile(fullfile(stage_results_dir, expect{k}))), 1, 0);
end

L = load(fullfile(mods_cfg.msg.forward_fields_base, 'leadfields_organised.mat'), 'scale_log');
n_fail = n_fail + report('no front-array file loaded', ...
    double(~any(contains(L.scale_log.file, '_front'))), 1, 0);

PM = load(fullfile(mods_cfg.msg.forward_fields_base, 'pert_metrics.mat'), 'PM'); PM = PM.PM;
lnm = PM.pert.cond.M.bem.lnmag(:, :, end, end);
n_fail = n_fail + report('cond files scaled to their reference', ...
    double(max(abs(lnm(:))) < log(10)), 1, 0);

NS = load(fullfile(stage_results_dir, '5_noise', 'data', 'noise_squid_msg.mat'), 'NS'); NS = NS.NS;
v2 = pt_shift_summary(PM.pert.source.M.bem, 're', numel(PM.axcfg_ids), numel(PM.ori_labels));
v5 = NS.pert.source.cm.re(2:end, 1, 1, end, end);
n_fail = n_fail + report('noise level 0 == stage 2 (max |diff|)', max(abs(v2 - v5)), 0, 1e-9);
c2 = PM.pert.source.M.bem.re(:, NS.c7_idx, end, end);
c5 = NS.pert.source.c7.re(2:end, 1, 1, end, end);
n_fail = n_fail + report('C7 level 0 == stage 2 at C7 (max |diff|)', max(abs(c2 - c5)), 0, 1e-9);
n_fail = n_fail + report('C7 source is the nearest to evoked_src_mm', ...
    abs(NS.c7_mm - evoked_src_mm), 0, src_spacing_mm / 2);

n_fail = n_fail + report('FEM conductivity draws scored', ...
    double(ismember('fem', PM.pert.cond.methods)), 1, 0);
n_fail = n_fail + report('FEM re-meshing floor scored, not a type', ...
    double(isfield(PM, 'remesh') && isfield(PM.remesh.M, 'fem') && ~isfield(PM.pert, 'remesh')), 1, 0);
NF = load(fullfile(stage_results_dir, '5_noise', 'data', 'noise_squid_msg_fem.mat'), 'NS'); NF = NF.NS;
n_fail = n_fail + report('FEM noise system has conductivity rows', ...
    double(isfield(NF.pert, 'cond')), 1, 0);

NB = load(fullfile(stage_results_dir, '5_noise', 'data', 'noise_squid_msg_bs.mat'), 'NS'); NB = NB.NS;
n_fail = n_fail + report('Biot-Savart system has no conductivity rows', ...
    double(~isfield(NB.pert, 'cond')), 1, 0);

T = readtable(fullfile(stage_results_dir, '2_within_msg', 'within_stats.csv'), 'TextType', 'char');
r = T(strcmp(T.family, 'dose_response') & strcmp(T.subset, 'Source space') & ...
      strcmp(T.metric, 'RE (%)') & strcmp(T.method, 'BEM') & ...
      strcmp(T.orientation, 'ALL') & strcmp(T.axis_set, 'Array'), :);
n_fail = n_fail + report('source dose response rho > 0.5', double(r.estimate(1) > 0.5), 1, 0);
n_fail = n_fail + report('source dose response interval excludes 0', double(r.separated(1)), 1, 0);

T7 = readtable(fullfile(stage_results_dir, '7_cord_activity', 'cord_divergence.csv'), 'TextType', 'char');
r7 = T7(strcmp(T7.comparison, 'FEM vs BEM') & strcmp(T7.window, 'N13'), :);
n_fail = n_fail + report('stage 7: FEM differs from BEM (N13 RE > 0)', ...
    double(~isempty(r7) && all(r7.re > 0)), 1, 0);

L5 = readtable(fullfile(stage_results_dir, '5_noise', 'stage5_long.csv'), 'TextType', 'char');
n_fail = n_fail + report('stage 5 long rows carry both sites', ...
    double(all(ismember({'C7', 'cord'}, unique(L5.site)))), 1, 0);

fprintf('\n%d failure(s). Outputs kept in %s\n', n_fail, tmp);


% =========================================================================
function L = synth(md, src_mm, msg_pos0, esg_pos0, sensor_shift, img_w, scale)
% Synthetic leadfield, per nA*m times scale, [n_chan x 3] per source (cols LR RC VD).
% Unbounded-medium dipole field plus a weighted mirror image of the source
% in the plane z = -35 mm, standing in for the volume conductor boundary.
    q = 1e-9;                                   % 1 nA*m in A*m
    img = @(p) [p(:, 1:2), -70 - p(:, 3)];       % mirror in z = -35 mm
    L = struct('leadfield', {cell(1, size(src_mm, 1))});
    for s = 1:size(src_mm, 1)
        if strcmp(md, 'msg')
            P = (msg_pos0 + sensor_shift) / 1000;
            blocks = cell(3, 1);
            for comp = 1:3, blocks{comp} = zeros(size(P, 1), 3); end
            for d = 1:3
                qv = zeros(1, 3); qv(d) = q;
                B = bfield(P, src_mm(s, :) / 1000, qv) + ...
                    img_w * bfield(P, img(src_mm(s, :)) / 1000, qv .* [1 1 -1]);
                for comp = 1:3, blocks{comp}(:, d) = B(:, comp); end
            end
            L.leadfield{s} = vertcat(blocks{:}) * scale;
        else
            blocks = cell(2, 1);
            for a = 1:2
                P = (esg_pos0{a} + sensor_shift) / 1000;
                V = zeros(size(P, 1), 3);
                for d = 1:3
                    qv = zeros(1, 3); qv(d) = q;
                    V(:, d) = vfield(P, src_mm(s, :) / 1000, qv) + ...
                        img_w * vfield(P, img(src_mm(s, :)) / 1000, qv .* [1 1 -1]);
                end
                blocks{a} = V;
            end
            L.leadfield{s} = vertcat(blocks{:}) * scale;
        end
    end
end

function L = scale_lf(L, f)
    L.leadfield = cellfun(@(x) x * f, L.leadfield, 'UniformOutput', false);
end

function B = bfield(P, r0, q)       % Biot-Savart, tesla
    R = P - r0; n = sqrt(sum(R.^2, 2)).^3;
    B = 1e-7 * cross(repmat(q, size(R, 1), 1), R, 2) ./ n;
end

function V = vfield(P, r0, q)       % unbounded-medium potential, volt (sigma 0.3 S/m)
    R = P - r0; n = sqrt(sum(R.^2, 2)).^3;
    V = (R * q') ./ (4 * pi * 0.3 * n);
end

function save_lf(d, f, var, L)
    if ~exist(d, 'dir'), mkdir(d); end
    S.(var) = L; %#ok<STRNU>
    save(fullfile(d, f), '-struct', 'S');
end

function f = report(name, got, want, tol)
    ok = abs(got - want) <= tol;
    if ok, s = 'OK'; else, s = '*** FAIL ***'; end
    fprintf('%-46s got %12.6g  want %12.6g  %s\n', name, got, want, s);
    f = double(~ok);
end

function R = pt_compute_leadfields(method, geom_dir, keys, out_dir, opts)
% pt_compute_leadfields - Lead fields for perturbed geometries, using msg_fwd
%
% Computes the lead fields msg_pert analyses, with msg_fwd's own solvers and
% in msg_fwd's file layout, so pt_load_leadfields reads them exactly as it
% reads files written by msg_fwd's runner scripts. The runners themselves
% are scripts with their own hard-coded folders and geometry lists; this
% function performs the same steps with the folders and keys msg_pert
% supplies:
%
%   'bem'       msg_fwd run_bem_leadfields: bem_build_boundaries (torso kept
%               at 0.5, bone ratio 40, conductivity pairing), HBF head model,
%               ft_prepare_leadfield with dipoleunit nA*m. Saved raw, as the
%               runner does; pt_unit_scale scales it on load.
%   'bslaw'     msg_fwd run_biot_savart_leadfields: bs_leadfield, already in
%               fT/nAm. MSG (coil) arrays only.
%   'bem_cond'  msg_fwd run_conductivity_perturbation: one geometry, every
%               compartment conductivity scaled by (1 + delta), delta drawn
%               U(0, pct) per bundle with a fixed seed (the draw
%               pt_compute_metrics replays). Deltas are compartments
%               innermost first: cord, bone, heart, lungs, torso. The
%               pairing is bem_build_boundaries': every outer conductivity is
%               torso (the bone's included), except the cord's when the bone
%               is continuous and the cord lies inside it, where it is bone.
%               Segmented bone leaves a gap round the cord, so the cord's
%               outer conductivity stays torso even where the bone encloses
%               it. The outer values follow the perturbed compartment, so a
%               zero perturbation reproduces the 'bem' lead field exactly.
%   'fem'       msg_fwd run_fem_leadfields: fem_build_volume_mesh (TetGen,
%               500 mm^3 bound, torso kept at 0.5), DUNEuro via
%               fem_calc_fwds, x1e6 to fT/nAm, convert_duneuro_to_fieldtrip.
%               MSG (coil) arrays only: the DUNEuro call is MEG only.
%   'fem_cond'  the 'fem' solve on the original geometry with every
%               compartment conductivity scaled by the same (1 + delta) draws
%               as 'bem_cond', so the two solvers see identical changes. The
%               FEM labels volumes, so there are no inner/outer pairs.
%   'fem_remesh' the original geometry meshed again with different TetGen
%               region seeds (opts.fem.remesh_seeds) and solved unchanged:
%               the FEM's re-meshing floor, against which FEM source shifts
%               are read.
%
% FEM MESHES
%   Volume meshes are cached in <out_dir>/meshes. A sensor shift or a
%   conductivity draw changes no surface, so it is solved on the mesh of the
%   unshifted geometry and its FEM difference contains no re-meshing. A
%   source shift moves the cord and bone, so it needs its own mesh; every
%   such mesh uses the same seed (opts.fem.seed), and 'fem_remesh' measures
%   how much of a FEM source-shift effect re-meshing alone could produce.
%
% Files that already exist are skipped, so an interrupted run resumes. A
% geometry that fails (TetGen, a non-nested BEM surface) is reported and
% skipped; the remaining geometries still run.
%
% USAGE:
%   R = pt_compute_leadfields('bem',      geom_dir, keys, out_dir)
%   R = pt_compute_leadfields('bslaw',    geom_dir, keys, out_dir, opts)
%   R = pt_compute_leadfields('bem_cond', geom_dir, ref_key, out_dir, opts)
%   R = pt_compute_leadfields('fem',      geom_dir, keys,    out_dir, opts)
%   R = pt_compute_leadfields('fem_cond' | 'fem_remesh', geom_dir, ref_key, out_dir, opts)
%
% INPUT:
%   method    - 'bem' | 'bslaw' | 'bem_cond' | 'fem' | 'fem_cond' | 'fem_remesh'
%   geom_dir  - folder holding geometries_<key>.mat
%   keys      - cellstr of geometry keys (without 'geometries_'); for
%               'bem_cond' the single reference key
%   out_dir   - method folder (config_pert_paths: msg_bem_dir, ...)
%   opts      - struct, every field optional:
%     .arrays     {'back'}   arrays to compute ('front', 'back',
%                            'experimental')
%     .overwrite  false      recompute files that exist
%     .cond       struct for 'bem_cond' and 'fem_cond': .seed, .bundle_pct,
%                            .n_shifts, .n_comp (config_pert: cond_seed,
%                            cond_bundle_pct, n_cond_shifts, n_cond_compartments)
%     .fem        struct for the FEM methods (config_pert: fem_*):
%                   .binpath       DUNEuro folder (config_pert_paths)
%                   .maxvol_mm3    500    .opt_scale 1
%                   .reduce_torso  true   .torso_keep 0.5
%                   .seed          1      region seed for every production mesh
%                   .remesh_seeds  []     seeds for 'fem_remesh'
%                   .cond          [0.33 0.33/40 0.62 0.05 0.23] S/m
%
% OUTPUT:
%   R - struct: .written .skipped .missing (counts of files / geometries),
%       .failed with .failed_keys and .failed_msgs: geometries whose solve
%       raised an error. Each is reported and skipped; the rest still run.
%
% FILES WRITTEN (msg_fwd naming):
%   bem       <out>/geometries_<key>/leadfield_<key>_bem_<array>.mat       leadfield_cord
%   bslaw     <out>/leadfield_geometries_<key>_bslaw_<array>.mat           leadfield_bs
%   bem_cond  <out>/geometries_<key>/leadfield_<key>_bem_cond_bundle<B>_shift<S>_<array>.mat
%                                                                          leadfield_cord
%   fem       <out>/geometries_<key>/cord_leadfield_<key>_fem_<array>.mat  leadfield_ft
%   fem_cond  <out>/geometries_<key>/cord_leadfield_<key>_fem_cond_bundle<B>_shift<S>_<array>.mat
%   fem_remesh <out>/geometries_<key>/cord_leadfield_<key>_fem_remesh<K>_<array>.mat
%   and compute_time_log.csv in out_dir (log_compute_time), as msg_fwd does.
%
% DEPENDENCIES:
%   msg_fwd/functions: bem_build_boundaries, bs_leadfield, log_compute_time,
%   fem_build_volume_mesh, convert_duneuro_to_fieldtrip
%   FieldTrip (via SPM) and HBF (via msg_coreg) for the BEM methods;
%   additionally ISO2Mesh, fem_calc_fwds and the DUNEuro binary (Windows)
%   for the FEM methods
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'arrays'),    opts.arrays    = {'back'}; end
if ~isfield(opts, 'overwrite'), opts.overwrite = false;    end
if ischar(keys), keys = {keys}; end

R = struct('written', 0, 'skipped', 0, 'missing', 0, 'failed', 0, ...
           'failed_keys', {{}}, 'failed_msgs', {{}});
if ~any(strcmp(method, {'bem', 'bslaw', 'bem_cond', 'fem', 'fem_cond', 'fem_remesh'}))
    error('pt_compute_leadfields:method', 'Unknown method "%s".', method);
end
if isempty(out_dir)
    warning('pt_compute_leadfields:noOutDir', ...
        'No output folder for %s (config_pert_paths) — nothing computed.', method);
    return
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
time_log = fullfile(out_dir, 'compute_time_log.csv');

if startsWith(method, 'fem')
    F = default_fem(getfield_or(opts, 'fem', struct()));
    mesh_dir = fullfile(out_dir, 'meshes');
    if ~exist(mesh_dir, 'dir'), mkdir(mesh_dir); end
    if strcmp(method, 'fem_remesh'), keys = keys(1); end
end

if any(strcmp(method, {'bem_cond', 'fem_cond'}))
    keys = keys(1);
    C = opts.cond;
    if ~isfield(C, 'n_comp'), C.n_comp = 5; end
    % The same draw as msg_fwd/run_conductivity_perturbation and the replay
    % in pt_compute_metrics: one rand(n_shifts, 5) per bundle, in order
    rng(C.seed);
    deltas = cell(1, numel(C.bundle_pct));
    for b = 1:numel(C.bundle_pct)
        deltas{b} = rand(C.n_shifts, C.n_comp) * C.bundle_pct(b);
    end
end

for g = 1:numel(keys)
    try
        key  = keys{g};
        gfile = fullfile(geom_dir, ['geometries_' key '.mat']);
        fprintf('  [%s %d/%d] %s\n', method, g, numel(keys), key);
        if ~isfile(gfile)
            warning('pt_compute_leadfields:noGeometry', 'Geometry not found: %s', gfile);
            R.missing = R.missing + 1;
            continue
        end

        % Output files for this geometry, decided before anything is loaded so
        % finished geometries cost nothing
        switch method
            case 'bem'
                sub   = fullfile(out_dir, ['geometries_' key]);
                fname = @(arr, b, s) fullfile(sub, sprintf('leadfield_%s_bem_%s.mat', key, arr));
                jobs  = [0 0];
            case 'bslaw'
                sub   = out_dir;
                fname = @(arr, b, s) fullfile(sub, sprintf('leadfield_geometries_%s_bslaw_%s.mat', key, arr));
                jobs  = [0 0];
            case 'bem_cond'
                sub   = fullfile(out_dir, ['geometries_' key]);
                fname = @(arr, b, s) fullfile(sub, sprintf( ...
                    'leadfield_%s_bem_cond_bundle%d_shift%d_%s.mat', key, b, s, arr));
                [bb, ss] = ndgrid(1:numel(C.bundle_pct), 1:C.n_shifts);
                jobs  = sortrows([bb(:), ss(:)]);
            case 'fem'
                sub   = fullfile(out_dir, ['geometries_' key]);
                fname = @(arr, b, s) fullfile(sub, sprintf('cord_leadfield_%s_fem_%s.mat', key, arr));
                jobs  = [0 0];
            case 'fem_cond'
                sub   = fullfile(out_dir, ['geometries_' key]);
                fname = @(arr, b, s) fullfile(sub, sprintf( ...
                    'cord_leadfield_%s_fem_cond_bundle%d_shift%d_%s.mat', key, b, s, arr));
                [bb, ss] = ndgrid(1:numel(C.bundle_pct), 1:C.n_shifts);
                jobs  = sortrows([bb(:), ss(:)]);
            case 'fem_remesh'
                sub   = fullfile(out_dir, ['geometries_' key]);
                fname = @(arr, b, s) fullfile(sub, sprintf( ...
                    'cord_leadfield_%s_fem_remesh%d_%s.mat', key, b, arr));
                jobs  = [(1:numel(F.remesh_seeds))', F.remesh_seeds(:)];   % [k, seed]
            otherwise
                error('pt_compute_leadfields:method', 'Unknown method "%s".', method);
        end
        todo = false(size(jobs, 1), numel(opts.arrays));
        for j = 1:size(jobs, 1)
            for a = 1:numel(opts.arrays)
                todo(j, a) = opts.overwrite || ~isfile(fname(opts.arrays{a}, jobs(j, 1), jobs(j, 2)));
            end
        end
        R.skipped = R.skipped + sum(~todo(:));
        if ~any(todo(:)), fprintf('    all files exist — skipped\n'); continue; end
        if ~exist(sub, 'dir'), mkdir(sub); end

        geoms = load(gfile);
        [sens, labels, is_elec] = find_arrays(geoms, opts.arrays, key);
        if isempty(sens), R.missing = R.missing + 1; continue; end
        pos_mm = geoms.sources_cent.pos;

        % ---- Biot-Savart -------------------------------------------------------
        if strcmp(method, 'bslaw')
            if is_elec
                warning('pt_compute_leadfields:bslawElec', ...
                    'Biot-Savart needs coil arrays; %s has electrodes — skipped.', key);
                R.missing = R.missing + 1;
                continue
            end
            for a = 1:numel(sens)
                if ~todo(1, strcmp(opts.arrays, labels{a})), continue; end
                grad = sens{a};
                leadfield_bs = struct();
                t_lf = tic;
                leadfield_bs.leadfield = bs_leadfield(grad, pos_mm);
                t_lf = toc(t_lf);
                leadfield_bs.label     = grad.label;
                leadfield_bs.pos       = pos_mm;
                leadfield_bs.unit      = 'mm';
                leadfield_bs.model     = 'biot_savart_infinite';
                leadfield_bs.geometry  = ['geometries_' key];
                leadfield_bs.array     = labels{a};
                leadfield_bs.mu0       = 4 * pi * 1e-7;
                leadfield_bs.units_out = 'fT/nAm';
                log_compute_time(time_log, struct('model', key, 'method', 'biot_savart', ...
                    'array', labels{a}, 'stage', 'leadfield', 'seconds', t_lf, ...
                    'n_sources', size(pos_mm, 1), 'n_channels', size(grad.tra, 1)));
                save(fname(labels{a}, 0, 0), 'leadfield_bs', '-v7.3');
                R.written = R.written + 1;
            end
            continue
        end

        % ---- FEM ---------------------------------------------------------------
        if startsWith(method, 'fem')
            if is_elec
                warning('pt_compute_leadfields:femElec', ...
                    ['The DUNEuro solve msg_fwd uses is MEG only; %s has ' ...
                     'electrodes — FEM skipped.'], key);
                R.missing = R.missing + 1;
                continue
            end
            src = struct('pos', pos_mm, 'inside', ones(size(pos_mm, 1), 1), 'unit', 'mm');
            src = ft_convert_units(src, 'm');
            grads_m = cellfun(@(x) ft_convert_units(x, 'm'), sens, 'UniformOutput', false);
            mesh_key = fem_mesh_key(key);
            for j = 1:size(jobs, 1)
                if ~any(todo(j, :)), continue; end
                switch method
                    case 'fem'
                        seed = F.seed; cond = F.cond; tag = '';
                    case 'fem_cond'
                        b = jobs(j, 1); s = jobs(j, 2);
                        d = deltas{b}(s, :);
                        seed = F.seed; cond = max(F.cond .* (1 + d), 1e-4);
                        tag  = sprintf('_cond_bundle%d_shift%d', b, s);
                    case 'fem_remesh'
                        seed = jobs(j, 2); cond = F.cond;
                        tag  = sprintf('_remesh%d', jobs(j, 1));
                end
                [tet, minfo] = fem_mesh(mesh_dir, geom_dir, mesh_key, seed, F, time_log, ...
                    size(pos_mm, 1));
                for a = 1:numel(grads_m)
                    ai = strcmp(opts.arrays, labels{a});
                    if ~todo(j, ai), continue; end
                    work = fullfile(out_dir, 'work', [key tag], labels{a});
                    if exist(work, 'dir'), rmdir(work, 's'); end   % DUNEuro needs a fresh folder
                    mkdir(work);
                    S = struct('dir', work, 'mesh', tet, 'grad', grads_m{a}, 'src', src, ...
                        'cond', cond, 'binpath', F.binpath);
                    t_solve = tic;
                    Lfem = fem_calc_fwds(S) * 1e6;                  % T/(A*m) -> fT/nAm
                    t_solve = toc(t_solve);
                    log_compute_time(time_log, struct('model', [key tag], 'method', 'fem', ...
                        'array', labels{a}, 'stage', 'solve', 'seconds', t_solve, ...
                        'n_sources', size(pos_mm, 1), 'n_channels', size(Lfem, 1), ...
                        'n_nodes', minfo.n_nodes, 'n_elements', minfo.n_tets));
                    leadfield_ft = convert_duneuro_to_fieldtrip(Lfem, src, grads_m{a}, S);
                    leadfield_ft.units_out  = 'fT/nAm';
                    leadfield_ft.geometry   = ['geometries_' key];
                    leadfield_ft.mesh_from  = mesh_key;
                    leadfield_ft.mesh_seed  = seed;
                    leadfield_ft.fem_cond   = cond;
                    if strcmp(method, 'fem_cond'), leadfield_ft.cond_delta_pct = d * 100; end
                    save(fname(labels{a}, jobs(j, 1), jobs(j, 2)), 'leadfield_ft', '-v7.3');
                    R.written = R.written + 1;
                end
            end
            continue
        end

        % ---- BEM ---------------------------------------------------------------
        [bnd, cond0, info] = bem_build_boundaries(geoms, ...
            struct('torso_keep', 0.5, 'cratio', 40, 'label', key));
        src = struct('pos', pos_mm, 'inside', true(size(pos_mm, 1), 1), 'unit', 'mm');
        src = ft_convert_units(src, 'm');
        sens_m = cellfun(@(x) ft_convert_units(x, 'm'), sens, 'UniformOutput', false);

        for j = 1:size(jobs, 1)
            if ~any(todo(j, :)), continue; end
            if strcmp(method, 'bem_cond')
                b = jobs(j, 1); s = jobs(j, 2);
                d = deltas{b}(s, :);
                % Each compartment's conductivity scaled; every outer value is
                % the conductivity of the compartment outside that boundary, so
                % it follows the same scaling: torso outside every boundary
                % (bone included, and the cord when a gap separates it from the
                % bone), bone outside the cord only for continuous bone with the
                % cord inside it (bem_build_boundaries' test), nothing outside
                % the torso
                ci = max(cond0(1, :) .* (1 + d), 1e-4);
                co = repmat(ci(end), 1, numel(ci));
                co(end) = 0;
                if info.bone_encloses_cord, co(1) = ci(2); end
                conductivity = [ci; co];
                fprintf('    bundle %d shift %d: +%s %%\n', b, s, ...
                    strjoin(arrayfun(@(x) sprintf('%.1f', x), d * 100, 'UniformOutput', false), '/'));
            else
                conductivity = cond0;
            end

            t_hm = tic;
            vol  = ft_prepare_headmodel(struct('method', 'hbf', ...
                'conductivity', conductivity, 'checkmesh', 'false'), bnd);
            t_hm = toc(t_hm);
            log_compute_time(time_log, struct('model', key, 'method', method, ...
                'stage', 'headmodel', 'seconds', t_hm, 'n_sources', size(pos_mm, 1), ...
                'n_nodes', info.n_vertices, 'n_elements', info.n_triangles));

            for a = 1:numel(sens_m)
                ai = strcmp(opts.arrays, labels{a});
                if ~todo(j, ai), continue; end
                cfg = struct('sourcemodel', src, 'headmodel', vol, 'reducerank', 'no', ...
                    'channel', 'all', 'normalize', 'no', 'dipoleunit', 'nA*m');
                if is_elec, cfg.elec = sens_m{a}; else, cfg.grad = sens_m{a}; end
                t_lf = tic;
                leadfield_cord = ft_prepare_leadfield(cfg);
                t_lf = toc(t_lf);
                log_compute_time(time_log, struct('model', key, 'method', method, ...
                    'array', labels{a}, 'stage', 'leadfield', 'seconds', t_lf, ...
                    'n_sources', size(pos_mm, 1), ...
                    'n_channels', size(leadfield_cord.leadfield{1}, 1), ...
                    'n_nodes', info.n_vertices, 'n_elements', info.n_triangles));
                leadfield_cord.geometry = ['geometries_' key];
                leadfield_cord.array    = labels{a};
                if strcmp(method, 'bem_cond')
                    leadfield_cord.model          = 'bem_cond';
                    leadfield_cord.cond_bundle    = b;
                    leadfield_cord.cond_shift     = s;
                    leadfield_cord.nominal_cond   = cond0;
                    leadfield_cord.perturbed_cond = conductivity;
                    leadfield_cord.cond_delta_pct = d * 100;
                end
                save(fname(labels{a}, jobs(j, 1), jobs(j, 2)), 'leadfield_cord', '-v7.3');
                R.written = R.written + 1;
            end
        end
    catch me
        % One geometry that cannot be solved (a bad mesh, a surface pushed
        % through another) must not stop the rest: report it and go on
        R.failed = R.failed + 1;
        R.failed_keys{end+1} = key;
        R.failed_msgs{end+1} = me.message;
        fprintf(2, '    FAILED %s: %s\n', key, me.message);
    end
end
end


function [sens, labels, is_elec] = find_arrays(geoms, wanted, key)
% Sensor arrays in msg_fwd's detection order: an experimental array, else
% front/back coils (3- then 2-axis), sensors, or 2-axis electrodes. Only the
% wanted arrays are returned.
    sens = {}; labels = {}; is_elec = false;
    if isfield(geoms, 'experimental_sensors')
        cand = {'experimental', geoms.experimental_sensors};
    else
        cand = cell(0, 2);
        for side = {'front', 'back'}
            for suffix = {'_coils_3axis', '_coils_2axis', '_sensors', '_sensors_2axis'}
                f = [side{1} suffix{1}];
                if isfield(geoms, f)
                    cand(end+1, :) = {side{1}, geoms.(f)}; %#ok<AGROW>
                    break
                end
            end
        end
    end
    for k = 1:size(cand, 1)
        if ~ismember(cand{k, 1}, wanted), continue; end
        labels{end+1} = cand{k, 1}; %#ok<AGROW>
        sens{end+1}   = cand{k, 2}; %#ok<AGROW>
    end
    if isempty(sens)
        warning('pt_compute_leadfields:noArray', ...
            'None of the arrays %s found in %s.', strjoin(wanted, ', '), key);
        return
    end
    s1 = sens{1};
    is_elec = (isfield(s1, 'elecpos') || isfield(s1, 'chanpos')) && ~isfield(s1, 'coilpos');
end


function F = default_fem(F)
    def = struct('binpath', '', 'maxvol_mm3', 500, 'opt_scale', 1, ...
        'reduce_torso', true, 'torso_keep', 0.5, 'seed', 1, ...
        'remesh_seeds', [], 'cond', [0.33, 0.33/40, 0.62, 0.05, 0.23]);
    for f = fieldnames(def)'
        if ~isfield(F, f{1}) || isempty(F.(f{1})) && ~strcmp(f{1}, 'remesh_seeds')
            F.(f{1}) = def.(f{1});
        end
    end
    exe = fullfile(F.binpath, 'bst_duneuro_meeg_win64.exe');
    if ~isfile(exe)
        error('pt_compute_leadfields:duneuro', ...
            ['DUNEuro binary not found: %s. Set duneuro_binpath in ' ...
             'config_pert_paths to the folder holding it.'], exe);
    end
    if isempty(which('fem_calc_fwds'))
        error('pt_compute_leadfields:femCalc', ...
            'fem_calc_fwds is not on the path (the FEM tools msg_fwd uses).');
    end
end

function v = getfield_or(s, f, v)
    if isfield(s, f), v = s.(f); end
end

function k = fem_mesh_key(key)
% A sensor shift changes no surface: it is solved on the unshifted anatomy's
% mesh, so its FEM difference holds no re-meshing
    k = regexprep(key, '_sensor_(original|bundle\d+_shift\d+)$', '_source_original');
end

function [tet, info] = fem_mesh(mesh_dir, geom_dir, mesh_key, seed, F, time_log, n_src)
% The cached volume mesh of one anatomy and seed, built on first use
    mfile = fullfile(mesh_dir, sprintf('volume_mesh_%s_seed%d.mat', mesh_key, seed));
    if isfile(mfile)
        m = load(mfile, 'tet', 'mesh_info');
        tet = m.tet; info = m.mesh_info;
        return
    end
    geoms = load(fullfile(geom_dir, ['geometries_' mesh_key '.mat']));
    [tet, info] = fem_build_volume_mesh(geoms, struct( ...
        'maxvol_mm3', F.maxvol_mm3, 'opt_scale', F.opt_scale, ...
        'reduce_torso', F.reduce_torso, 'torso_keep', F.torso_keep, ...
        'seed', seed, 'label', mesh_key, 'dump_dir', mesh_dir));
    log_compute_time(time_log, struct('model', mesh_key, 'method', 'fem', ...
        'stage', 'mesh', 'seconds', info.time_s.total, 'n_sources', n_src, ...
        'n_nodes', info.n_nodes, 'n_elements', info.n_tets, ...
        'note', sprintf('maxvol %g mm^3, seed %d', F.maxvol_mm3, seed)));
    mesh_info = info; %#ok<NASGU>
    save(mfile, 'tet', 'mesh_info', '-v7.3');
end

% pt_compute_metrics - All four comparison metrics for every perturbation
%
% For the active modality, compares every perturbed leadfield with its
% unperturbed reference and stores RE, r2, RDM and lnMAG at every source
% position, for every forward model method, dipole orientation and sensor
% axis set. Every later stage reads this one file, so all figures and tables
% are computed from identical numbers.
%
% USAGE:
%   pt_modality('set', 'msg');   % or 'esg' (run_perturbation_analysis does this)
%   pt_compute_metrics
%
% INPUT:
%   <forward_fields_base>/leadfields_organised.mat   from pt_load_leadfields
%
% OUTPUT:
%   <forward_fields_base>/pert_metrics.mat, holding struct PM:
%     .modality        'msg' | 'esg'
%     .dist_mm         [1 x n_src] cord distance of each kept source
%     .src_range       source indices kept (first and last are dropped)
%     .ori_labels      {'VD','RC','LR','ALL'}
%     .axcfg_ids       sensor axis per axis set; 0 = whole array
%     .axcfg_names     e.g. {'X','Y','Z','Array'}
%     .axcfg_slot      comparison slot of each axis set (0 = array)
%     .n_trunc         sensors per axis used
%     .pert.<type>     for type = source | sensor | cond:
%        .keys         {1 x n} perturbed geometry keys (without method prefix
%                      for source/sensor; full keys for cond)
%        .bundle_idx   [n x 1]
%        .shift_idx    [n x 1]
%        .magnitude    [n x 1] |shift| in mm, or mean % conductivity change
%        .mag_unit     'mm' | '%'
%        .methods      {1 x m} methods with results
%        .ref_key      reference geometry key
%        .M.<method>   struct with fields re, rsq, rdm, lnmag, each
%                      [n x n_src x n_axcfg x n_ori]; missing leadfields NaN
%
% METRICS:
%   From pt_metrics_block, identical in definition to msg_fwd's lf_metrics.
%   The unperturbed geometry is always the reference (the RE denominator).
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

config_pert;
pt_add_functions;

fprintf('pt_compute_metrics — modality %s\n\n', upper(active_modality));

lf_file = fullfile(forward_fields_base, 'leadfields_organised.mat');
if ~isfile(lf_file)
    error('pt_compute_metrics: %s not found. Run pt_load_leadfields first.', lf_file);
end
load(lf_file, 'leadfields', 'loaded_models');

mopts = metric_defaults();
oris  = metric_orientations;
n_ori = numel(oris);


% =========================================================================
% PERTURBATION DEFINITIONS
% =========================================================================

src_vec = vertcat(source_shift_vectors{:});
sen_vec = vertcat(sensor_shift_vectors{:});
sdims   = mods_cfg.(active_modality).sensor_shift_dims;

% Mean % conductivity change per realisation: replays the RNG stream of
% msg_fwd/run_conductivity_perturbation, so cond_seed and
% n_cond_compartments must match that run.
rng(cond_seed);
cond_pct = nan(n_cond_total, 1);
for b = 1:n_cond_bundles
    deltas = rand(n_cond_shifts, n_cond_compartments) * cond_bundle_pct(b);
    for s = 1:n_cond_shifts
        cond_pct((b-1)*n_cond_shifts + s) = mean(deltas(s, :)) * 100;
    end
end

P = struct();

P.source.keys       = sensitivity_keys;
P.source.ref_key    = sensitivity_ref_key;
P.source.bundle_idx = source_bundle_idx(:);
P.source.shift_idx  = source_shift_idx(:);
P.source.magnitude  = sqrt(sum(src_vec.^2, 2));
P.source.mag_unit   = 'mm';
P.source.shift_vec  = src_vec;
P.source.methods    = fwd_methods;

P.sensor.keys       = sensor_sensitivity_keys;
P.sensor.ref_key    = sensor_sensitivity_ref_key;
P.sensor.bundle_idx = sensor_sensitivity_bundle_idx(:);
P.sensor.shift_idx  = sensor_sensitivity_shift_idx(:);
P.sensor.magnitude  = sqrt(sum(sen_vec(:, sdims).^2, 2));
P.sensor.mag_unit   = 'mm';
P.sensor.shift_vec  = sen_vec .* ismember(1:3, sdims);
P.sensor.methods    = fwd_methods;

P.cond.keys         = cellfun(@(k) ['bem_cond_' cond_sensitivity_ref_key ...
                        regexprep(k, '^.*(_bundle\d+_shift\d+)$', '$1')], ...
                        cond_sensitivity_keys, 'UniformOutput', false);
P.cond.ref_key      = cond_sensitivity_ref_key;
P.cond.bundle_idx   = cond_bundle_idx(:);
P.cond.shift_idx    = cond_shift_idx(:);
P.cond.magnitude    = cond_pct;
P.cond.mag_unit     = '%';
P.cond.shift_vec    = [];
P.cond.methods      = {'bem'};

if ~have_bem_cond
    P = rmfield(P, 'cond');
end


% =========================================================================
% DIMENSIONS
% =========================================================================

% Any loaded reference defines the source grid and sensor count
ref_any = '';
for k = 1:numel(loaded_models)
    if endsWith(loaded_models{k}, [base_geom_name '_source_original'])
        ref_any = loaded_models{k};
        break
    end
end
if isempty(ref_any)
    error('pt_compute_metrics: no unperturbed reference leadfield was loaded.');
end

n_sources = leadfields.(ref_any).n_sources;
src_range = 2:(n_sources - 1);
n_src     = numel(src_range);
dist_mm   = src_range * src_spacing_mm;
n_axes    = leadfields.(ref_any).n_sensor_axes;

% Common sensor count across every loaded model
n_trunc = inf;
for k = 1:numel(loaded_models)
    n_trunc = min(n_trunc, leadfields.(loaded_models{k}).n_sensors_per_axis);
end

axcfg_ids   = [1:n_axes, 0];
axcfg_names = [axis_names(:)', {'Array'}];
axcfg_slot  = [mods_cfg.(active_modality).axis_slot, 0];
n_cfg       = numel(axcfg_ids);

fprintf('  %d sources (%d kept), %d sensor axes, %d sensors per axis\n', ...
    n_sources, n_src, n_axes, n_trunc);
fprintf('  Axis sets: %s\n', strjoin(axcfg_names, ', '));
fprintf('  Orientations: %s\n\n', strjoin(oris, ', '));


% =========================================================================
% COMPUTE
% =========================================================================

types = fieldnames(P)';
for t = 1:numel(types)
    ty = types{t};
    T  = P.(ty);
    n  = numel(T.keys);
    T.M = struct();
    done_methods = {};

    for m = 1:numel(T.methods)
        method  = T.methods{m};
        ref_key = [method '_' T.ref_key];
        if ~isfield(leadfields, ref_key)
            fprintf('  [%s/%s] reference %s not loaded — skipped\n', ty, method, ref_key);
            continue
        end
        if strcmp(ty, 'cond')
            full_keys = T.keys;
        else
            full_keys = cellfun(@(k) [method '_' k], T.keys, 'UniformOutput', false);
        end
        present = cellfun(@(k) isfield(leadfields, k), full_keys);
        fprintf('  [%s/%s] %d of %d perturbations loaded\n', ty, method, sum(present), n);
        if ~any(present), continue; end

        S = struct('re', nan(n, n_src, n_cfg, n_ori), 'rsq', nan(n, n_src, n_cfg, n_ori), ...
                   'rdm', nan(n, n_src, n_cfg, n_ori), 'lnmag', nan(n, n_src, n_cfg, n_ori));

        for c = 1:n_cfg
            for o = 1:n_ori
                A = pt_lf_matrix(leadfields.(ref_key), oris{o}, axcfg_ids(c), ...
                    n_trunc, src_range);
                for i = find(present(:))'
                    B  = pt_lf_matrix(leadfields.(full_keys{i}), oris{o}, ...
                        axcfg_ids(c), n_trunc, src_range);
                    Mi = pt_metrics_block(A, B, mopts);
                    S.re(i, :, c, o)    = Mi.re;
                    S.rsq(i, :, c, o)   = Mi.rsq;
                    S.rdm(i, :, c, o)   = Mi.rdm;
                    S.lnmag(i, :, c, o) = Mi.lnmag;
                end
            end
        end

        % A perturbation of this size cannot change the field magnitude by
        % an order of magnitude; if it appears to, the units are wrong.
        big = abs(median(S.lnmag(:, :, end, end), 2, 'omitnan')) > log(10);
        if any(big)
            warning('pt_compute_metrics:scale', ...
                ['[%s/%s] %d perturbation(s) differ from the reference by more ' ...
                 'than 10x in magnitude — a units mismatch, not a result. ' ...
                 'Check leadfield_scale_log.csv.'], ty, method, sum(big));
        end

        T.M.(method) = S;
        done_methods{end+1} = method; %#ok<SAGROW>
    end

    T.methods = done_methods;
    P.(ty) = T;
end

PM = struct();
PM.modality    = active_modality;
PM.display     = mods_cfg.(active_modality).display;
PM.unit        = mods_cfg.(active_modality).unit;
PM.dist_mm     = dist_mm;
PM.src_range   = src_range;
PM.ori_labels  = oris;
PM.axcfg_ids   = axcfg_ids;
PM.axcfg_names = axcfg_names;
PM.axcfg_slot  = axcfg_slot;
PM.radial_cfg  = radial_axis;
PM.n_trunc     = n_trunc;
PM.pert        = P;
PM.metric_opts = mopts;
PM.created     = char(datetime('now'));

outfile = fullfile(forward_fields_base, 'pert_metrics.mat');
save(outfile, 'PM', '-v7.3');
fprintf('\nSaved: %s\n', outfile);

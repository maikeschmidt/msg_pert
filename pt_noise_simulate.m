% pt_noise_simulate - Stage 5a: perturbations measured through realistic sensor noise
%
% Stages 1-4 compare noise-free forward fields. Here each field is "measured"
% by a real sensor system: an evoked response is generated through the
% (perturbed) leadfield, trial-averaged sensor noise is added, and the field
% pattern recovered from the noisy data is scored against the noise-free
% ORIGINAL field with the same four metrics as every other stage.
%
% THE MEASUREMENT MODEL
%   For one source, trial-averaged sensor data are Y = g w' + E, where g is
%   the leadfield of the geometry actually measured (perturbed or not), w the
%   evoked source waveform in nA*m, and E white noise with per-sample s.d.
%       sigma = density * sqrt(bandwidth) / sqrt(n_trials)
%   The field pattern is recovered by projecting the data onto the known
%   waveform (a matched filter):
%       g_hat = Y w / (w'w) = g + n,   n ~ N(0, sigma^2 / ||w||^2)
%   and compared with the unperturbed noise-free field a:
%       noise only             g = a         error is due to noise alone
%       perturbation + noise   g = perturbed error is model error AND noise
%   At zero noise g_hat = g, so level 0 reproduces stages 2-4 exactly.
%   Because the noise is zero-mean and independent of g,
%       E||g_hat - a||^2 = ||g - a||^2 + E||n||^2
%   i.e. squared RE from the perturbation and from the noise add.
%
% NOISE
%   Levels are multiples of each system's own floor (config_pert:
%   noise_factors, with 0 prepended). One set of noise draws per realisation
%   is shared across every geometry and every level (common random numbers),
%   so differences between geometries are not blurred by different draws.
%
% USAGE:
%   pt_noise_simulate            % then pt_noise_analyse
%
% INPUT:
%   leadfields_organised.mat and pert_metrics.mat for each modality used by
%   noise_systems (config_pert)
%
% OUTPUT (to <stage_results_dir>/5_noise/data/):
%   noise_<system>.mat   struct NS:
%     .system       the noise_systems entry
%     .levels       [1 x n_lev] noise factors (first = 0)
%     .sigma        [1 x n_lev] per-sample noise s.d. after trial averaging
%     .sd_ghat      [1 x n_lev] noise s.d. of the recovered field pattern
%     .snr          [n_lev x n_cfg x n_ori] median over cord of rms(a)/sd_ghat
%     .dist_mm, .axcfg_names, .ori_labels
%     .pert.<type>  .bundle_idx .magnitude .mag_unit, and
%        .cm.<metric>   [n+1 x n_lev x n_real x n_cfg x n_ori]
%                        cord median per realisation; row 1 = noise only
%        .src.<metric>  [n+1 x n_src x n_lev x n_cfg x n_ori]
%                        mean over realisations at each cord position
%     metrics: re, rsq, rdm, lnmag, absgain
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

pt_modality('clear');
config_pert;
pt_add_functions;

fprintf('Stage 5a — perturbations measured through sensor noise\n\n');

data_dir = fullfile(stage_results_dir, '5_noise', 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

% Source waveform and its norm (sets the noise on the recovered pattern)
t_ax  = 0 : 1/noise_fs : noise_duration - 1/noise_fs;
w     = noise_peak_nAm * sin(2*pi*noise_freq*(t_ax - noise_latency)) ...
        .* exp(-((t_ax - noise_latency).^2) / (2*noise_env_sd^2));
wnorm = norm(w);

levels = [0, noise_factors(:)'];
n_lev  = numel(levels);
n_real = noise_n_real;
mopts  = metric_defaults();
mnames = {'re', 'rsq', 'rdm', 'lnmag', 'absgain'};

fprintf('  Waveform: %g nA*m peak, %g Hz, ||w|| = %.3g nA*m*sqrt(samples)\n', ...
    noise_peak_nAm, noise_freq, wnorm);
fprintf('  %d trials, levels x floor: %s, %d realisations\n\n', noise_n_trials, ...
    mat2str(levels), n_real);

LFc = struct(); PMc = struct();   % caches, one per modality

for k = 1:numel(noise_systems)
    sys = noise_systems(k);
    md  = sys.modality;
    fprintf('[%d/%d] %s (%s, %s)\n', k, numel(noise_systems), sys.label, ...
        upper(md), pt_method_label(sys.method));

    if ~isfield(LFc, md)
        f1 = fullfile(mods_cfg.(md).forward_fields_base, 'leadfields_organised.mat');
        f2 = fullfile(mods_cfg.(md).forward_fields_base, 'pert_metrics.mat');
        if ~isfile(f1) || ~isfile(f2)
            warning('  %s or %s missing — %s skipped.', f1, f2, sys.label);
            continue
        end
        tmp = load(f1, 'leadfields'); LFc.(md) = tmp.leadfields;
        tmp = load(f2, 'PM');         PMc.(md) = tmp.PM;
    end
    LF = LFc.(md); PM = PMc.(md);

    bw_eff  = min(sys.bandwidth_hz, noise_fs / 2);
    sigma   = sys.density * levels * sqrt(bw_eff) / sqrt(noise_n_trials);
    sd_ghat = sigma / wnorm;
    fprintf('  effective bandwidth %g Hz; sigma at 1x = %.3g; pattern noise s.d. at 1x = %.3g %s\n', ...
        bw_eff, sys.density * sqrt(bw_eff) / sqrt(noise_n_trials), ...
        sys.density * sqrt(bw_eff) / sqrt(noise_n_trials) / wnorm, mods_cfg.(md).unit);

    n_cfg = numel(PM.axcfg_ids);
    n_ori = numel(PM.ori_labels);
    n_src = numel(PM.src_range);

    NS = struct();
    NS.system = sys; NS.levels = levels; NS.sigma = sigma; NS.sd_ghat = sd_ghat;
    NS.dist_mm = PM.dist_mm; NS.axcfg_names = PM.axcfg_names;
    NS.ori_labels = PM.ori_labels; NS.modality = md; NS.n_real = n_real;
    NS.snr = nan(n_lev, n_cfg, n_ori);
    NS.pert = struct();

    types = pert_types(cellfun(@(ty) isfield(PM.pert, ty) && ...
        isfield(PM.pert.(ty).M, sys.method), pert_types));

    for t = 1:numel(types)
        ty = types{t}; P = PM.pert.(ty);
        n  = numel(P.keys);
        NS.pert.(ty).bundle_idx = P.bundle_idx;
        NS.pert.(ty).magnitude  = P.magnitude;
        NS.pert.(ty).mag_unit   = P.mag_unit;
        for f = mnames
            NS.pert.(ty).cm.(f{1})  = nan(n + 1, n_lev, n_real, n_cfg, n_ori);
            NS.pert.(ty).src.(f{1}) = nan(n + 1, n_src, n_lev, n_cfg, n_ori);
        end
    end

    for c = 1:n_cfg
        for o = 1:n_ori
            % Common random numbers: one draw per realisation, reused for
            % every geometry and level of this system/axis set/orientation
            rng(noise_seed + 7919*k + 101*c + o);
            Z = [];

            for t = 1:numel(types)
                ty = types{t}; P = PM.pert.(ty);
                ref_key = [sys.method '_' P.ref_key];
                if ~isfield(LF, ref_key), continue; end
                A = pt_lf_matrix(LF.(ref_key), PM.ori_labels{o}, PM.axcfg_ids(c), ...
                    PM.n_trunc, PM.src_range);
                if isempty(Z), Z = randn([size(A), n_real]); end
                if t == 1
                    NS.snr(:, c, o) = median(sqrt(mean(A.^2, 1))) ./ sd_ghat(:);
                end

                if strcmp(ty, 'cond')
                    keys = P.keys;
                else
                    keys = cellfun(@(x) [sys.method '_' x], P.keys, 'UniformOutput', false);
                end

                for i = 0:numel(keys)
                    if i == 0
                        G = A;
                    elseif isfield(LF, keys{i})
                        G = pt_lf_matrix(LF.(keys{i}), PM.ori_labels{o}, ...
                            PM.axcfg_ids(c), PM.n_trunc, PM.src_range);
                    else
                        continue
                    end
                    acc = struct('re', zeros(1, n_src), 'rsq', zeros(1, n_src), ...
                        'rdm', zeros(1, n_src), 'lnmag', zeros(1, n_src), ...
                        'absgain', zeros(1, n_src));
                    for L = 1:n_lev
                        for f = mnames, acc.(f{1})(:) = 0; end
                        for r = 1:n_real
                            M = pt_metrics_block(A, G + sd_ghat(L) * Z(:, :, r), mopts);
                            M.absgain = abs(exp(M.lnmag) - 1) * 100;
                            for f = mnames
                                NS.pert.(ty).cm.(f{1})(i+1, L, r, c, o) = ...
                                    median(M.(f{1}), 'omitnan');
                                acc.(f{1}) = acc.(f{1}) + M.(f{1});
                            end
                            if sd_ghat(L) == 0, break; end   % level 0 is deterministic
                        end
                        n_used = r;
                        for f = mnames
                            NS.pert.(ty).src.(f{1})(i+1, :, L, c, o) = acc.(f{1}) / n_used;
                            if n_used < n_real   % fill the deterministic level
                                NS.pert.(ty).cm.(f{1})(i+1, L, :, c, o) = ...
                                    NS.pert.(ty).cm.(f{1})(i+1, L, 1, c, o);
                            end
                        end
                    end
                end
            end
        end
        fprintf('  axis set %s done\n', PM.axcfg_names{c});
    end

    outfile = fullfile(data_dir, sprintf('noise_%s.mat', sys.short));
    save(outfile, 'NS', '-v7.3');
    fprintf('  saved %s\n\n', outfile);
end

fprintf('Stage 5a complete. Next: pt_noise_analyse\n');

% sim_simulate_noise - Simulate sensor data with realistic noise and score it
%                      against the perfect (noise-free) forward field
%
% Steps 2-4 of the realistic-measurement analysis:
%   2. Simulate a 1 nA*m, 90 Hz sinusoidal source at EVERY source point on the
%      cord, projected through the perfect forward model.
%   3. Add white sensor noise for three real systems (SQUID MSG, OP-MSG, ESG),
%      sweeping each system's noise floor up and down by the same multiplicative
%      steps around its published baseline.
%   4. Score the noisy data against the noise-free data with r-squared, per
%      source point and per noise level.
%
% USAGE:
%   sim_simulate_noise
%
% OUTPUT (to <sim_out_dir>/sim_noise_rsq.mat):
%   rsq_mean   [n_sys x n_ori x n_src x n_lev]  mean r^2 over realisations
%   rsq_sd     [n_sys x n_ori x n_src x n_lev]  s.d. across realisations
%   snr        [n_sys x n_ori x n_src x n_lev]  amplitude SNR, rms(signal)/sigma
%   sigma_abs  [n_sys x n_lev]  time-domain noise s.d. in each system's units
%   plus the source axis, labels, and a copy of the sweep settings
%
% NOISE MODEL:
%   Each system's noise floor is a white spectral density (units/sqrt(Hz)).
%   Broadband white noise sampled at sim_fs occupies a bandwidth of sim_fs/2,
%   so the time-domain standard deviation is
%       sigma = density * sqrt(sim_fs / 2)
%   No band-pass filter is applied around the 90 Hz signal. Narrow-band
%   filtering would shrink the noise by sqrt(bandwidth ratio) and inflate every
%   r-squared value, so leaving it broadband keeps this a conservative bound
%   and avoids baking a filter choice into the result.
%
% WHY THE r-SQUARED IS COMPUTED IN CLOSED FORM:
%   The simulated signal is rank-1: sensor data = leadfield vector x waveform.
%   Materialising the full [n_channels x n_time] matrix for every
%   (system, orientation, source, level, realisation) combination would mean
%   tens of thousands of large matrix operations. Because the signal is rank-1,
%   the Pearson correlation between the clean and noisy data reduces exactly to
%   a few scalar sums, which collapse the inner loop to O(n_channels). This is
%   an exact identity, not an approximation — see the derivation inline below.
%
% COMMON RANDOM NUMBERS:
%   One unit-normal noise draw is reused across all noise levels, orientations,
%   and any systems sharing a forward model, rescaled by sigma. Systems that
%   measure the same field (SQUID and OP-MSG) therefore see the SAME noise
%   pattern at different amplitudes, which isolates the effect of the noise
%   floor from realisation-to-realisation luck.
%
% DEPENDENCIES:
%   config_sim, pt_add_functions, sim_load_leadfield()
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

config_sim;
pt_add_functions;

fprintf('sim_simulate_noise\n');
fprintf('  Source:  %.1f nA*m, %.0f Hz sine, %.2f s @ %d Hz\n', ...
    sim_dipole_nAm, sim_freq, sim_duration, sim_fs);
fprintf('  Array:   %s\n', sim_array);
fprintf('  Levels:  %s (x baseline)\n', mat2str(sim_noise_factors));
fprintf('  Realisations per level: %d\n\n', sim_n_realisations);


% =========================================================================
% LOAD THE FORWARD MODELS USED BY THE SYSTEMS
% =========================================================================
% Only load models that a system actually references, and only once each.

n_sys   = numel(sim_systems);
n_ori   = numel(sim_orientations);
n_lev   = numel(sim_noise_factors);
n_model = numel(sim_models);

lf_by_model = cell(n_model, 1);
models_used = unique([sim_systems.model]);

for m = models_used
    fprintf('  Loading %-20s (%s array) ... ', sim_models(m).label, sim_array);
    lf_by_model{m} = sim_load_leadfield(sim_models(m).(sim_array), ...
        sim_models(m).var, sim_models(m).scale, sim_models(m).is_meg, ...
        sim_models(m).n_axes);
    fprintf('%d sources, %d channels total\n', ...
        lf_by_model{m}.n_sources, ...
        lf_by_model{m}.n_sensor_axes * lf_by_model{m}.n_sensors_per_axis);
end

% All models must cover the same source grid, or "the same source point" means
% different things in different rows of the comparison.
n_src_each = cellfun(@(l) l.n_sources, lf_by_model(models_used));
if numel(unique(n_src_each)) > 1
    error(['Models disagree on source count (%s). The MSG and ESG forward\n' ...
           'models must be built on the same source grid to be comparable.'], ...
           mat2str(n_src_each(:)'));
end
n_src = n_src_each(1);

src_mm = (0:n_src-1) * src_spacing_mm;
fprintf('\n  Source grid: %d points, 0 to %.0f mm\n\n', n_src, src_mm(end));


% =========================================================================
% FLATTEN LEADFIELDS: all sensor axes concatenated into one channel vector
% =========================================================================
% A real recording reads every channel at once — a triaxial MSG array measures
% all three axes simultaneously, and ESG reads every electrode. So the noise is
% scored against the FULL channel set, not one sensor axis at a time.

lfvec = cell(n_model, n_ori);   % {model, ori} -> [n_channels x n_src]

for m = models_used
    lf = lf_by_model{m};
    for o = 1:n_ori
        ori = sim_orientations{o};
        n_ch = lf.n_sensor_axes * lf.n_sensors_per_axis;
        V = zeros(n_ch, n_src);
        for s = 1:n_src
            V(:, s) = vertcat(lf.(ori){:, s});   % stack sensor axes
        end
        lfvec{m, o} = V;
    end
end


% =========================================================================
% NOISE AMPLITUDES
% =========================================================================
% sigma = spectral density * sqrt(bandwidth), bandwidth = Nyquist = fs/2

bandwidth = sim_fs / 2;
sigma_abs = zeros(n_sys, n_lev);   % time-domain s.d., in each system's own units

fprintf('  Noise floors (time-domain s.d. over %.0f Hz bandwidth):\n', bandwidth);
for k = 1:n_sys
    sigma_abs(k, :) = sim_systems(k).noise_baseline * sim_noise_factors ...
                      * sqrt(bandwidth);
    fprintf('    %-10s baseline %6.2f %s  ->  sigma at 1x = %8.2f\n', ...
        sim_systems(k).label, sim_systems(k).noise_baseline, ...
        sim_systems(k).noise_unit, sigma_abs(k, sim_noise_factors == 1));
end
fprintf('\n');


% =========================================================================
% SIMULATE
% =========================================================================
% Closed-form Pearson r^2 for a rank-1 signal.
%
% Let the clean sensor data be the outer product  x = g * w'  where
%   g = leadfield vector  [n_ch x 1]     (fT/nAm or uV/nAm)
%   w = source waveform   [1 x n_time]   (nA*m)
% and let the noisy data be  y = x + sigma*Z,  Z ~ N(0,1) elementwise.
%
% Pearson r over all n = n_ch*n_time elements needs only these sums:
%   Sx  = sum(x)    = sum(g) * sum(w)
%   Sxx = sum(x.^2) = sum(g.^2) * sum(w.^2)          (outer-product identity)
%   Sxz = sum(x.*z) = sigma * g' * (Z * w')          (contract time, then space)
%   Sz  = sigma * sum(Z(:))
%   Szz = sigma^2 * sum(Z(:).^2)
% from which
%   Sy  = Sx + Sz
%   Sxy = Sxx + Sxz
%   Syy = Sxx + 2*Sxz + Szz
%   r   = (n*Sxy - Sx*Sy) / sqrt((n*Sxx - Sx^2) * (n*Syy - Sy^2))
%
% The only cost that scales with n_time is the single contraction u = Z*w',
% done once per (model, source, realisation) and reused for every orientation,
% every noise level, and every system sharing that model.

w      = sim_waveform(:);          % [n_time x 1]
n_time = numel(w);
Sw     = sum(w);
Sww    = sum(w.^2);

rsq_all = nan(n_sys, n_ori, n_src, n_lev, sim_n_realisations);
snr     = nan(n_sys, n_ori, n_src, n_lev);

rng(sim_noise_seed);

fprintf('  Simulating %d sources x %d levels x %d realisations...\n', ...
    n_src, n_lev, sim_n_realisations);
t_start = tic;

for m = models_used

    n_ch = size(lfvec{m, 1}, 1);
    n_el = n_ch * n_time;                       % elements in the data matrix
    sys_of_model = find([sim_systems.model] == m);

    % Precompute the signal-only sums per orientation and source
    Sg   = zeros(n_ori, n_src);
    Sgg  = zeros(n_ori, n_src);
    for o = 1:n_ori
        Sg(o, :)  = sum(lfvec{m, o},    1);
        Sgg(o, :) = sum(lfvec{m, o}.^2, 1);
    end

    for s = 1:n_src
        for r = 1:sim_n_realisations

            % One unit-normal draw, reused across orientations, levels, and
            % every system that measures this model (common random numbers).
            Z  = randn(n_ch, n_time);
            u  = Z * w;              % [n_ch x 1]  — the only O(n_time) step
            SZ  = sum(Z(:));
            SZZ = sum(Z(:).^2);

            for o = 1:n_ori
                g   = lfvec{m, o}(:, s);
                Sx  = Sg(o, s)  * Sw;
                Sxx = Sgg(o, s) * Sww;
                gu  = g' * u;        % = sum over elements of (g*w') .* Z

                % Degenerate source (leadfield identically zero): r^2 undefined
                if Sxx <= 0
                    continue
                end

                for k = sys_of_model
                    for L = 1:n_lev
                        sg = sigma_abs(k, L);

                        Sxz = sg   * gu;
                        Sz  = sg   * SZ;
                        Szz = sg^2 * SZZ;

                        Sy  = Sx  + Sz;
                        Sxy = Sxx + Sxz;
                        Syy = Sxx + 2*Sxz + Szz;

                        den = (n_el*Sxx - Sx^2) * (n_el*Syy - Sy^2);
                        if den <= 0
                            continue
                        end
                        rr = (n_el*Sxy - Sx*Sy) / sqrt(den);

                        rsq_all(k, o, s, L, r) = rr^2;

                        if r == 1
                            % Amplitude SNR: rms of the clean signal over the
                            % noise s.d. Same for every realisation.
                            snr(k, o, s, L) = sqrt(Sxx / n_el) / sg;
                        end
                    end
                end
            end
        end
    end

    fprintf('    %-20s done (%d channels)\n', sim_models(m).label, n_ch);
end

rsq_mean = mean(rsq_all, 5, 'omitnan');
rsq_sd   = std(rsq_all, 0, 5, 'omitnan');

fprintf('  Elapsed: %.1f s\n\n', toc(t_start));


% =========================================================================
% SAVE
% =========================================================================

sys_labels   = {sim_systems.label};
sys_shorts   = {sim_systems.short};
sys_colors   = vertcat(sim_systems.color);
sys_baseline = [sim_systems.noise_baseline];
sys_units    = {sim_systems.noise_unit};

outfile = fullfile(sim_out_dir, 'sim_noise_rsq.mat');
save(outfile, 'rsq_mean', 'rsq_sd', 'snr', 'sigma_abs', ...
    'sim_noise_factors', 'sim_n_realisations', 'sim_noise_seed', ...
    'sim_orientations', 'sim_ori_display', 'sim_array', ...
    'sim_fs', 'sim_freq', 'sim_duration', 'sim_dipole_nAm', ...
    'src_mm', 'n_src', 'src_spacing_mm', ...
    'sys_labels', 'sys_shorts', 'sys_colors', 'sys_baseline', 'sys_units', ...
    '-v7.3');

fprintf('Saved: %s\n', outfile);


% =========================================================================
% QUICK SUMMARY AT BASELINE
% =========================================================================

base_L = find(sim_noise_factors == 1, 1);
if ~isempty(base_L)
    fprintf('\n  Mean r^2 across the whole cord, at each system''s baseline noise:\n');
    for o = 1:n_ori
        fprintf('    %-16s', sim_ori_display{o});
        for k = 1:n_sys
            fprintf('  %-10s %.3f', sim_systems(k).label, ...
                mean(squeeze(rsq_mean(k, o, :, base_L)), 'omitnan'));
        end
        fprintf('\n');
    end
end

fprintf('\nNext: run sim_plot_noise_curves\n');

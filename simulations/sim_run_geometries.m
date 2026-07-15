% sim_run_geometries - Evoked-response + sensor-noise simulation, looped over
%                      every geometry variant in config_sim
%
% Runs the (working) single-geometry noise simulation once per entry in
% sim_geometries, into its own output subfolder. For each variant and each
% applicable sensor system it simulates the evoked response through that
% geometry's leadfield, adds trial-averaged noise across the sweep, and scores
% r^2 against the noise-free field per source and noise level.
%
% This is the geometry-looping generalisation of the original steps 1-4: point
% it at a list of leadfield files (the original plus one representative shift per
% source / sensor / conductivity bundle) and it produces a comparable noise
% curve for each.
%
% USAGE:
%   sim_run_geometries          % then sim_plot_comparison
%
% OUTPUT (per variant, to <sim_out_dir>/<variant.name>/):
%   sim_noise_<system>.mat   r^2 and SNR for that system on that geometry
% and a combined index at <sim_out_dir>/sim_geometry_index.mat listing every
% (variant, system) that produced output, for sim_plot_comparison.
%
% MODALITY / SKIPPING:
%   A system runs on a variant only if that variant's leadfield exists for the
%   system's model. Source, sensor, and conductivity variants all exist for both
%   MSG and ESG (conductivity under each modality's own cond_root). Any missing
%   (variant, system) combination is skipped automatically — no special-casing.
%
% DEPENDENCIES:
%   config_sim, sim_lf_path(), sim_load_leadfield(), sim_evoked_noise_rsq()
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

fprintf('sim_run_geometries\n');
fprintf('  %d geometry variants x %d systems, array = %s\n', ...
    numel(sim_geometries), numel(sim_systems), sim_array);
fprintf('  Evoked burst %g nA*m @ %g Hz, %d trials, %d noise levels\n\n', ...
    sim_dipole_nAm, sim_freq, sim_n_trials, numel(sim_noise_factors));


% -------------------------------------------------------------------------
% Noise s.d. per system per level (bandwidth + trial averaging)
% -------------------------------------------------------------------------
nyquist    = sim_fs / 2;
trial_gain = sqrt(sim_n_trials);
n_lev      = numel(sim_noise_factors);
n_sys      = numel(sim_systems);
n_ori      = numel(sim_orientations);

sigma_abs = zeros(n_sys, n_lev);
bw_eff    = zeros(1, n_sys);
for k = 1:n_sys
    bw_eff(k)      = min(sim_systems(k).bandwidth_hz, nyquist);
    sigma_abs(k,:) = sim_systems(k).noise_baseline * sim_noise_factors ...
                     * sqrt(bw_eff(k)) / trial_gain;
end

w = sim_waveform(:);

% Index of everything successfully simulated, for the comparison plot
index = struct('variant', {}, 'name', {}, 'group', {}, 'bundle', {}, ...
               'system', {}, 'sys_short', {}, 'file', {});


% -------------------------------------------------------------------------
% Loop over geometry variants
% -------------------------------------------------------------------------
for gi = 1:numel(sim_geometries)
    geo = sim_geometries(gi);

    out_sub = fullfile(sim_out_dir, geo.name);
    if ~exist(out_sub, 'dir'); mkdir(out_sub); end

    fprintf('[%d/%d] %-16s (%s, %s)\n', gi, numel(sim_geometries), ...
        geo.name, geo.group, geo.short);

    for k = 1:n_sys
        model = sim_models(sim_systems(k).model);

        lf_path = sim_lf_path(model, geo, sim_array);
        if ~isfile(lf_path)
            fprintf('    %-10s  no leadfield (%s) — skipped\n', ...
                sim_systems(k).label, geo.kind);
            continue
        end

        lf = sim_load_leadfield(lf_path, model.var, model.scale, ...
            model.is_meg, model.n_axes);
        n_src = lf.n_sources;

        % Flatten to full channel vector per source, per orientation, and score
        rsq_mean = nan(n_ori, n_src, n_lev);
        rsq_sd   = nan(n_ori, n_src, n_lev);
        snr      = nan(n_ori, n_src, n_lev);

        for o = 1:n_ori
            ori  = sim_orientations{o};
            V    = zeros(lf.n_sensor_axes * lf.n_sensors_per_axis, n_src);
            for s = 1:n_src
                V(:, s) = vertcat(lf.(ori){:, s});
            end
            [rm, rs, sn] = sim_evoked_noise_rsq(V, w, sigma_abs(k, :), ...
                sim_n_realisations);
            rsq_mean(o, :, :) = rm;
            rsq_sd(o, :, :)   = rs;
            snr(o, :, :)      = sn;
        end

        src_mm = (0:n_src-1) * src_spacing_mm;
        outfile = fullfile(out_sub, sprintf('sim_noise_%s.mat', sim_systems(k).short));
        save(outfile, 'rsq_mean', 'rsq_sd', 'snr', 'src_mm', 'n_src', ...
            'sim_noise_factors', 'sigma_abs', 'bw_eff', 'sim_orientations', ...
            'sim_ori_display', 'sim_array', 'sim_n_trials', 'sim_dipole_nAm', ...
            'sim_freq', 'geo', '-v7.3');

        % Peak signal vs baseline noise, so a scale mismatch is visible per run
        [~, base_col] = min(abs(sim_noise_factors - 1));
        peak_g = 0;
        for o = 1:n_ori
            for s = 1:n_src
                peak_g = max(peak_g, max(abs(vertcat(lf.(sim_orientations{o}){:, s}))));
            end
        end
        fprintf('    %-10s  peak|g|=%.3g  sigma@1x=%.3g  meanR^2@1x=%.3f\n', ...
            sim_systems(k).label, peak_g, sigma_abs(k, base_col), ...
            mean(rsq_mean(:, :, base_col), 'all', 'omitnan'));

        index(end+1) = struct('variant', gi, 'name', geo.name, ...
            'group', geo.group, 'bundle', geo.bundle, ...
            'system', sim_systems(k).label, 'sys_short', sim_systems(k).short, ...
            'file', outfile); %#ok<SAGROW>
    end
end

if isempty(index)
    error(['No leadfields were found for any variant. Check the roots and ' ...
           'geometry stems in config_sim.']);
end

idxfile = fullfile(sim_out_dir, 'sim_geometry_index.mat');
save(idxfile, 'index', 'sim_noise_factors', 'sim_orientations', ...
    'sim_ori_display', 'sim_bundle_colors', 'sim_array', 'sim_n_trials', ...
    'sim_dipole_nAm', 'sim_freq', '-v7.3');

fprintf('\nSaved %d (variant x system) results.\nIndex: %s\n', numel(index), idxfile);
fprintf('Next: run sim_plot_comparison\n');

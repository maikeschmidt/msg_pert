% sim_plot_noise_topoplot - Topoplots of a chosen source at a chosen noise level
%
% Step 5 of the realistic-measurement analysis. Picks one source point and one
% noise level and shows, side by side, what each system would ACTUALLY measure
% versus the perfect forward field.
%
% This is the figure that makes the r-squared curves concrete: it shows whether
% a given r^2 still corresponds to a recognisable field pattern or to a noise
% blob, which the summary statistic alone cannot tell you.
%
% USAGE:
%   sim_plot_noise_topoplot
%
% CONFIGURATION (in config_sim):
%   sim_focus_src_mm       - cord distance (mm) of the source to show
%   sim_focus_noise_factor - noise level, as a multiple of system baseline
%   sim_array              - 'front' or 'back'
%
% OUTPUT (to <sim_save_dir>/noise_topoplots/):
%   noisy_topoplot_<system>_src<MM>mm_<factor>x_<array>.png / .fig
%     One figure per system. Rows = [noise-free; noisy], columns = orientation.
%     One figure per sensor axis.
%
% WHAT IS PLOTTED:
%   The sensor data at the time sample where the source waveform peaks, i.e.
%   the instant of maximum signal. Noise is a single random draw at the
%   configured level — this is one plausible measurement, not an average, which
%   is the point: it shows the noise a real recording would actually contain.
%
% COLOUR LIMITS:
%   Shared between the clean and noisy panels of the same orientation, so the
%   noisy panel is not silently rescaled to look cleaner than it is.
%
% DEPENDENCIES:
%   config_sim, pt_add_functions
%   sim_load_leadfield(), sim_sensor_positions()  — msg_pert/functions/
%   plot_topoplot_publication()                   — msg_fwd/functions/
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

save_dir = fullfile(sim_save_dir, 'noise_topoplots');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

n_sys = numel(sim_systems);
n_ori = numel(sim_orientations);

% Time sample at the peak of the source waveform — the best-case instant
[~, t_peak] = max(abs(sim_waveform));
w_peak      = sim_waveform(t_peak);

% Noise s.d. at the configured level (see sim_simulate_noise for the derivation).
% This must include the sqrt(n_trials) reduction from trial averaging, or the
% topoplot would show single-trial noise while the r-squared curves describe
% trial-averaged data — the two figures would silently disagree.
bandwidth  = sim_fs / 2;
trial_gain = sqrt(sim_n_trials);

fprintf('sim_plot_noise_topoplot\n');
fprintf('  Source: %.0f mm    Noise: %gx baseline    Array: %s\n', ...
    sim_focus_src_mm, sim_focus_noise_factor, sim_array);
fprintf('  Plotting at t = %.1f ms (waveform peak)\n\n', ...
    1000 * sim_time(t_peak));

rng(sim_noise_seed);

for k = 1:n_sys

    m  = sim_systems(k).model;
    md = sim_models(m);

    fprintf('  %s (%s)\n', sim_systems(k).label, md.label);

    lf  = sim_load_leadfield(md.(sim_array), md.var, md.scale, md.is_meg, md.n_axes);
    pos = sim_sensor_positions(md.geom_file, sim_array, md.is_meg, lf.n_sensor_axes);

    src_idx = round(sim_focus_src_mm / src_spacing_mm) + 1;
    if src_idx > lf.n_sources
        error(['Source index %d exceeds the %d sources in %s.\n' ...
               'Check sim_focus_src_mm (%.0f mm) and src_spacing_mm (%.1f mm).'], ...
               src_idx, lf.n_sources, md.label, sim_focus_src_mm, src_spacing_mm);
    end

    sigma = sim_systems(k).noise_baseline * sim_focus_noise_factor ...
            * sqrt(bandwidth) / trial_gain;

    if md.is_meg
        unit_str = 'fT';
    else
        unit_str = '\muV';
    end

    if numel(md.axis_names) < lf.n_sensor_axes
        error(['Model "%s" has %d sensor axes but only %d entries in .axis_names.\n' ...
               'Set sim_models(%d).axis_names in config_sim.m.'], ...
               md.label, lf.n_sensor_axes, numel(md.axis_names), m);
    end

    % Loops over the axes this model ACTUALLY has — 3 for triaxial MSG,
    % 2 for the ESG tangential/radial electrode sets.
    for ax = 1:lf.n_sensor_axes

        n_sens = lf.n_sensors_per_axis;

        fig = figure('Color', 'w', 'Units', 'inches', ...
            'Position', [1, 1, 3*n_ori + 1.5, 6]);
        tl = tiledlayout(2, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');

        for o = 1:n_ori
            ori = sim_orientations{o};

            % Clean sensor data at the waveform peak:
            %   leadfield (unit/nAm) x source amplitude (nAm) = unit
            clean = lf.(ori){ax, src_idx} * w_peak;

            % One plausible noisy measurement at this instant
            noisy = clean + sigma * randn(n_sens, 1);

            % Shared colour limit so the noisy panel is not flattered by a
            % rescaled colour axis
            cmax = max(max(abs(clean)), max(abs(noisy)));
            if cmax == 0; cmax = 1; end
            clim = [-cmax, cmax];

            nexttile(tl, o);
            plot_topoplot_publication(pos{ax}, clean, clim, md.is_meg);
            title(sim_ori_display{o}, 'FontWeight', 'bold', 'FontSize', 11);
            if o == 1
                text(-0.12, 0.5, sprintf('Noise-free\n(%s)', unit_str), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Rotation', 90, ...
                    'FontWeight', 'bold', 'FontSize', 10);
            end

            nexttile(tl, n_ori + o);
            plot_topoplot_publication(pos{ax}, noisy, clim, md.is_meg);
            if o == 1
                text(-0.12, 0.5, sprintf('Measured\n(%s)', unit_str), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Rotation', 90, ...
                    'FontWeight', 'bold', 'FontSize', 10);
            end
        end

        title(tl, sprintf('%s  |  source %.0f mm  |  %s  |  %s array', ...
            sim_systems(k).label, sim_focus_src_mm, ...
            md.axis_names{ax}, sim_array), ...
            'FontSize', 13, 'FontWeight', 'bold');
        subtitle(tl, sprintf(['noise %gx baseline (%g %s) over %g Hz, averaged over ' ...
                              '%d trials  ->  \\sigma = %.3g %s;  evoked burst ' ...
                              '%g nA\\cdotm @ %g Hz, peak %.0f ms'], ...
            sim_focus_noise_factor, sim_systems(k).noise_baseline, ...
            sim_systems(k).noise_unit, bandwidth, sim_n_trials, ...
            sigma, unit_str, sim_dipole_nAm, sim_freq, sim_evoked_latency*1000), ...
            'FontSize', 9, 'Color', [0.4 0.4 0.4]);

        ax_slug = lower(regexprep(md.axis_names{ax}, '[^A-Za-z0-9]', ''));
        fname = sprintf('noisy_topoplot_%s_src%03dmm_%gx_%s_%s', ...
            sim_systems(k).short, round(sim_focus_src_mm), ...
            sim_focus_noise_factor, ax_slug, sim_array);
        exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 600);
        saveas(fig, fullfile(save_dir, [fname '.fig']));
        close(fig);

        fprintf('    Saved: %s\n', fname);
    end
end

fprintf('\nsim_plot_noise_topoplot complete.\nFigures: %s\n', save_dir);

% pt_cord_activity - Stage 7: N11 volley and N13 as a recording, BEM vs FEM
%
% Stage 5 measures one source at a time. Here the spinal response to median
% nerve stimulation is simulated as a recording: every active source summed,
% over time, through each forward model and sensor system:
%
%   N11  the afferent volley entering the cord at C7 and ascending in the
%        dorsal columns, a travelling quadrupole (leading and trailing
%        longitudinal dipoles of opposite sign)
%   N13  the stationary dorsal-horn response at C7, transverse, swept from a
%        single source to a wide patch along the cord (cord_n13_fwhm_mm)
%
% and asks, as the N13 spreads:
%
%   A0. Spread effect     how much does spreading change the recording at
%                         all, within one model? (bounds every effect below)
%   A. Model divergence   how different are the BEM-, FEM- and Biot-Savart-
%                         predicted recordings? (MRI-derived BEM reference)
%   B. Detectability      how well is each component recovered through each
%                         system's noise, and up to what noise level is it
%                         still detected, per forward model?
%   C. Perturbations      how much do the source, sensor and conductivity
%                         perturbations change the recording, per model,
%                         and do BEM and FEM differ in that sensitivity?
%
% A simulation cannot say which model is right: each recording is generated
% by the model being tested. What it can say is whether BEM and FEM diverge
% more for spread sources, and whether spread sources are more or less
% detectable and robust under each. Which one fits measured data better is a
% question for msg_analysis.
%
% SCORING
%   Every comparison uses the four msg_fwd metrics on the recording within a
%   component window (sensors x time, flattened): the N11 window spans the
%   volley's travel, the N13 window its peak +- 2 s.d. Detection: a
%   component counts as detected in one noise realisation when the mean
%   square of the measured data in its window exceeds the noise_detect_pct
%   percentile, over realisations, of the same statistic for noise alone in
%   that window (not the pre-stimulus baseline: the zero-phase bandwidth
%   filter smears the N11 back into it for narrow-band systems). The
%   critical noise
%   level is where the detection rate falls below noise_detect_rate
%   (pt_critical_level, as in stage 5), over stage 7's own wider range of
%   noise levels (cord_noise_factors: the pooled window test stays sensitive
%   far above the floor, so stage 5's range would leave every component
%   detected at its top level). Noise draws are shared across
%   extents and across BEM, FEM and Biot-Savart for the same sensors.
%
% USAGE:
%   pt_cord_activity
%
% INPUT:
%   leadfields_organised.mat and pert_metrics.mat for each modality
%
% OUTPUT (to <stage_results_dir>/7_cord_activity/):
%   source_model.png                  the simulated activity, cord x time
%   spread_effect.png, cord_spread_effect.*   A0
%   recordings_<modality>.png         noise-free recordings per model, and
%                                     one noisy realisation per system
%   divergence.png, cord_divergence.* A
%   snr_by_system.png                 B: SNR at 1x per system and model
%   detection_vs_noise.png            B: detection rate against noise level
%   detectability_snr.png, cord_snr.* B: exact SNR vs N13 extent, relative to
%                                     the focal N13 and to the BEM
%   detectability_<window>.png        B: critical noise level vs N13 extent
%   cord_detectability.csv            B: every system, extent, level, window
%   cord_critical_levels.*            B
%   cord_bem_vs_fem_noise.*           B: FEM minus BEM recovery error, paired
%   perturbation_<modality>.png       C: perturbation effect vs N13 extent
%   cord_perturbation.csv             C: every realisation summarised
%   cord_perturbation_bem_vs_fem.*    C: FEM minus BEM, paired by realisation
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
rng(stats_seed);

fprintf('Stage 7 — cord activity: N11 volley and N13 extent, BEM vs FEM\n\n');

out_dir = fullfile(stage_results_dir, '7_cord_activity');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

mopts  = metric_defaults();
fwhms  = cord_n13_fwhm_mm(:)';
n_ext  = numel(fwhms);
e_head = find(fwhms == cord_n13_head_fwhm, 1);
if isempty(e_head), e_head = min(2, n_ext); end
levels = [0, cord_noise_factors(:)'];   % stage 7's own, wider range (config_pert)
n_lev  = numel(levels);
L_ref  = find(abs(levels - noise_reference_factor) < 1e-12, 1);
ext_lbl = arrayfun(@(f) sprintf('%g mm', f), fwhms, 'UniformOutput', false);
ext_lbl(fwhms == 0) = {'focal'};
mlab   = struct('bem', 'BEM', 'fem', 'FEM', 'bslaw', 'Biot-Savart');
mstyle = struct('bem', '-', 'fem', '--', 'bslaw', ':');
mcol   = struct('bem', fwd_method_colors(3, :), 'fem', fwd_method_colors(4, :), ...
                'bslaw', fwd_method_colors(1, :));


% =========================================================================
% LOAD, AND BUILD THE SOURCE MODEL PER MODALITY
% =========================================================================

t  = (0:round(cord_duration * cord_fs) - 1) / cord_fs;
D  = struct();
for m = pert_modalities
    md = m{1};
    f1 = fullfile(mods_cfg.(md).forward_fields_base, 'leadfields_organised.mat');
    f2 = fullfile(mods_cfg.(md).forward_fields_base, 'pert_metrics.mat');
    if ~isfile(f1) || ~isfile(f2)
        warning('Stage 7: %s or %s missing — %s skipped.', f1, f2, upper(md));
        continue
    end
    tmp = load(f1, 'leadfields'); D.(md).LF = tmp.leadfields;
    tmp = load(f2, 'PM');         D.(md).PM = tmp.PM;
    [D.(md).Q, W] = cord_sources(D.(md).PM.dist_mm, t, fwhms, struct( ...
        'entry', cord_n11_entry_mm, 'end_mm', cord_n11_end_mm, ...
        'onset', cord_n11_onset, 'v_mm_s', cord_n11_velocity * 1000, ...
        'lag', cord_n11_lag_mm, 'nAm', cord_n11_nAm, 'ramp', cord_n11_ramp, ...
        'c7', evoked_src_mm, 'lat', noise_latency, 'sd', noise_env_sd, ...
        'freq', noise_freq, 'peak', noise_peak_nAm));
    D.(md).ref = [base_geom_name '_source_original'];
end
mods = fieldnames(D)';
if isempty(mods), error('Stage 7: no modality could be loaded.'); end

% Component windows (the same for every modality: they depend on time only)
win = struct('name', {'N11', 'N13'}, 't', {W.n11, W.n13});
n_win = numel(win);
for w = 1:n_win
    win(w).idx = t >= win(w).t(1) & t <= win(w).t(2);
end
fprintf('  N11 window %.1f-%.1f ms, N13 window %.1f-%.1f ms, %d extents\n\n', ...
    W.n11 * 1e3, W.n13 * 1e3, n_ext);

% Recording of one lead field set: [channels x time], whole array
rec = @(E, Q, PMx) pt_lf_matrix(E, cord_n11_ori, 0, PMx.n_trunc, PMx.src_range) * Q.n11 + ...
                   pt_lf_matrix(E, cord_n13_ori, 0, PMx.n_trunc, PMx.src_range) * Q.n13;


% =========================================================================
% SOURCE MODEL AND RECORDINGS (figures)
% =========================================================================

fprintf('[0] Source model and recordings\n');
md0 = mods{1}; Q0 = D.(md0).Q(e_head); dist0 = D.(md0).PM.dist_mm;
fig = figure('Color', 'w', 'Position', [40 40 1100 420]);
tl  = tiledlayout(1, 2, 'TileSpacing', 'compact');
comp = {'n11', sprintf('N11 volley (%s, nA*m)', cord_n11_ori); ...
        'n13', sprintf('N13, %s extent (%s, nA*m)', ext_lbl{e_head}, cord_n13_ori)};
for k = 1:2
    ax = nexttile(tl);
    imagesc(ax, t * 1e3, dist0, Q0.(comp{k, 1}));
    lim = max(abs(Q0.(comp{k, 1})(:))); if lim == 0, lim = 1; end
    caxis(ax, [-lim lim]); colormap(ax, diverging_map()); colorbar(ax);
    xline(ax, win(k).t * 1e3, ':k');
    yline(ax, evoked_src_mm, '--k', evoked_src_name, 'LabelHorizontalAlignment', 'left');
    xlabel(ax, 'Time (ms)'); ylabel(ax, 'Distance along cord (mm, rostral at top)');
    title(ax, comp{k, 2});
end
title(tl, 'Simulated cord activity after median nerve stimulation', 'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'source_model', fig_resolution);

% Recordings as field strength over time: the RMS over every channel of the
% trial-averaged recording (a global field power), which is how an evoked
% response is read in real data and can be set beside an empirical average.
% Rows = sensor technology, columns = forward model; in each panel the
% noise-free recording, three noisy realisations at 1x the floor, and the
% expected RMS of the noise alone.
techs = unique(regexprep({noise_systems.short}, '_(bs|fem)$', ''), 'stable');
for md = mods
    md = md{1}; X = D.(md); %#ok<FXSET>
    tk = techs(cellfun(@(s) any(strcmp({noise_systems.short}, s) & ...
        strcmp({noise_systems.modality}, md)), techs));
    meths = intersect({'bem', 'fem', 'bslaw'}, methods_with(X.LF, X.ref), 'stable');
    if isempty(tk) || isempty(meths), continue; end
    unit_s = mods_cfg.(md).unit(1:find(mods_cfg.(md).unit == '/', 1) - 1);
    fig = figure('Color', 'w', 'Position', [40 40 380 * numel(meths) + 120, 300 * numel(tk) + 140]);
    tl  = tiledlayout(numel(tk), numel(meths), 'TileSpacing', 'compact');
    for q = 1:numel(tk)
        base = noise_systems(strcmp({noise_systems.short}, tk{q}));
        bw   = min([base.bandwidth_hz, cord_fs / 2, cord_max_bandwidth_hz]);
        sd   = base.density * sqrt(bw) / sqrt(noise_n_trials);   % noise s.d. per sample at 1x
        for k = 1:numel(meths)
            ax = nexttile(tl); hold(ax, 'on');
            k_sys = find(strcmp({noise_systems.modality}, md) & ...
                strcmp(regexprep({noise_systems.short}, '_(bs|fem)$', ''), tk{q}) & ...
                strcmp({noise_systems.method}, meths{k}), 1);
            if isempty(k_sys)
                text(ax, 0.5, 0.5, sprintf('%s not simulated\nwith %s', base.label, mlab.(meths{k})), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5], ...
                    'Interpreter', 'none');
                axis(ax, 'off'); continue
            end
            sys = noise_systems(k_sys);
            Y0  = lowpass_fft(rec(X.LF.([meths{k} '_' X.ref]), X.Q(e_head), X.PM), bw, cord_fs);
            for w = 1:n_win
                patch(ax, 1e3 * win(w).t([1 2 2 1]), [0 0 1 1] * 1e9, [0.92 0.92 0.92], ...
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
            rng(noise_seed + q);
            for r = 1:3
                En = lowpass_fft(base.density * sqrt(cord_fs / 2) / sqrt(noise_n_trials) * ...
                    randn(size(Y0)), bw, cord_fs);
                plot(ax, t * 1e3, sqrt(mean((Y0 + En).^2, 1)), '-', 'Color', [sys.color 0.45], ...
                    'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
            plot(ax, t * 1e3, sqrt(mean(Y0.^2, 1)), '-', 'Color', 'k', 'LineWidth', 2, ...
                'DisplayName', 'noise-free');
            yline(ax, sd, '--', 'Color', sys.color, 'LineWidth', 1.4, 'HandleVisibility', 'off');
            ylim(ax, [0, max([sqrt(mean(Y0.^2, 1)), 2.5 * sd]) * 1.1]);
            xlim(ax, 1e3 * [t(1), t(end)]);
            grid(ax, 'on'); set(ax, 'TickDir', 'out', 'FontSize', 9, 'Layer', 'top');
            if q == 1, title(ax, mlab.(meths{k})); end
            if k == 1, ylabel(ax, sprintf('%s\nRMS over channels (%s)', base.label, unit_s)); end
            if q == numel(tk), xlabel(ax, 'Time after stimulus (ms)'); end
            text(ax, 1e3 * win(1).t(1), 0, ' N11', 'VerticalAlignment', 'bottom', 'FontSize', 8);
            text(ax, 1e3 * win(2).t(1), 0, ' N13', 'VerticalAlignment', 'bottom', 'FontSize', 8);
        end
    end
    title(tl, {sprintf(['%s recordings, N13 %s wide, %g trials — RMS over channels'], ...
        upper(md), ext_lbl{e_head}, noise_n_trials), ...
        'black = noise-free; coloured = 3 noisy realisations at 1\times floor; dashed = noise alone; grey = N11, N13 windows'}, ...
        'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['recordings_' md], fig_resolution);
end


% =========================================================================
% A0. DOES SPREADING CHANGE THE RECORDING AT ALL? (each model, vs its focal N13)
% =========================================================================
% Every later comparison is a ratio of two recordings of the same activity,
% so it can only depend on extent as far as spreading changes the recording
% itself. This measures that directly: each extent against the focal N13,
% within one model. Values of a few percent mean every extent effect below
% will be small, and that is the physics (a deep source, blurred along the
% cord by much less than its depth), not a fault.

fprintf('[A0] Spread vs focal\n');
e_foc = find(fwhms == 0, 1);
if isempty(e_foc), e_foc = 1; end
sp_rows = {};
for md = mods
    md = md{1}; X = D.(md); %#ok<FXSET>
    for mt = intersect({'bem', 'fem', 'bslaw'}, methods_with(X.LF, X.ref), 'stable')
        E  = X.LF.([mt{1} '_' X.ref]);
        Y0 = rec(E, X.Q(e_foc), X.PM);
        for e = setdiff(1:n_ext, e_foc)
            Ye = rec(E, X.Q(e), X.PM);
            for w = 1:n_win
                M = vec_metrics(Y0(:, win(w).idx), Ye(:, win(w).idx), mopts);
                sp_rows(end+1, :) = {upper(md), mlab.(mt{1}), mt{1}, fwhms(e), ...
                    win(w).name, M.re, M.rsq, M.rdm, M.gain}; %#ok<SAGROW>
            end
        end
    end
end
Tsp = cell2table(sp_rows, 'VariableNames', {'modality', 'model', 'method', ...
    'n13_fwhm_mm', 'window', 're', 'rsq', 'rdm', 'gain'});
Tsp.group = strcat(Tsp.modality, {' — '}, Tsp.model);
pt_write_table(Tsp, fullfile(out_dir, 'cord_spread_effect'), struct( ...
    'cols', {{ 'window', 'Window', 'Window', '%s'; ...
               'n13_fwhm_mm', 'N13 FWHM mm', 'N13 FWHM (mm)', '%g'; ...
               're', 'RE %', 'RE (\%)', '%.2f'; ...
               'rsq', 'r2', '$r^2$', '%.4f'; ...
               'rdm', 'RDM', 'RDM', '%.4f'; ...
               'gain', 'gain %', 'Gain (\%)', '%+.2f' }}, ...
    'title', 'STAGE 7 — HOW MUCH SPREADING THE N13 CHANGES THE RECORDING (each model vs its focal N13)', ...
    'notes', {{'Reference = the same model''s recording with a focal N13; same total moment at every extent.', ...
               'Small values here bound every extent effect in the divergence, detectability and perturbation outputs.'}}, ...
    'caption', ['Change in the simulated recording when the N13 is spread along ' ...
               'the cord, relative to a focal N13 of the same total moment, per ' ...
               'forward model.'], ...
    'label', 'tab:cord-spread', 'group_by', 'group'));

if ~isempty(sp_rows)
    ext_sp = setdiff(1:n_ext, e_foc);
    pm_ = {'re', 'RE vs focal (%)'; 'rdm', 'RDM vs focal'; 'gain', 'Gain vs focal (%)'};
    fig = figure('Color', 'w', 'Position', [40 40 420 * n_win, 300 * size(pm_, 1) + 120]);
    tl  = tiledlayout(size(pm_, 1), n_win, 'TileSpacing', 'compact');
    grp = unique(strcat(Tsp.modality, '|', Tsp.method), 'stable');
    for r = 1:size(pm_, 1)
        for w = 1:n_win
            ax = nexttile(tl); hold(ax, 'on');
            for q = 1:numel(grp)
                gm = strsplit(grp{q}, '|');
                sel = strcmp(Tsp.modality, gm{1}) & strcmp(Tsp.method, gm{2}) & ...
                      strcmp(Tsp.window, win(w).name);
                plot(ax, 1:sum(sel), Tsp.(pm_{r, 1})(sel), mstyle.(gm{2}), ...
                    'Color', mcol.(gm{2}), 'LineWidth', 1.6, ...
                    'DisplayName', sprintf('%s (%s)', gm{1}, mlab.(gm{2})));
            end
            ext_axis(ax, ext_lbl(ext_sp), r == size(pm_, 1));
            if r == 1, title(ax, sprintf('%s window', win(w).name)); end
            if w == 1, ylabel(ax, pm_{r, 2}); end
            if r == 1 && w == n_win, legend(ax, 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, ['How much spreading the N13 changes the recording — each model against ' ...
        'its focal N13 (noise-free, whole array)'], 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, 'spread_effect', fig_resolution);
end


% =========================================================================
% A. MODEL DIVERGENCE VS N13 EXTENT (noise-free, MRI-derived BEM reference)
% =========================================================================

fprintf('[A] Model divergence\n');
div_rows = {};
for md = mods
    md = md{1}; X = D.(md); %#ok<FXSET>
    if ~isfield(X.LF, ['bem_' X.ref]), continue; end
    for cm = intersect({'fem', 'bslaw'}, methods_with(X.LF, X.ref), 'stable')
        for e = 1:n_ext
            Ya = rec(X.LF.(['bem_' X.ref]), X.Q(e), X.PM);
            Yb = rec(X.LF.([cm{1} '_' X.ref]), X.Q(e), X.PM);
            for w = 1:n_win
                M = vec_metrics(Ya(:, win(w).idx), Yb(:, win(w).idx), mopts);
                div_rows(end+1, :) = {upper(md), sprintf('%s vs BEM', mlab.(cm{1})), ...
                    cm{1}, fwhms(e), win(w).name, M.re, M.rsq, M.rdm, M.gain}; %#ok<SAGROW>
            end
        end
    end
end
Tdiv = cell2table(div_rows, 'VariableNames', {'modality', 'comparison', 'method', ...
    'n13_fwhm_mm', 'window', 're', 'rsq', 'rdm', 'gain'});
Tdiv.group = strcat(Tdiv.modality, {' — '}, Tdiv.comparison);
pt_write_table(Tdiv, fullfile(out_dir, 'cord_divergence'), struct( ...
    'cols', {{ 'window', 'Window', 'Window', '%s'; ...
               'n13_fwhm_mm', 'N13 FWHM mm', 'N13 FWHM (mm)', '%g'; ...
               're', 'RE %', 'RE (\%)', '%.2f'; ...
               'rsq', 'r2', '$r^2$', '%.4f'; ...
               'rdm', 'RDM', 'RDM', '%.4f'; ...
               'gain', 'gain %', 'Gain (\%)', '%+.2f' }}, ...
    'title', 'STAGE 7 — MODEL DIVERGENCE VS N13 EXTENT (noise-free recordings, whole array)', ...
    'notes', {{'Reference = MRI-derived BEM recording of the same activity. Metrics on sensors x time within each window.', ...
               'The N11 is the same at every extent; its rows show what spreading the N13 does to the overlap.'}}, ...
    'caption', ['Difference between the recordings predicted by each forward model ' ...
               'and the MRI-derived BEM, for the N11 volley and the N13 at each ' ...
               'spatial extent.'], ...
    'label', 'tab:cord-divergence', 'group_by', 'group'));

if ~isempty(div_rows)
    pm_ = {'re', 'RE (%)'; 'rdm', 'RDM'; 'gain', 'Gain (%)'};
    fig = figure('Color', 'w', 'Position', [40 40 420 * n_win, 300 * size(pm_, 1) + 120]);
    tl  = tiledlayout(size(pm_, 1), n_win, 'TileSpacing', 'compact');
    cmps = unique(Tdiv.comparison, 'stable');
    for r = 1:size(pm_, 1)
        for w = 1:n_win
            ax = nexttile(tl); hold(ax, 'on');
            for c = 1:numel(cmps)
                for md = unique(Tdiv.modality, 'stable')'
                    sel = strcmp(Tdiv.comparison, cmps{c}) & strcmp(Tdiv.window, win(w).name) & ...
                          strcmp(Tdiv.modality, md{1});
                    if ~any(sel), continue; end
                    mt = Tdiv.method{find(sel, 1)};
                    plot(ax, 1:n_ext, Tdiv.(pm_{r, 1})(sel), ['o' mstyle.(mt)], ...
                        'Color', mcol.(mt), 'MarkerFaceColor', mcol.(mt), 'LineWidth', 1.6, ...
                        'DisplayName', sprintf('%s %s', md{1}, cmps{c}));
                end
            end
            ext_axis(ax, ext_lbl, r == size(pm_, 1));
            if r == 1, title(ax, sprintf('%s window', win(w).name)); end
            if w == 1, ylabel(ax, pm_{r, 2}); end
            if r == 1 && w == n_win, legend(ax, 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, ['Model divergence vs N13 extent — each model against the MRI-derived ' ...
        'BEM (noise-free, whole array)'], 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, 'divergence', fig_resolution);
end


% =========================================================================
% B. DETECTABILITY VS N13 EXTENT, PER FORWARD MODEL AND SYSTEM
% =========================================================================

fprintf('[B] Detectability\n');
all_short = {noise_systems.short};
det_rows = {}; crit_rows = {}; pair_rows = {};
RE_store = struct();
for k = 1:numel(noise_systems)
    sys = noise_systems(k); md = sys.modality;
    if ~isfield(D, md), continue; end
    X = D.(md); key = [sys.method '_' X.ref];
    if ~isfield(X.LF, key)
        fprintf('  %s: %s not loaded — skipped\n', sys.label, key);
        continue
    end
    fprintf('  %s\n', sys.label);
    bw    = min([sys.bandwidth_hz, cord_fs / 2, cord_max_bandwidth_hz]);
    sd0   = sys.density * sqrt(cord_fs / 2) / sqrt(noise_n_trials);   % per sample, before low-pass
    % Shared noise draws: same stream for every system with these sensors
    k_rng = find(strcmp(all_short, regexprep(sys.short, '_(bs|fem)$', '')), 1);
    rng(noise_seed + 104729 * k_rng);
    n_ch  = size(rec(X.LF.(key), X.Q(1), X.PM), 1);
    Z     = zeros(n_ch, numel(t), noise_n_real);
    for r = 1:noise_n_real
        Z(:, :, r) = lowpass_fft(sd0 * randn(n_ch, numel(t)), bw, cord_fs);
    end
    REk = nan(n_ext, n_lev, noise_n_real, n_win);
    for e = 1:n_ext
        Y0 = lowpass_fft(rec(X.LF.(key), X.Q(e), X.PM), bw, cord_fs);
        for w = 1:n_win
            iw  = win(w).idx;
            snr = sqrt(mean(Y0(:, iw).^2, 'all')) / (sys.density * sqrt(bw) / sqrt(noise_n_trials));
            rate = nan(1, n_lev);
            for L = 1:n_lev
                stat = nan(1, noise_n_real); null = stat;
                for r = 1:noise_n_real
                    En = levels(L) * Z(:, :, r);
                    Ym = Y0 + En;
                    M  = vec_metrics(Y0(:, iw), Ym(:, iw), mopts);
                    REk(e, L, r, w) = M.re;
                    stat(r) = mean(Ym(:, iw).^2, 'all');
                    null(r) = mean(En(:, iw).^2, 'all');   % noise alone, same window
                end
                if levels(L) == 0
                    rate(L) = 1;
                else
                    rate(L) = mean(stat > pctl(null, noise_detect_pct));
                end
                if levels(L) > 0, snr_L = snr / levels(L); else, snr_L = Inf; end
                det_rows(end+1, :) = {sys.label, sys.short, upper(md), mlab.(sys.method), ...
                    fwhms(e), win(w).name, levels(L), snr_L, median(REk(e, L, :, w), 'omitnan'), ...
                    pctl(squeeze(REk(e, L, :, w)), 25), pctl(squeeze(REk(e, L, :, w)), 75), ...
                    rate(L)}; %#ok<SAGROW>
            end
            [crit, ctxt] = pt_critical_level(levels, rate, noise_detect_rate);
            crit_rows(end+1, :) = {sys.label, sys.short, upper(md), mlab.(sys.method), ...
                fwhms(e), win(w).name, snr, crit, ctxt}; %#ok<SAGROW>
        end
    end
    RE_store.(sys.short) = REk;
end
Tdet = cell2table(det_rows, 'VariableNames', {'system', 'short', 'modality', 'model', ...
    'n13_fwhm_mm', 'window', 'noise_factor', 'snr', 're_median', 're_q1', 're_q3', 'detect_rate'});
writetable(Tdet, fullfile(out_dir, 'cord_detectability.csv'));

Tcrit = cell2table(crit_rows, 'VariableNames', {'system', 'short', 'modality', 'model', ...
    'n13_fwhm_mm', 'window', 'snr_at_floor', 'critical_factor', 'critical_text'});
Tcrit.group = strcat(Tcrit.window, {' — '}, Tcrit.system);
pt_write_table(Tcrit, fullfile(out_dir, 'cord_critical_levels'), struct( ...
    'cols', {{ 'n13_fwhm_mm', 'N13 FWHM mm', 'N13 FWHM (mm)', '%g'; ...
               'snr_at_floor', 'SNR at 1x', 'SNR at 1$\times$', '%.3g'; ...
               'critical_text', 'critical x floor', 'Critical $\times$ floor', '%s' }}, ...
    'title', 'STAGE 7 — CRITICAL NOISE LEVEL PER COMPONENT VS N13 EXTENT (whole array)', ...
    'notes', {{sprintf(['Highest noise level (x each system''s floor) at which >= %.0f%% of noise ' ...
               'realisations still exceed the %gth percentile of noise alone in the window.'], ...
               noise_detect_rate * 100, noise_detect_pct), ...
               'SNR = rms of the noise-free recording in the window / noise s.d. per sample at 1x, after the system bandwidth.'}}, ...
    'caption', ['Noise level up to which the N11 and N13 remain detectable, per ' ...
               'sensor system and forward model, as the N13 spreads along the cord.'], ...
    'label', 'tab:cord-critical', 'group_by', 'group'));

% FEM vs BEM recovery error, paired by noise draw, at the reference level
for k = 1:numel(noise_systems)
    sys = noise_systems(k);
    if ~strcmp(sys.method, 'fem') || ~isfield(RE_store, sys.short), continue; end
    bem_short = regexprep(sys.short, '_fem$', '');
    if ~isfield(RE_store, bem_short), continue; end
    for e = 1:n_ext
        for w = 1:n_win
            a = squeeze(RE_store.(sys.short)(e, L_ref, :, w));
            b = squeeze(RE_store.(bem_short)(e, L_ref, :, w));
            R = pt_compare_two(a, b, true, stats_n_boot, stats_ci_level);
            pair_rows(end+1, :) = pt_compare_row(R, 'fem_vs_bem', ...
                sprintf('%s, %s', win(w).name, ext_lbl{e}), ...
                sprintf('%s: FEM vs BEM', regexprep(sys.label, ' \((FEM|BEM)\)$', ''))); %#ok<SAGROW>
        end
    end
end
if ~isempty(pair_rows)
    Tpair = cell2table(pair_rows, 'VariableNames', pt_compare_row());
    pt_write_table(Tpair, fullfile(out_dir, 'cord_bem_vs_fem_noise'), struct( ...
        'cols', {{ 'comparison', 'System', 'System', '%s'; ...
                   'subset', 'Window, extent', 'Window, extent', '%s'; ...
                   'median_a', 'FEM RE', 'FEM RE', '%.3g'; ...
                   'median_b', 'BEM RE', 'BEM RE', '%.3g'; ...
                   'estimate', 'FEM-BEM [95% CI]', 'FEM$-$BEM [95\% CI]', 'eci'; ...
                   'separated', 'CI excl. 0', 'CI excl.\ 0', 'yn' }}, ...
        'title', sprintf('STAGE 7 — RECOVERY ERROR THROUGH NOISE, FEM VS BEM (%gx floor)', noise_reference_factor), ...
        'notes', {{'RE of the measured recording against the noise-free recording of the same model, in each window.', ...
                   'Paired by noise draw (identical draws for both models). Negative = FEM recording recovered more accurately.'}}, ...
        'caption', sprintf(['Recovery error of the N11 and N13 through sensor noise at ' ...
                   '%g$\\times$ the noise floor, FEM against BEM, per N13 extent.'], ...
                   noise_reference_factor), ...
        'label', 'tab:cord-fem-bem-noise'));
end

% SNR against N13 extent: the exact detectability measure. It is computed
% from the noise-free recording and the noise level, with no noise draws,
% so it resolves the small changes that spreading produces; the critical
% level, built from noise_n_real draws per level, moves only in steps and
% stays the same when the SNR changes by a few percent.
snr_rows = {};
for k = 1:numel(noise_systems)
    sys = noise_systems(k);
    Tk  = Tcrit(strcmp(Tcrit.short, sys.short), :);
    if isempty(Tk), continue; end
    Tb  = Tcrit(strcmp(Tcrit.short, regexprep(sys.short, '_(bs|fem)$', '')), :);
    for w = 1:n_win
        Tw_ = Tk(strcmp(Tk.window, win(w).name), :);
        f0  = Tw_.snr_at_floor(Tw_.n13_fwhm_mm == fwhms(e_foc));
        if isempty(f0), f0 = NaN; end
        for i = 1:height(Tw_)
            sb = Tb.snr_at_floor(strcmp(Tb.window, win(w).name) & ...
                                 Tb.n13_fwhm_mm == Tw_.n13_fwhm_mm(i));
            if isempty(sb), sb = NaN; end
            snr_rows(end+1, :) = {sys.label, sys.short, Tw_.model{i}, win(w).name, ...
                Tw_.n13_fwhm_mm(i), Tw_.snr_at_floor(i), 100 * Tw_.snr_at_floor(i) / f0, ...
                Tw_.snr_at_floor(i) / sb(1)}; %#ok<SAGROW>
        end
    end
end
Tsnr = cell2table(snr_rows, 'VariableNames', {'system', 'short', 'model', 'window', ...
    'n13_fwhm_mm', 'snr_at_floor', 'snr_pct_of_focal', 'snr_ratio_to_bem'});
Tsnr.group = strcat(Tsnr.window, {' — '}, Tsnr.system);
pt_write_table(Tsnr, fullfile(out_dir, 'cord_snr'), struct( ...
    'cols', {{ 'n13_fwhm_mm', 'N13 FWHM mm', 'N13 FWHM (mm)', '%g'; ...
               'snr_at_floor', 'SNR at 1x', 'SNR at 1$\times$', '%.4g'; ...
               'snr_pct_of_focal', '% of focal', '\% of focal', '%.1f'; ...
               'snr_ratio_to_bem', 'vs BEM', 'vs BEM', '%.3f' }}, ...
    'title', 'STAGE 7 — SIGNAL-TO-NOISE OF EACH COMPONENT VS N13 EXTENT (whole array)', ...
    'notes', {{'SNR = rms of the noise-free recording in the window / noise s.d. per sample at 1x the floor, after the system bandwidth. Exact: no noise draws.', ...
               '% of focal: against the same system with a focal N13. vs BEM: against the BEM system with the same sensors (1 = identical).'}}, ...
    'caption', ['Signal-to-noise of the N11 and N13 per sensor system and forward ' ...
               'model, as the N13 spreads along the cord.'], ...
    'label', 'tab:cord-snr', 'group_by', 'group'));

% Against extent: relative to the focal N13, and each model relative to the
% BEM, where a few percent is visible (absolute SNR is in snr_by_system)
pm_ = {'snr_pct_of_focal', 'SNR, % of focal N13'; 'snr_ratio_to_bem', 'SNR relative to BEM'};
fig = figure('Color', 'w', 'Position', [40 40 460 * n_win + 160, 340 * size(pm_, 1) + 120]);
tl  = tiledlayout(size(pm_, 1), n_win, 'TileSpacing', 'compact');
for r = 1:size(pm_, 1)
    for w = 1:n_win
        ax = nexttile(tl); hold(ax, 'on');
        for k = 1:numel(noise_systems)
            sys = noise_systems(k);
            sel = strcmp(Tsnr.short, sys.short) & strcmp(Tsnr.window, win(w).name);
            if ~any(sel), continue; end
            plot(ax, 1:sum(sel), Tsnr.(pm_{r, 1})(sel), mstyle.(sys.method), ...
                'Color', sys.color, 'LineWidth', 1.6, 'DisplayName', sys.label);
        end
        ext_axis(ax, ext_lbl, r == size(pm_, 1));
        if r == 1, title(ax, sprintf('%s window', win(w).name)); end
        if w == 1, ylabel(ax, pm_{r, 2}); end
        if r == 1 && w == n_win, legend(ax, 'Location', 'bestoutside', 'Box', 'off'); end
    end
end
title(tl, {'Signal-to-noise vs N13 extent — exact, no noise draws (whole array)', ...
    'top: against the same system with a focal N13; bottom: against the BEM system with the same sensors (1 = identical)'}, ...
    'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'detectability_snr', fig_resolution);

% SNR at 1x per sensor technology and model, at the headline extent: the
% direct comparison of systems
sel_e = Tsnr.n13_fwhm_mm == fwhms(e_head);
techs_all = unique(regexprep(Tsnr.short(sel_e), '_(bs|fem)$', ''), 'stable');
mod_order = {'bem', 'fem', 'bslaw'};
fig = figure('Color', 'w', 'Position', [40 40 460 * n_win + 160, 420]);
tl  = tiledlayout(1, n_win, 'TileSpacing', 'compact');
for w = 1:n_win
    ax = nexttile(tl); hold(ax, 'on');
    Mv = nan(numel(techs_all), numel(mod_order));
    for q = 1:numel(techs_all)
        for j = 1:numel(mod_order)
            sh = techs_all{q};
            if ~strcmp(mod_order{j}, 'bem'), sh = [sh '_' strrep(mod_order{j}, 'bslaw', 'bs')]; end
            v = Tsnr.snr_at_floor(sel_e & strcmp(Tsnr.short, sh) & strcmp(Tsnr.window, win(w).name));
            if ~isempty(v), Mv(q, j) = v(1); end
        end
    end
    hb = bar(ax, Mv, 'grouped');
    for j = 1:numel(hb)
        hb(j).FaceColor = mcol.(mod_order{j});
        hb(j).HandleVisibility = 'off';
        ok = ~isnan(Mv(:, j));
        if any(ok)
            % A shape per model on each bar, so the bars read without colour
            scatter(ax, hb(j).XEndPoints(ok), hb(j).YEndPoints(ok), 45, 'k', pt_marker(j), ...
                'filled', 'MarkerFaceColor', 'w', 'DisplayName', mlab.(mod_order{j}));
        end
    end
    tl_lbl = cellfun(@(sh) noise_systems(strcmp({noise_systems.short}, sh)).label, techs_all, ...
        'UniformOutput', false);
    set(ax, 'XTick', 1:numel(techs_all), 'XTickLabel', regexprep(tl_lbl, ' \(BEM\)$', ''), ...
        'TickDir', 'out', 'FontSize', 9);
    grid(ax, 'on');
    title(ax, sprintf('%s window', win(w).name));
    if w == 1, ylabel(ax, 'SNR at 1\times floor'); end
    if w == n_win, legend(ax, 'Location', 'best', 'Box', 'off'); end
end
title(tl, sprintf(['Signal-to-noise per sensor system and forward model — N13 %s wide, ' ...
    '%g trials'], ext_lbl{e_head}, noise_n_trials), 'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'snr_by_system', fig_resolution);

% Detection rate against noise level at the headline extent: how far each
% system is from detecting each component
lev_lbl7 = arrayfun(@(x) sprintf('%gx', x), levels, 'UniformOutput', false);
lev_lbl7{1} = '0';
tick_i = unique([1, find(abs(log2(max(levels, eps)) - round(log2(max(levels, eps)))) < 1e-9)]);
fig = figure('Color', 'w', 'Position', [40 40 520 * n_win + 200, 440]);
tl  = tiledlayout(1, n_win, 'TileSpacing', 'compact');
for w = 1:n_win
    ax = nexttile(tl); hold(ax, 'on');
    for k = 1:numel(noise_systems)
        sys = noise_systems(k);
        sel = strcmp(Tdet.short, sys.short) & strcmp(Tdet.window, win(w).name) & ...
              Tdet.n13_fwhm_mm == fwhms(e_head);
        if ~any(sel), continue; end
        plot(ax, 1:sum(sel), Tdet.detect_rate(sel), mstyle.(sys.method), 'Color', sys.color, ...
            'LineWidth', 1.6, 'DisplayName', sys.label);
    end
    xline(ax, find(abs(levels - 1) < 1e-12, 1), ':k', '1\times floor', 'HandleVisibility', 'off');
    yline(ax, noise_detect_rate, '--k', sprintf('%g%% detected', 100 * noise_detect_rate), ...
        'HandleVisibility', 'off');
    set(ax, 'XTick', tick_i, 'XTickLabel', lev_lbl7(tick_i), 'XLim', [0.5, n_lev + 0.5], ...
        'YLim', [-0.02 1.02], 'TickDir', 'out', 'FontSize', 9);
    grid(ax, 'on');
    title(ax, sprintf('%s window', win(w).name));
    xlabel(ax, 'Noise (\times system floor)');
    if w == 1, ylabel(ax, 'Fraction of realisations detected'); end
    if w == n_win, legend(ax, 'Location', 'bestoutside', 'Box', 'off'); end
end
title(tl, {sprintf('Detection against noise level — N13 %s wide, %g trials', ext_lbl{e_head}, ...
    noise_n_trials), 'solid = BEM, dashed = FEM, dotted = Biot-Savart'}, 'FontWeight', 'bold');
pt_save_figure(fig, out_dir, 'detection_vs_noise', fig_resolution);

% Figure per window: critical level vs N13 extent, systems x models
lv = levels(levels > 0);
for w = 1:n_win
    Tw = Tcrit(strcmp(Tcrit.window, win(w).name), :);
    if isempty(Tw), continue; end
    fig = figure('Color', 'w', 'Position', [40 40 760 480]);
    ax  = axes(fig); hold(ax, 'on');
    for k = 1:numel(noise_systems)
        sys = noise_systems(k);
        sel = strcmp(Tw.short, sys.short);
        if ~any(sel), continue; end
        y = level_position(Tw.critical_factor(sel), lv);
        plot(ax, 1:sum(sel), y, ['o' mstyle.(sys.method)], 'Color', sys.color, ...
            'MarkerFaceColor', sys.color, 'LineWidth', 1.6, 'DisplayName', sys.label);
    end
    % The simulated noise levels at evenly spaced positions, labelled with
    % their value, as on every stage-5 noise axis; 'none' and '> max' sit
    % half a step outside
    set(ax, 'YTick', [0.5, 1:numel(lv), numel(lv) + 0.5], 'YLim', [0.3, numel(lv) + 0.7], ...
        'YTickLabel', [{'none'}, arrayfun(@(x) sprintf('%gx', x), lv, 'UniformOutput', false), ...
                       {sprintf('> %gx', lv(end))}]);
    yline(ax, level_position(1, lv), ':k', 'realistic floor', 'HandleVisibility', 'off');
    ext_axis(ax, ext_lbl, true); ylabel(ax, 'Critical noise level (\times floor)');
    legend(ax, 'Location', 'bestoutside', 'Box', 'off');
    title(ax, {sprintf('%s: noise level up to which it is still detected, vs N13 extent', win(w).name), ...
        'solid = BEM, dashed = FEM, dotted = Biot-Savart; above 1 = visible in a real recording'});
    pt_save_figure(fig, out_dir, ['detectability_' win(w).name], fig_resolution);
end


% =========================================================================
% C. PERTURBATION SENSITIVITY VS N13 EXTENT (noise-free)
% =========================================================================

fprintf('[C] Perturbations\n');
pert_rows = {}; pb_rows = {};
sopts = struct('bundle_names', {bundle_names}, 'n_boot', stats_n_boot, 'level', stats_ci_level);
for md = mods
    md = md{1}; X = D.(md); PMx = X.PM; %#ok<FXSET>
    tys = pert_types(isfield(PMx.pert, pert_types));
    Sx  = struct();
    sets = {};                                   % {type, method, P}
    for ty = tys
        P = PMx.pert.(ty{1});
        meths = P.methods;
        if ismember(ty{1}, {'source', 'sensor'}) && isfield(X.LF, ['bslaw_' P.ref_key])
            meths = union(meths, {'bslaw'}, 'stable');
        end
        for m = meths, sets(end+1, :) = {ty{1}, m{1}, P}; end %#ok<AGROW>
    end
    if isfield(PMx, 'remesh'), sets(end+1, :) = {'remesh', 'fem', PMx.remesh}; end

    for q = 1:size(sets, 1)
        [ty, mt, P] = sets{q, :};
        rk = [mt '_' P.ref_key];
        if ~isfield(X.LF, rk), continue; end
        n = numel(P.keys);
        V = struct('re', nan(n, n_ext, n_win), 'rdm', nan(n, n_ext, n_win), ...
                   'absgain', nan(n, n_ext, n_win), 'rsq', nan(n, n_ext, n_win));
        for e = 1:n_ext
            Ya = rec(X.LF.(rk), X.Q(e), PMx);
            for i = 1:n
                kk = [mt '_' P.keys{i}];
                if ~isfield(X.LF, kk), continue; end
                Yb = rec(X.LF.(kk), X.Q(e), PMx);
                for w = 1:n_win
                    M = vec_metrics(Ya(:, win(w).idx), Yb(:, win(w).idx), mopts);
                    V.re(i, e, w) = M.re; V.rdm(i, e, w) = M.rdm;
                    V.absgain(i, e, w) = abs(M.gain); V.rsq(i, e, w) = M.rsq;
                end
            end
        end
        Sx.(ty).(mt) = V; Sx.(ty).bundle_idx = P.bundle_idx;
        for e = 1:n_ext
            for w = 1:n_win
                if strcmp(ty, 'remesh'), subs = 0; else, subs = [1:numel(bundle_names), 0]; end
                for b = subs
                    if b == 0, sel = true(n, 1); sub = 'pooled';
                    else, sel = P.bundle_idx == b; sub = bundle_names{b}; end
                    Dd = pt_describe(V.re(sel, e, w), stats_n_boot, stats_ci_level, stats_cover_pct);
                    pert_rows(end+1, :) = {upper(md), ty, mlab.(mt), sub, fwhms(e), ...
                        win(w).name, Dd.n, Dd.median, Dd.ci_lo, Dd.ci_hi, ...
                        median(V.rdm(sel, e, w), 'omitnan'), ...
                        median(V.absgain(sel, e, w), 'omitnan')}; %#ok<SAGROW>
                end
            end
        end
    end

    % FEM vs BEM, paired by realisation (same shift / same conductivity draw)
    for ty = tys
        if ~isfield(Sx, ty{1}) || ~all(isfield(Sx.(ty{1}), {'bem', 'fem'})), continue; end
        bi = Sx.(ty{1}).bundle_idx;
        for e = 1:n_ext
            for w = 1:n_win
                T = pt_stats_between(Sx.(ty{1}).fem.re(:, e, w), Sx.(ty{1}).bem.re(:, e, w), ...
                    bi, bi, true, sopts, 'FEM vs BEM');
                T.modality     = repmat({upper(md)}, height(T), 1);
                T.perturbation = repmat({pert_type_display.(ty{1})}, height(T), 1);
                T.n13_fwhm_mm  = repmat(fwhms(e), height(T), 1);
                T.window       = repmat({win(w).name}, height(T), 1);
                pb_rows{end+1} = T; %#ok<SAGROW>
            end
        end
    end
    D.(md).S = Sx;
end

Tp = cell2table(pert_rows, 'VariableNames', {'modality', 'perturbation', 'model', 'bundle', ...
    'n13_fwhm_mm', 'window', 'n', 're', 're_ci_lo', 're_ci_hi', 'rdm_median', 'absgain_median'});
writetable(Tp, fullfile(out_dir, 'cord_perturbation.csv'));

if ~isempty(pb_rows)
    Tpb = vertcat(pb_rows{:});
    Tpb = movevars(Tpb, {'modality', 'perturbation', 'window', 'n13_fwhm_mm'}, 'Before', 1);
    Tpb.group = strcat(Tpb.perturbation, {' — '}, Tpb.window);
    pt_write_table(Tpb, fullfile(out_dir, 'cord_perturbation_bem_vs_fem'), struct( ...
        'cols', {{ 'n13_fwhm_mm', 'N13 FWHM mm', 'N13 FWHM (mm)', '%g'; ...
                   'subset', 'Bundle', 'Bundle', '%s'; ...
                   'median_a', 'FEM RE', 'FEM RE', '%.3g'; ...
                   'median_b', 'BEM RE', 'BEM RE', '%.3g'; ...
                   'estimate', 'FEM-BEM [95% CI]', 'FEM$-$BEM [95\% CI]', 'eci'; ...
                   'separated', 'CI excl. 0', 'CI excl.\ 0', 'yn' }}, ...
        'title', 'STAGE 7 — PERTURBATION SENSITIVITY OF THE RECORDING, FEM VS BEM (noise-free)', ...
        'notes', {{'RE of the perturbed recording against the unperturbed recording of the same model, in each window.', ...
                   'Paired by realisation. Positive = FEM more sensitive. FEM source shifts include re-meshing (see cord_perturbation, model FEM, perturbation remesh).'}}, ...
        'caption', ['Sensitivity of the simulated N11 and N13 recordings to source, ' ...
                   'sensor and conductivity perturbations, FEM against BEM, per N13 extent.'], ...
        'label', 'tab:cord-pert-fem-bem', 'group_by', 'group', ...
        'tex_rows', strcmp(Tpb.subset, 'pooled')));
end

% Figure per modality: rows = windows, cols = types; median RE (pooled) vs extent
for md = mods
    md = md{1}; %#ok<FXSET>
    if ~isfield(D.(md), 'S'), continue; end
    Sx = D.(md).S; tys = pert_types(isfield(Sx, pert_types));
    if isempty(tys), continue; end
    fig = figure('Color', 'w', 'Position', [40 40 400 * numel(tys) + 100, 340 * n_win + 120]);
    tl  = tiledlayout(n_win, numel(tys), 'TileSpacing', 'compact');
    for w = 1:n_win
        for q = 1:numel(tys)
            ty = tys{q}; ax = nexttile(tl); hold(ax, 'on');
            if strcmp(ty, 'source') && isfield(Sx, 'remesh') && isfield(Sx.remesh, 'fem')
                pt_plot_band(ax, 1:n_ext, Sx.remesh.fem.re(:, :, w), [0.6 0.6 0.6], ':', ...
                    'FEM re-meshing');
            end
            for mt = intersect({'bem', 'fem', 'bslaw'}, fieldnames(Sx.(ty))', 'stable')
                pt_plot_band(ax, 1:n_ext, Sx.(ty).(mt{1}).re(:, :, w), mcol.(mt{1}), ...
                    ['-' 'o'], mlab.(mt{1}));
            end
            ext_axis(ax, ext_lbl, w == n_win);
            if w == 1, title(ax, pert_type_display.(ty)); end
            if q == 1, ylabel(ax, sprintf('%s window: RE (%%)', win(w).name)); end
            if w == 1 && q == numel(tys), legend(ax, 'Location', 'best', 'Box', 'off'); end
        end
    end
    title(tl, sprintf(['%s — perturbation effect on the recording vs N13 extent (noise-free, ' ...
        'whole array)\nmedian and IQR over all realisations'], upper(md)), 'FontWeight', 'bold');
    pt_save_figure(fig, out_dir, ['perturbation_' md], fig_resolution);
end

fprintf('\nStage 7 complete: %s\n', out_dir);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [Q, W] = cord_sources(dist, t, fwhms, p)
% Moments [n_src x n_t] in nA*m per orientation set: Q(e).n11 (the volley,
% identical at every extent) and Q(e).n13 (the N13 at extent fwhms(e)).
% W holds the component windows (s).
    n_src = numel(dist); T = numel(t);
    sp    = median(diff(dist));
    if isempty(p.end_mm), x_end = min(dist); else, x_end = p.end_mm; end
    travel = max(p.entry - x_end, 0) / p.v_mm_s;          % s, leading dipole
    lag_t  = p.lag / p.v_mm_s;                            % s, trailing behind
    n11 = zeros(n_src, T);
    for k = 1:T
        for side = 1:2
            if side == 1, t0 = p.onset; sgn = 1;
            else,         t0 = p.onset + lag_t; sgn = -1;
                          if p.lag <= 0, continue; end
            end
            tau = t(k) - t0;
            if tau < 0 || tau > travel, continue; end
            a = taper(tau, travel, p.ramp);
            x = p.entry - p.v_mm_s * tau;                 % ascending: towards smaller distance
            wgt = max(0, 1 - abs(dist(:) - x) / sp);     % linear between neighbouring sources
            n11(:, k) = n11(:, k) + sgn * p.nAm * a * wgt;
        end
    end
    w13 = p.peak * cos(2*pi*p.freq*(t - p.lat)) .* exp(-((t - p.lat).^2) / (2 * p.sd^2));
    Q = struct('n11', {}, 'n13', {});
    for e = 1:numel(fwhms)
        if fwhms(e) <= 0
            [~, i0] = min(abs(dist - p.c7));
            g = zeros(n_src, 1); g(i0) = 1;
        else
            s = fwhms(e) / (2 * sqrt(2 * log(2)));
            g = exp(-((dist(:) - p.c7).^2) / (2 * s^2));
            g = g / sum(g);                               % same total moment at every extent
        end
        Q(e).n11 = n11;
        Q(e).n13 = g * w13;
    end
    W.n11 = [p.onset, p.onset + travel + max(lag_t, 0)];
    W.n13 = [p.lat - 2 * p.sd, p.lat + 2 * p.sd];
end

function a = taper(tau, dur, ramp)
% Raised-cosine onset and offset over ramp seconds
    a = 1;
    if ramp <= 0, return; end
    if tau < ramp, a = 0.5 * (1 - cos(pi * tau / ramp)); end
    if dur - tau < ramp, a = min(a, 0.5 * (1 - cos(pi * max(dur - tau, 0) / ramp))); end
end

function p = level_position(c, lv)
% Position of a critical noise level on an axis that places the simulated
% levels lv at 1..n: interpolated between neighbouring levels, 0.5 when
% never detected, n + 0.5 when still detected at the largest level.
    p = nan(size(c));
    for i = 1:numel(c)
        if isinf(c(i)),       p(i) = numel(lv) + 0.5;
        elseif c(i) <= 0,     p(i) = 0.5;
        elseif c(i) <= lv(1), p(i) = 1;
        else, p(i) = interp1(log(lv), 1:numel(lv), log(min(c(i), lv(end))));
        end
    end
end

function Y = lowpass_fft(Y, fc, fs)
% Zero-phase low-pass along time: a sensor bandwidth applied to signal and
% noise alike. No toolbox needed.
    if ~isfinite(fc) || fc >= fs / 2, return; end
    T  = size(Y, 2);
    Tp = 2 * T;                                   % zero-pad: no wrap-around
    f  = (0:Tp-1) * fs / Tp; f = min(f, fs - f);
    F  = fft(Y, Tp, 2);
    F(:, f > fc) = 0;
    Y  = real(ifft(F, [], 2));
    Y  = Y(:, 1:T);
end

function M = vec_metrics(A, B, mopts)
% The four msg_fwd metrics on two recordings, flattened
    R = pt_metrics_block(A(:), B(:), mopts);
    M = struct('re', R.re, 'rsq', R.rsq, 'rdm', R.rdm, 'gain', (exp(R.lnmag) - 1) * 100);
end

function m = methods_with(LF, ref)
% Methods whose unperturbed lead field is loaded
    m = {};
    for c = {'bem', 'fem', 'bslaw', 'sphere'}
        if isfield(LF, [c{1} '_' ref]), m{end+1} = c{1}; end %#ok<AGROW>
    end
end

function ext_axis(ax, lbl, show_x)
    set(ax, 'XTick', 1:numel(lbl), 'XTickLabel', lbl, 'XLim', [0.6, numel(lbl) + 0.4], ...
        'TickDir', 'out', 'FontSize', 9);
    grid(ax, 'on');
    if show_x, xlabel(ax, 'N13 extent (FWHM along the cord)'); end
end

function cm = diverging_map()
    n = 128; s = linspace(0, 1, n)';
    cm = [[s, s, ones(n, 1)]; [ones(n, 1), flipud(s), flipud(s)]];
end

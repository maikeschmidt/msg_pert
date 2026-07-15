% sim_perturbation_noise - Noise-degradation curves for every perturbation
%
% Combines the two analyses in this repository. The perturbation pipeline (pt_*)
% asks how the forward field changes when the model is wrong; the simulation
% module asks how sensor noise degrades a measurement. This script asks both at
% once: for each perturbed forward field, how does adding realistic sensor noise
% degrade what you can recover?
%
% For every perturbation configuration (source shift, sensor shift, tissue
% conductivity), the perturbed forward field is treated as the true field, the
% evoked source is projected through it, trial-averaged sensor noise is added
% across the usual sweep, and r^2 is computed TWO ways:
%
%   self     — vs the noise-free version of the SAME perturbed field.
%              Isolates robustness to noise for that perturbed configuration.
%   perfect  — vs the unperturbed ("perfect") field.
%              Combines model error (the perturbation) and sensor noise. r^2 is
%              below 1 even at zero noise, because the perturbation alone
%              already differs from perfect.
%
% Results are averaged across the 8 shifts within each bundle, giving one small
% / medium / large curve per perturbation type — matching how the rest of the
% perturbation pipeline groups its results.
%
% USAGE:
%   sim_perturbation_noise           % then sim_plot_perturbation_noise
%
% SCOPE:
%   Uses leadfields_organised.mat from pt_load_leadfields, so it runs on ONE
%   modality at a time (whatever that run loaded — MSG or ESG). The applicable
%   noise systems are selected automatically: MSG runs use the magnetic systems
%   (SQUID, OP-MSG), ESG runs use the electric system (ESG). Perturbations are
%   BEM only (conductivity perturbation exists only for BEM), so sim_pert_method
%   is fixed to 'bem'.
%
% OUTPUT (to <sim_out_dir>/sim_pert_noise.mat):
%   For each perturbation type present, cord-averaged and per-source r^2 grouped
%   by bundle, for both references. See the save block for the variable list.
%
% DEPENDENCIES:
%   config_pert  (perturbation keys, bundles), config_sim (noise systems,
%   signal, trials, bandwidth), pt_add_functions, leadfields_organised.mat,
%   sim_noise_rsq_refs()
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
config_sim;
pt_add_functions;

% SETTINGS
sim_pert_method = 'bem';   % perturbations are BEM (conductivity is BEM-only)
sim_pert_n_real = 10;      % noise realisations per config (fewer than the
                           % perfect-field sim: bundle averaging over 8 shifts
                           % already smooths the result)

fprintf('sim_perturbation_noise\n');

% -------------------------------------------------------------------------
% Load organised perturbation leadfields
% -------------------------------------------------------------------------
lf_file = fullfile(forward_fields_base, 'leadfields_organised.mat');
if ~isfile(lf_file)
    error('Not found: %s\nRun pt_load_leadfields first.', lf_file);
end
S = load(lf_file, 'leadfields');
leadfields = S.leadfields;


% -------------------------------------------------------------------------
% Define the three perturbation families (keys, references, bundles)
% -------------------------------------------------------------------------
% Each entry lists: the unperturbed reference key, the per-shift perturbed
% keys, and the bundle index / display / colour for each shift.

perts = struct('name', {}, 'ref', {}, 'shift_keys', {}, ...
               'bundle_idx', {}, 'bundle_display', {}, 'bundle_colors', {});

% Source shifts
perts(1).name           = 'source';
perts(1).ref            = [sim_pert_method '_' sensitivity_ref_key];
perts(1).shift_keys     = cellfun(@(k) [sim_pert_method '_' k], ...
                                  sensitivity_keys, 'UniformOutput', false);
perts(1).bundle_idx     = source_bundle_idx;
perts(1).bundle_display = source_bundle_display;
perts(1).bundle_colors  = source_bundle_colors;

% Sensor shifts
perts(2).name           = 'sensor';
perts(2).ref            = [sim_pert_method '_' sensor_sensitivity_ref_key];
perts(2).shift_keys     = cellfun(@(k) [sim_pert_method '_' k], ...
                                  sensor_sensitivity_keys, 'UniformOutput', false);
perts(2).bundle_idx     = sensor_sensitivity_bundle_idx;
perts(2).bundle_display = sensor_bundle_display;
perts(2).bundle_colors  = sensor_bundle_colors;

% Conductivity perturbations (BEM-cond keys have their own naming)
cond_shift_keys = cell(1, n_cond_total);
for b = 1:n_cond_bundles
    for s = 1:n_cond_shifts
        idx = (b-1)*n_cond_shifts + s;
        cond_shift_keys{idx} = sprintf('bem_cond_%s_bundle%d_shift%d', ...
            cond_sensitivity_ref_key, b, s);
    end
end
perts(3).name           = 'cond';
perts(3).ref            = [sim_pert_method '_' cond_sensitivity_ref_key];
perts(3).shift_keys     = cond_shift_keys;
perts(3).bundle_idx     = cond_bundle_idx;
perts(3).bundle_display = cond_bundle_display;
perts(3).bundle_colors  = cond_bundle_colors;


% -------------------------------------------------------------------------
% Pick the noise systems that match the loaded modality
% -------------------------------------------------------------------------
% Modality is read from any available reference field.
ref_probe = '';
for p = 1:numel(perts)
    if isfield(leadfields, perts(p).ref)
        ref_probe = perts(p).ref; break
    end
end
if isempty(ref_probe)
    error(['No perturbation reference found in leadfields_organised.mat for ' ...
           'method "%s". Was this a %s run?'], sim_pert_method, sim_pert_method);
end
is_meg_run = leadfields.(ref_probe).is_meg;

sys_mask   = ([sim_models([sim_systems.model]).is_meg] == is_meg_run);
sys_use    = find(sys_mask);
if isempty(sys_use)
    error('No noise systems match the loaded modality (is_meg = %d).', is_meg_run);
end
n_sys      = numel(sys_use);
sys_labels = {sim_systems(sys_use).label};
sys_shorts = {sim_systems(sys_use).short};
sys_colors = vertcat(sim_systems(sys_use).color);

if is_meg_run
    modality_str = 'MSG (magnetic)';
else
    modality_str = 'ESG (electric)';
end
fprintf('  Modality: %s   Systems: %s\n', modality_str, strjoin(sys_labels, ', '));


% -------------------------------------------------------------------------
% Effective noise s.d. per system per level (bandwidth + trial averaging)
% -------------------------------------------------------------------------
nyquist    = sim_fs / 2;
trial_gain = sqrt(sim_n_trials);
n_lev      = numel(sim_noise_factors);
sigma_mat  = zeros(n_sys, n_lev);
for i = 1:n_sys
    k = sys_use(i);
    bw = min(sim_systems(k).bandwidth_hz, nyquist);
    sigma_mat(i, :) = sim_systems(k).noise_baseline * sim_noise_factors ...
                      * sqrt(bw) / trial_gain;
end

w     = sim_waveform(:);
n_ori = numel(sim_orientations);


% -------------------------------------------------------------------------
% Main loop over perturbation families
% -------------------------------------------------------------------------
results = struct();   % results.(name).self / .perfect etc.

rng(sim_noise_seed);

for p = 1:numel(perts)
    P = perts(p);

    if ~isfield(leadfields, P.ref)
        fprintf('  [%s] reference not loaded — skipping.\n', P.name);
        continue
    end

    % Keep only shift keys that were actually loaded
    present = cellfun(@(k) isfield(leadfields, k), P.shift_keys);
    keys    = P.shift_keys(present);
    b_idx   = P.bundle_idx(present);
    if isempty(keys)
        fprintf('  [%s] no shift leadfields loaded — skipping.\n', P.name);
        continue
    end

    ref_lf   = leadfields.(P.ref);
    n_src    = ref_lf.n_sources;
    n_bundle = max(P.bundle_idx);

    fprintf('  [%s] %d configs, %d sources, %d orientations\n', ...
        P.name, numel(keys), n_src, n_ori);

    % Accumulators: [n_sys x n_ori x n_bundle x n_src x n_lev]
    acc_self  = zeros(n_sys, n_ori, n_bundle, n_src, n_lev);
    acc_perf  = zeros(n_sys, n_ori, n_bundle, n_src, n_lev);
    n_in_bund = zeros(1, n_bundle);

    t0 = tic;
    for ki = 1:numel(keys)
        b     = b_idx(ki);
        lf    = leadfields.(keys{ki});
        n_in_bund(b) = n_in_bund(b) + 1;

        for o = 1:n_ori
            ori = sim_orientations{o};
            for src = 1:n_src
                g = vertcat(lf.(ori){:, src});          % perturbed field
                h = vertcat(ref_lf.(ori){:, src});      % perfect field

                % refs{1} = self (g), refs{2} = perfect (h)
                rq = sim_noise_rsq_refs(g, {g, h}, w, sigma_mat, sim_pert_n_real);
                rq = mean(rq, 4, 'omitnan');            % avg realisations -> [2 x n_sys x n_lev]

                acc_self(:, o, b, src, :) = acc_self(:, o, b, src, :) ...
                    + reshape(rq(1, :, :), [n_sys 1 1 1 n_lev]);
                acc_perf(:, o, b, src, :) = acc_perf(:, o, b, src, :) ...
                    + reshape(rq(2, :, :), [n_sys 1 1 1 n_lev]);
            end
        end
        fprintf('    %-40s (bundle %d)\n', keys{ki}, b);
    end

    % Average across shifts within each bundle
    for b = 1:n_bundle
        if n_in_bund(b) > 0
            acc_self(:, :, b, :, :) = acc_self(:, :, b, :, :) / n_in_bund(b);
            acc_perf(:, :, b, :, :) = acc_perf(:, :, b, :, :) / n_in_bund(b);
        end
    end

    R = struct();
    R.persrc_self  = acc_self;                       % [sys ori bundle src lev]
    R.persrc_perf  = acc_perf;
    % Reshape (not squeeze) to drop only the source dim — squeeze would also
    % collapse the system dim when there is a single system (the ESG run).
    R.cordmean_self = reshape(mean(acc_self, 4, 'omitnan'), [n_sys, n_ori, n_bundle, n_lev]);
    R.cordmean_perf = reshape(mean(acc_perf, 4, 'omitnan'), [n_sys, n_ori, n_bundle, n_lev]);
    R.n_in_bundle   = n_in_bund;
    R.n_bundle      = n_bundle;
    R.n_src         = n_src;
    R.bundle_display = P.bundle_display;
    R.bundle_colors  = P.bundle_colors;
    results.(P.name) = R;

    fprintf('    done (%.1f s)\n', toc(t0));
end


% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
pert_names = fieldnames(results);
if isempty(pert_names)
    error('No perturbation families produced results — nothing to save.');
end
n_src_out = results.(pert_names{1}).n_src;
src_mm    = (0:n_src_out-1) * src_spacing_mm;

outfile = fullfile(sim_out_dir, 'sim_pert_noise.mat');
save(outfile, 'results', 'pert_names', ...
    'sim_noise_factors', 'sim_pert_n_real', 'sim_n_trials', ...
    'sim_orientations', 'sim_ori_display', ...
    'sys_labels', 'sys_shorts', 'sys_colors', ...
    'is_meg_run', 'src_mm', 'src_spacing_mm', ...
    'sim_dipole_nAm', 'sim_freq', 'sim_evoked_latency', '-v7.3');

fprintf('\nSaved: %s\n', outfile);
fprintf('Next: run sim_plot_perturbation_noise\n');

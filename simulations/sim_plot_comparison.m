% sim_plot_comparison - Overlay noise curves across geometry variants
%
% Reads the per-variant results produced by sim_run_geometries and, for each
% sensor system, overlays the r^2-vs-noise curves so the geometry variants can
% be compared directly. One figure per system: rows are perturbation families
% (source / sensor / conductivity), columns are dipole orientation, and within
% each panel the unperturbed baseline is drawn against the small / medium /
% large representative shifts.
%
% USAGE:
%   sim_plot_comparison
%
% OUTPUT (to <sim_save_dir>/comparison/):
%   comparison_<system>.png / .fig     one per system
%   comparison_table.tsv                cord-mean r^2 per variant/system/level
%
% Each curve is the mean r^2 across the whole cord at each noise level. The
% baseline (unperturbed) curve is repeated in every family's panels as the
% reference every shift is compared against.
%
% DEPENDENCIES:
%   config_sim, sim_geometry_index.mat + per-variant sim_noise_<system>.mat
%   (all from sim_run_geometries)
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

idxfile = fullfile(sim_out_dir, 'sim_geometry_index.mat');
if ~isfile(idxfile)
    error('Not found: %s\nRun sim_run_geometries first.', idxfile);
end
load(idxfile);   %#ok<LOAD>

save_dir = fullfile(sim_save_dir, 'comparison');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

nf     = sim_noise_factors;
n_ori  = numel(sim_orientations);
n_lev  = numel(nf);

group_title = containers.Map( ...
    {'source', 'sensor', 'cond'}, ...
    {'Source-position shift', 'Sensor-position shift', 'Conductivity change'});
bundle_order = {'small', 'medium', 'large'};

sys_shorts = unique({index.sys_short}, 'stable');

fprintf('sim_plot_comparison\n');
fprintf('  %d systems, %d variants indexed\n\n', numel(sys_shorts), numel(index));

% Accumulate rows for a summary table
tbl = {sprintf('%-12s\t%-10s\t%-8s\t%-16s\t%-10s\t%-10s', ...
    'System', 'Group', 'Bundle', 'Orientation', 'NoiseFactor', 'MeanRsq')};


% -------------------------------------------------------------------------
% One figure per system
% -------------------------------------------------------------------------
for si = 1:numel(sys_shorts)
    ss    = sys_shorts{si};
    ents  = index(strcmp({index.sys_short}, ss));
    sys_label = ents(1).system;

    % Baseline (unperturbed) for this system, if present
    base_ent = ents(strcmp({ents.group}, 'baseline'));
    if ~isempty(base_ent)
        base_cm = load_cordmean(base_ent(1), n_ori);
    else
        base_cm = [];
    end

    % Perturbation families present for this system
    groups = unique({ents.group}, 'stable');
    groups = groups(~strcmp(groups, 'baseline'));
    if isempty(groups)
        fprintf('  %s: only baseline present — skipping figure.\n', sys_label);
        continue
    end
    n_grp = numel(groups);

    fig = figure('Color', 'w', 'Units', 'inches', ...
        'Position', [1, 1, 4.6*n_ori, 3.2*n_grp]);
    tl = tiledlayout(n_grp, n_ori, 'TileSpacing', 'compact', 'Padding', 'tight');

    for gi = 1:n_grp
        grp      = groups{gi};
        grp_ents = ents(strcmp({ents.group}, grp));
        if group_title.isKey(grp); gtitle = group_title(grp); else; gtitle = grp; end

        for o = 1:n_ori
            ax = nexttile(tl, (gi-1)*n_ori + o);
            hold(ax, 'on');

            leg_h = gobjects(0); leg_t = {};

            % Baseline reference
            if ~isempty(base_cm)
                bcol = sim_bundle_colors('none');
                h0 = plot(ax, nf, base_cm(o, :), '--', 'Color', bcol, ...
                    'LineWidth', pub_line_width);
                leg_h(end+1) = h0; leg_t{end+1} = 'baseline (unperturbed)'; %#ok<SAGROW>
            end

            % Bundle variants in small -> large order
            for bo = 1:numel(bundle_order)
                be = grp_ents(strcmp({grp_ents.bundle}, bundle_order{bo}));
                if isempty(be); continue; end
                cm  = load_cordmean(be(1), n_ori);
                col = sim_bundle_colors(bundle_order{bo});
                hh  = plot(ax, nf, cm(o, :), '-', 'Color', col, ...
                    'LineWidth', pub_line_width, 'Marker', 'o', ...
                    'MarkerSize', pub_marker_size, 'MarkerFaceColor', col);
                leg_h(end+1) = hh; leg_t{end+1} = bundle_order{bo}; %#ok<SAGROW>

                for L = 1:n_lev
                    tbl{end+1} = sprintf('%-12s\t%-10s\t%-8s\t%-16s\t%-10g\t%-10.4f', ...
                        sys_label, grp, bundle_order{bo}, sim_ori_display{o}, ...
                        nf(L), cm(o, L)); %#ok<SAGROW>
                end
            end

            xline(ax, 1, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
            set(ax, 'XScale', 'log', 'XTick', nf, ...
                'XTickLabel', arrayfun(@(f) sprintf('%gx', f), nf, ...
                                       'UniformOutput', false));
            xlabel(ax, 'Sensor noise (\times baseline)');
            ylabel(ax, 'r^2  vs noise-free field');
            ylim(ax, [0, 1.02]); grid(ax, 'on'); box(ax, 'on');

            if gi == 1; title(ax, sim_ori_display{o}, 'FontWeight', 'bold'); end
            if o == 1
                text(-0.30, 0.5, gtitle, 'Units', 'normalized', 'Rotation', 90, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontWeight', 'bold', 'FontSize', 10);
                legend(ax, leg_h, leg_t, 'Location', 'southwest', 'FontSize', 8);
            end
        end
    end

    title(tl, sprintf('Noise robustness across geometry variants  |  %s', sys_label), ...
        'FontSize', 13, 'FontWeight', 'bold');
    subtitle(tl, sprintf(['%s array; evoked %g nA\\cdotm @ %g Hz, %d trials; ' ...
                          'mean r^2 across cord'], ...
        sim_array, sim_dipole_nAm, sim_freq, sim_n_trials), ...
        'FontSize', 9, 'Color', [0.4 0.4 0.4]);

    fname = sprintf('comparison_%s', ss);
    exportgraphics(fig, fullfile(save_dir, [fname '.png']), 'Resolution', 300);
    saveas(fig, fullfile(save_dir, [fname '.fig']));
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

% Summary table
fid = fopen(fullfile(save_dir, 'comparison_table.tsv'), 'w');
fprintf(fid, '%s\n', tbl{:});
fclose(fid);
fprintf('  Saved: comparison_table.tsv\n');

fprintf('\nsim_plot_comparison complete.\nFigures: %s\n', save_dir);


% -------------------------------------------------------------------------
% Local function: cord-mean r^2 [n_ori x n_lev] for one indexed variant
% -------------------------------------------------------------------------
function cm = load_cordmean(entry, n_ori)
    d  = load(entry.file, 'rsq_mean');           % [n_ori x n_src x n_lev]
    cm = squeeze(mean(d.rsq_mean, 2, 'omitnan')); % -> [n_ori x n_lev]
    if size(cm, 1) ~= n_ori                       % guard single-orientation case
        cm = reshape(cm, n_ori, []);
    end
end

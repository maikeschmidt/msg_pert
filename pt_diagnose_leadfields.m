% pt_diagnose_leadfields - Check whether perturbed leadfields actually differ
%                          from the reference (diagnose r² = 100% everywhere)
%
% The heatmaps compute r² directly as corrcoef(reference, perturbed) per source.
% If they read 100% for every perturbation, the loaded perturbed leadfields are
% identical (or proportional) to the reference. This script confirms that by
% loading leadfields_organised.mat for the ACTIVE modality and, for each
% perturbation key, reporting:
%   - max |perturbed - reference|   (0 => byte-for-byte identical data)
%   - median r² vs reference        (1.00 => proportional / no change)
%   - the sensor count per axis     (a tiny count would force corrcoef = 1)
%
% USAGE:
%   pt_modality('set','esg');   % or 'msg'
%   pt_diagnose_leadfields
%
% INTERPRETATION:
%   max-diff = 0            -> the perturbed FILE holds the same data as the
%                             reference: the perturbation was not applied when
%                             the leadfields were generated (a msg_fwd / geometry
%                             issue), or the loader matched the wrong file.
%   max-diff > 0 but r²≈1  -> the field genuinely barely changes (smooth ESG
%                             fields tolerate the shift), OR the shift was too
%                             small — not an analysis bug.
%   sensors/axis <= 2      -> corrcoef of 2-element vectors is always ±1; the
%                             axis split is wrong (check sensor_n_axes).
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

lf_file = fullfile(forward_fields_base, 'leadfields_organised.mat');
if ~isfile(lf_file)
    error('Not found: %s\nRun pt_load_leadfields first (for this modality).', lf_file);
end
S = load(lf_file, 'leadfields');
leadfields = S.leadfields;

fprintf('pt_diagnose_leadfields  |  modality = %s  (%d axes, is_meg=%d)\n', ...
    active_modality, sensor_n_axes, sensor_is_meg);
fprintf('  %s\n\n', lf_file);

% Perturbation families to check: {reference key, {perturbed keys}}
checks = {};
checks{end+1} = {'SOURCE', ['bem_' sensitivity_ref_key], ...
    cellfun(@(k) ['bem_' k], sensitivity_keys, 'UniformOutput', false)};
checks{end+1} = {'SENSOR', ['bem_' sensor_sensitivity_ref_key], ...
    cellfun(@(k) ['bem_' k], sensor_sensitivity_keys, 'UniformOutput', false)};
cond_keys = cell(1, n_cond_total);
for b = 1:n_cond_bundles
    for s = 1:n_cond_shifts
        cond_keys{(b-1)*n_cond_shifts+s} = ...
            sprintf('bem_cond_%s_bundle%d_shift%d', cond_sensitivity_ref_key, b, s);
    end
end
checks{end+1} = {'COND', ['bem_' cond_sensitivity_ref_key], cond_keys};

ori = orientation_labels{1};   % check the first orientation

for c = 1:numel(checks)
    name     = checks{c}{1};
    ref_key  = checks{c}{2};
    pert_keys = checks{c}{3};

    fprintf('=== %s (ref: %s) ===\n', name, ref_key);
    if ~isfield(leadfields, ref_key)
        fprintf('  reference key not loaded — skipping.\n\n');
        continue
    end

    ref = leadfields.(ref_key);
    n_ax  = ref.n_sensor_axes;
    n_src = ref.n_sources;
    mid   = round(n_src/2);
    fprintf('  reference: %d axes, %d sensors/axis, %d sources\n', ...
        n_ax, ref.n_sensors_per_axis, n_src);

    n_present = 0;
    for k = 1:numel(pert_keys)
        pk = pert_keys{k};
        if ~isfield(leadfields, pk); continue; end
        n_present = n_present + 1;

        % Compare at a mid-cord source, axis 1
        va = ref.(ori){1, mid};
        vb = leadfields.(pk).(ori){1, mid};
        md = max(abs(vb - va));
        if numel(va) >= 2
            cc = corrcoef(va, vb); r2 = cc(1,2)^2;
        else
            r2 = NaN;
        end

        % Only print the first few and any that are suspiciously identical
        if k <= 3 || md == 0
            fprintf('    %-45s  maxdiff=%.3g  r²=%.4f  (n=%d)\n', ...
                pk, md, r2, numel(va));
        end
    end
    if n_present == 0
        fprintf('  NO perturbed keys loaded for %s — check paths/naming.\n', name);
    else
        fprintf('  %d/%d perturbed keys present.\n', n_present, numel(pert_keys));
    end
    fprintf('\n');
end

fprintf(['Read the maxdiff column: 0 means the perturbed leadfield FILE is ' ...
         'identical to the\nreference (generation issue), not an analysis bug.\n']);

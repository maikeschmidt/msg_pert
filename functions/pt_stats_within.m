function T = pt_stats_within(vals, bund, mags, opts)
% pt_stats_within - Magnitude and perturbation-type comparisons, one modality
%
% Runs the three families of within-modality comparisons on one metric,
% where each observation is one perturbation realisation (one random shift)
% summarised over the cord. Every row carries an estimate and a bootstrap
% interval; no significance test is applied (see pt_compare_two for why).
%
%   dose_response  per perturbation type, Spearman correlation between the
%                  metric and the perturbation magnitude (|shift| in mm, or
%                  mean % conductivity change), bootstrap interval.
%   bundle_pairs   per perturbation type, small vs medium, medium vs large,
%                  small vs large. Unpaired: different random shifts.
%   type_pairs     per bundle and pooled, every pair of perturbation types.
%                  Paired where opts.paired_types says the realisations match
%                  one-to-one (source and sensor shifts drawn from the same
%                  shift vectors), unpaired otherwise.
%
% USAGE:
%   T = pt_stats_within(vals, bund, mags, opts)
%
% INPUT:
%   vals  - struct, one field per perturbation type: [n x 1] values
%   bund  - struct, same fields: [n x 1] bundle index (1..n_bundles)
%   mags  - struct, same fields: [n x 1] perturbation magnitude
%   opts  - struct:
%             .types         perturbation types, in reporting order
%             .type_display  struct: type -> display name
%             .bundle_names  {1 x n_bundles}
%             .paired_types  {k x 2} pairs of types matched by index
%             .n_boot        bootstrap draws
%             .level         interval coverage, e.g. 0.95
%
% OUTPUT:
%   T     - table: family, subset, comparison, design, n_a, n_b, median_a,
%           median_b, estimate, estimate_name, ci_lo, ci_hi, separated
%           estimate is the median difference (A - B) for comparisons and
%           Spearman rho for dose response; separated = interval excludes 0.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

rows  = {};
types = opts.types(isfield(vals, opts.types));
nb    = numel(opts.bundle_names);

% ---- dose response ------------------------------------------------------
for t = 1:numel(types)
    ty = types{t};
    [rho, ci] = pt_rank_corr(mags.(ty), vals.(ty), opts.n_boot, opts.level);
    ok = ~isnan(vals.(ty)) & ~isnan(mags.(ty));
    rows(end+1, :) = {'dose_response', opts.type_display.(ty), ...
        'metric vs magnitude', 'rank correlation', sum(ok), sum(ok), ...
        NaN, NaN, rho, 'Spearman rho', ci(1), ci(2), ...
        ~isnan(ci(1)) && (ci(1) > 0 || ci(2) < 0)}; %#ok<AGROW>
end

% ---- bundle pairs -------------------------------------------------------
pairs = nchoosek(1:nb, 2);
for t = 1:numel(types)
    ty = types{t};
    for k = 1:size(pairs, 1)
        b1 = pairs(k, 1); b2 = pairs(k, 2);
        R = pt_compare_two(vals.(ty)(bund.(ty) == b2), ...
                           vals.(ty)(bund.(ty) == b1), false, opts.n_boot, opts.level);
        rows(end+1, :) = pt_compare_row(R, 'bundle_pairs', opts.type_display.(ty), ...
            sprintf('%s vs %s', opts.bundle_names{b2}, opts.bundle_names{b1})); %#ok<AGROW>
    end
end

% ---- type pairs ---------------------------------------------------------
if numel(types) >= 2
    tp = nchoosek(1:numel(types), 2);
    subsets = [num2cell(1:nb), {0}];          % 0 = pooled over bundles
    for si = 1:numel(subsets)
        b = subsets{si};
        if b == 0, sub_name = 'pooled'; else, sub_name = opts.bundle_names{b}; end
        for k = 1:size(tp, 1)
            ta = types{tp(k, 1)}; tb = types{tp(k, 2)};
            if b == 0
                ia = true(size(vals.(ta))); ib = true(size(vals.(tb)));
            else
                ia = bund.(ta) == b;        ib = bund.(tb) == b;
            end
            paired = is_paired(ta, tb, opts.paired_types) && sum(ia) == sum(ib);
            R = pt_compare_two(vals.(ta)(ia), vals.(tb)(ib), paired, opts.n_boot, opts.level);
            rows(end+1, :) = pt_compare_row(R, 'type_pairs', sub_name, ...
                sprintf('%s vs %s', opts.type_display.(ta), ...
                opts.type_display.(tb))); %#ok<AGROW>
        end
    end
end

T = cell2table(rows, 'VariableNames', pt_compare_row());
end


function tf = is_paired(a, b, pairs)
    tf = false;
    for k = 1:size(pairs, 1)
        if (strcmp(pairs{k,1}, a) && strcmp(pairs{k,2}, b)) || ...
           (strcmp(pairs{k,1}, b) && strcmp(pairs{k,2}, a))
            tf = true;
            return
        end
    end
end

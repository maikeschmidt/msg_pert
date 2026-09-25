function R = pt_compare_two(x, y, paired, n_perm)
% pt_compare_two - One two-condition comparison: test, effect size, medians
%
% The single place msg_pert decides how two sets of per-realisation values
% are compared, so every stage reports the same statistics.
%
%   paired   sign-flip permutation test on x - y (st_signflip_test),
%            effect size = matched-pairs rank-biserial (st_rank_biserial),
%            median difference = median(x - y). Pairs with a NaN are dropped.
%   unpaired label-shuffling permutation test (pt_perm_test2),
%            effect size = Cliff's delta (pt_cliffs_delta),
%            median difference = median(x) - median(y).
%
% Both effect sizes lie in [-1, 1] and are positive when x tends to be
% larger than y, so they can share one colour scale in a figure.
%
% USAGE:
%   R = pt_compare_two(x, y, paired)
%   R = pt_compare_two(x, y, paired, n_perm)
%
% INPUT:
%   x, y    - per-realisation values of the two conditions
%   paired  - true when x(i) and y(i) come from the same realisation
%   n_perm  - permutations (default 10000)
%
% OUTPUT:
%   R       - struct: .test .n_a .n_b .median_a .median_b .median_diff
%                     .effect .effect_name .p   (two-sided)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 4 || isempty(n_perm), n_perm = 10000; end

x = x(:); y = y(:);

if paired
    if numel(x) ~= numel(y)
        error('pt_compare_two:pairLength', ...
            'Paired comparison needs equal lengths (%d vs %d).', numel(x), numel(y));
    end
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);
    d = x - y;
    R.test        = 'sign-flip (paired)';
    R.n_a         = numel(x);
    R.n_b         = numel(y);
    R.median_a    = median_or_nan(x);
    R.median_b    = median_or_nan(y);
    R.median_diff = median_or_nan(d);
    R.effect      = st_rank_biserial(d);
    R.effect_name = 'rank-biserial';
    R.p           = st_signflip_test(d, n_perm, 'both');
else
    x = x(~isnan(x)); y = y(~isnan(y));
    R.test        = 'permutation (unpaired)';
    R.n_a         = numel(x);
    R.n_b         = numel(y);
    R.median_a    = median_or_nan(x);
    R.median_b    = median_or_nan(y);
    R.median_diff = R.median_a - R.median_b;
    R.effect      = pt_cliffs_delta(x, y);
    R.effect_name = 'Cliff delta';
    R.p           = pt_perm_test2(x, y, max(n_perm, 12870), 'both');
end
end


function m = median_or_nan(v)
    if isempty(v), m = NaN; else, m = median(v); end
end

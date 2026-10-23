function R = pt_compare_two(x, y, paired, n_boot, level)
% pt_compare_two - One two-condition comparison: medians and an interval
%
% The single place msg_pert decides how two sets of per-realisation values
% are compared, so every stage reports the same numbers.
%
% Every value here is computed exactly from a simulated geometry, so there is
% no measurement noise and no null hypothesis for a significance test to
% test. What is uncertain is how well the random draws represent the
% perturbations that could occur, which is a question about spread, answered
% with an interval (the msg_fwd convention, see st_warp_summary).
%
%   paired   difference = median(x - y), interval = percentile bootstrap of
%            that median, resampling realisations (st_boot_ci_median).
%            Pairs with a NaN are dropped.
%   unpaired difference = median(x) - median(y), interval = percentile
%            bootstrap resampling each condition independently.
%
% The two conditions are called separated when the interval of the
% difference excludes zero. That is a description of the interval, not a
% test: it says the difference is larger than the realisation-to-realisation
% spread would produce by resampling.
%
% USAGE:
%   R = pt_compare_two(x, y, paired)
%   R = pt_compare_two(x, y, paired, n_boot, level)
%
% INPUT:
%   x, y    - per-realisation values of the two conditions
%   paired  - true when x(i) and y(i) come from the same realisation
%   n_boot  - bootstrap draws (default 10000)
%   level   - interval coverage (default 0.95)
%
% OUTPUT:
%   R       - struct: .design .n_a .n_b .median_a .median_b .median_diff
%                     .ci_lo .ci_hi .separated
%             The interval is NaN with fewer than 3 values in either
%             condition, and separated is then false.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 4 || isempty(n_boot), n_boot = 10000; end
if nargin < 5 || isempty(level),  level  = 0.95;  end

x = x(:); y = y(:);

if paired
    if numel(x) ~= numel(y)
        error('pt_compare_two:pairLength', ...
            'Paired comparison needs equal lengths (%d vs %d).', numel(x), numel(y));
    end
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);
    d = x - y;
    R.design      = 'paired';
    R.n_a         = numel(x);
    R.n_b         = numel(y);
    R.median_a    = median_or_nan(x);
    R.median_b    = median_or_nan(y);
    R.median_diff = median_or_nan(d);
    ci            = st_boot_ci_median(d, n_boot, level);
else
    x = x(~isnan(x)); y = y(~isnan(y));
    R.design      = 'unpaired';
    R.n_a         = numel(x);
    R.n_b         = numel(y);
    R.median_a    = median_or_nan(x);
    R.median_b    = median_or_nan(y);
    R.median_diff = R.median_a - R.median_b;
    ci            = boot_ci_diff(x, y, n_boot, level);
end
R.ci_lo     = ci(1);
R.ci_hi     = ci(2);
R.separated = ~isnan(ci(1)) && (ci(1) > 0 || ci(2) < 0);
end


function m = median_or_nan(v)
    if isempty(v), m = NaN; else, m = median(v); end
end

function ci = boot_ci_diff(x, y, n_boot, level)
% Percentile bootstrap of median(x) - median(y), the two resampled
% independently. Same percentile convention as st_boot_ci_median (pctl).
    ci = [NaN; NaN];
    na = numel(x); nb = numel(y);
    if na < 3 || nb < 3, return; end
    d = median(x(randi(na, na, n_boot)), 1) - median(y(randi(nb, nb, n_boot)), 1);
    a = (1 - level) / 2;
    ci = [pctl(d, a * 100); pctl(d, (1 - a) * 100)];
end

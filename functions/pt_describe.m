function D = pt_describe(v, n_boot, level, cover_pct)
% pt_describe - Median, IQR, coverage percentile and their intervals
%
% The same summary msg_fwd reports for its replicate geometries
% (st_warp_summary): the median with a bootstrap interval, the IQR, and a
% coverage percentile with its own bootstrap interval, so a msg_pert number
% can be read alongside a msg_fwd one. Percentiles use msg_fwd's pctl.
%
% USAGE:
%   D = pt_describe(v)
%   D = pt_describe(v, n_boot, level, cover_pct)
%
% INPUT:
%   v         - values (NaNs removed). In msg_pert these are usually one
%               value per perturbation realisation, so the intervals answer
%               "how much would this number move for a different set of
%               random draws of this size" — not a between-subject interval.
%   n_boot    - bootstrap draws (default 10000)
%   level     - interval coverage (default 0.95)
%   cover_pct - coverage percentile (default 95)
%
% OUTPUT:
%   D         - struct: .n .median .ci_lo .ci_hi .q1 .q3 .cover
%               .cover_ci_lo .cover_ci_hi .min .max
%               Intervals are NaN with fewer than 3 values.
%
% READING THE COVERAGE PERCENTILE
%   Worded descriptively: "in 95% of the realisations, RE was at or below
%   X% (95% CI ...)". It is not a predictive probability for a new
%   perturbation, and with 8 realisations per bundle it rests on the largest
%   one or two values, so quote the interval with it.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 2 || isempty(n_boot),    n_boot    = 10000; end
if nargin < 3 || isempty(level),     level     = 0.95;  end
if nargin < 4 || isempty(cover_pct), cover_pct = 95;    end

v = v(:); v = v(~isnan(v));
D.n = numel(v);
if D.n == 0
    [D.median, D.ci_lo, D.ci_hi, D.q1, D.q3, D.cover, D.cover_ci_lo, ...
        D.cover_ci_hi, D.min, D.max] = deal(NaN);
    return
end
D.median = median(v);
ci = st_boot_ci_median(v, n_boot, level);
D.ci_lo = ci(1);
D.ci_hi = ci(2);
D.q1    = pctl(v, 25);
D.q3    = pctl(v, 75);
D.cover = pctl(v, cover_pct);
[D.cover_ci_lo, D.cover_ci_hi] = deal(NaN);
if D.n >= 3
    % pctl's interpolation rule applied to every resample at once: with no
    % NaNs every resample has n values, so the interpolation point is shared
    S   = sort(v(randi(D.n, D.n, n_boot)), 1);
    pos = max(1, min(D.n, (cover_pct/100) * D.n + 0.5));
    lo  = floor(pos); hi = ceil(pos);
    b   = S(lo, :) + (pos - lo) * (S(hi, :) - S(lo, :));
    a = (1 - level) / 2;
    D.cover_ci_lo = pctl(b, a * 100);
    D.cover_ci_hi = pctl(b, (1 - a) * 100);
end
D.min = min(v);
D.max = max(v);
end

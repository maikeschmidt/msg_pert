function D = pt_describe(v, n_boot)
% pt_describe - Median, bootstrap CI, IQR and range of one set of values
%
% USAGE:
%   D = pt_describe(v)
%   D = pt_describe(v, n_boot)
%
% INPUT:
%   v       - values (NaNs removed). In msg_pert these are usually one value
%             per perturbation realisation, so the CI answers "how much would
%             the median move for a different set of random shifts of this
%             size" — not a between-subject interval.
%   n_boot  - bootstrap draws for the 95% percentile CI (default 10000)
%
% OUTPUT:
%   D       - struct: .n .median .ci_lo .ci_hi .q1 .q3 .min .max
%             CI is NaN with fewer than 3 values (st_boot_ci_median).
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 2 || isempty(n_boot), n_boot = 10000; end

v = v(:); v = v(~isnan(v));
D.n = numel(v);
if D.n == 0
    [D.median, D.ci_lo, D.ci_hi, D.q1, D.q3, D.min, D.max] = deal(NaN);
    return
end
D.median = median(v);
ci = st_boot_ci_median(v, n_boot, 0.95);
D.ci_lo = ci(1);
D.ci_hi = ci(2);
q = pt_quantile(v, [25 75], 1);
D.q1  = q(1);
D.q3  = q(2);
D.min = min(v);
D.max = max(v);
end

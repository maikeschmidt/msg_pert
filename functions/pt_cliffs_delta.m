function d = pt_cliffs_delta(x, y)
% pt_cliffs_delta - Cliff's delta effect size for two independent groups
%
% The probability that a value from x exceeds one from y, minus the reverse:
%   delta = P(x > y) - P(x < y)
% It is the independent-groups counterpart of the matched-pairs rank-biserial
% correlation (st_rank_biserial) and is read on the same scale.
%
% USAGE:
%   d = pt_cliffs_delta(x, y)
%
% INPUT:
%   x, y  - observations (NaNs removed)
%
% OUTPUT:
%   d     - Cliff's delta in [-1, 1]; +1 when every x exceeds every y.
%           Rough guides: |d| < 0.147 negligible, < 0.33 small, < 0.474
%           medium, otherwise large (Romano et al. 2006). NaN if a group is
%           empty.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

x = x(:); x = x(~isnan(x));
y = y(:); y = y(~isnan(y));
if isempty(x) || isempty(y)
    d = NaN;
    return
end
s = sign(x - y');
d = mean(s(:));
end

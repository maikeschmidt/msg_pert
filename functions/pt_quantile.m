function q = pt_quantile(x, p, dim)
% pt_quantile - Percentiles without the Statistics Toolbox
%
% Linear-interpolation percentiles matching MATLAB's prctile convention,
% ignoring NaNs, along one dimension.
%
% USAGE:
%   q = pt_quantile(x, p)          % along the first non-singleton dimension
%   q = pt_quantile(x, p, dim)
%
% INPUT:
%   x    - numeric array
%   p    - percentile(s) in [0, 100]
%   dim  - dimension to work along
%
% OUTPUT:
%   q    - x with dimension dim replaced by numel(p) percentiles. All-NaN
%          slices give NaN.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 3
    dim = find(size(x) > 1, 1);
    if isempty(dim), dim = 1; end
end

nd   = max(ndims(x), dim);
perm = [dim, setdiff(1:nd, dim)];
xp   = permute(x, perm);
sz   = size(xp);
xp   = reshape(xp, sz(1), []);

q = nan(numel(p), size(xp, 2));
for c = 1:size(xp, 2)
    v = sort(xp(~isnan(xp(:, c)), c));
    n = numel(v);
    if n == 0, continue; end
    if n == 1, q(:, c) = v; continue; end
    for k = 1:numel(p)
        pos = p(k) / 100 * n + 0.5;
        pos = max(1, min(n, pos));
        lo  = floor(pos); hi = ceil(pos);
        q(k, c) = v(lo) + (pos - lo) * (v(hi) - v(lo));
    end
end

sz(1) = numel(p);
q = ipermute(reshape(q, sz), perm);
end

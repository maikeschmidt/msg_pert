function q = pt_quantile(x, p, dim)
% pt_quantile - Percentiles of an array along one dimension
%
% Applies msg_fwd's pctl to every slice, so every percentile in msg_pert is
% computed by the same definition as every percentile in msg_fwd and
% msg_models (linear interpolation matching MATLAB's prctile, NaNs ignored).
% This function only adds the array handling pctl leaves to its callers.
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
    for k = 1:numel(p)
        q(k, c) = pctl(xp(:, c), p(k));
    end
end

sz(1) = numel(p);
q = ipermute(reshape(q, sz), perm);
end

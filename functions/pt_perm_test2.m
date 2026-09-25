function [p, obs] = pt_perm_test2(x, y, n_perm, tail)
% pt_perm_test2 - Two-sample (independent groups) permutation test
%
% Tests whether two independent groups differ in mean, by reshuffling the
% group labels. The unpaired counterpart of msg_fwd's st_signflip_test, with
% the same conventions: mean-difference statistic, exact enumeration when the
% number of relabellings is small, otherwise Monte Carlo with the Phipson &
% Smyth (2010) +1 correction so p is never reported as exactly zero.
%
% USAGE:
%   p = pt_perm_test2(x, y)
%   [p, obs] = pt_perm_test2(x, y, n_perm, tail)
%
% INPUT:
%   x, y    - observations of the two groups (NaNs removed)
%   n_perm  - Monte Carlo draws, and the enumeration limit (default 10000):
%             with nchoosek(nx+ny, nx) <= n_perm every relabelling is used
%             and p is exact. Two groups of 8 need 12870, so pass >= 12870
%             for exact bundle-vs-bundle tests.
%   tail    - 'both' (default) | 'right' (H1 mean(x) > mean(y)) | 'left'
%
% OUTPUT:
%   p       - permutation p-value (NaN with fewer than 3 per group)
%   obs     - observed mean(x) - mean(y)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 3 || isempty(n_perm), n_perm = 10000; end
if nargin < 4 || isempty(tail),   tail   = 'both'; end

x = x(:); x = x(~isnan(x));
y = y(:); y = y(~isnan(y));
nx = numel(x); ny = numel(y);

p = NaN; obs = NaN;
if nx < 3 || ny < 3, return; end

obs = mean(x) - mean(y);
z   = [x; y];
n   = nx + ny;
tot = sum(z);

n_exact = nchoosek_safe(n, nx);
if n_exact <= n_perm
    idx   = nchoosek(1:n, nx);
    sx    = sum(reshape(z(idx), size(idx)), 2);
    null  = sx / nx - (tot - sx) / ny;
    exact = true;
else
    [~, ord] = sort(rand(n_perm, n), 2);     % one random relabelling per row
    sx    = sum(reshape(z(ord(:, 1:nx)), n_perm, nx), 2);
    null  = sx / nx - (tot - sx) / ny;
    exact = false;
end

tol = eps(max(1, abs(obs))) * 10;
switch lower(tail)
    case 'both',  cnt = sum(abs(null) >= abs(obs) - tol);
    case 'right', cnt = sum(null >= obs - tol);
    case 'left',  cnt = sum(null <= obs + tol);
    otherwise
        error('pt_perm_test2:badTail', 'Unknown tail "%s".', tail);
end

if exact
    p = cnt / numel(null);
else
    p = (cnt + 1) / (numel(null) + 1);
end
p = min(p, 1);
end


function c = nchoosek_safe(n, k)
% Number of relabellings without the precision warning nchoosek gives for
% large n; anything above ~1e15 is simply "too many to enumerate".
    c = exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1));
    c = round(c);
end

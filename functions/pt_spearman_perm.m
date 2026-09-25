function [rho, p] = pt_spearman_perm(x, y, n_perm)
% pt_spearman_perm - Spearman rank correlation with a permutation p-value
%
% Used for dose-response questions: does a metric grow monotonically with
% the size of the perturbation? Ranks are tie-averaged, the correlation is
% Pearson's r on the ranks, and the two-sided p-value comes from shuffling
% y against x. No Statistics Toolbox is needed.
%
% USAGE:
%   [rho, p] = pt_spearman_perm(x, y)
%   [rho, p] = pt_spearman_perm(x, y, n_perm)
%
% INPUT:
%   x, y    - paired observations (pairs with any NaN are removed)
%   n_perm  - permutations (default 10000)
%
% OUTPUT:
%   rho     - Spearman rank correlation
%   p       - two-sided permutation p-value, Phipson & Smyth corrected
%             (NaN with fewer than 4 pairs)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 3 || isempty(n_perm), n_perm = 10000; end

x = x(:); y = y(:);
ok = ~isnan(x) & ~isnan(y);
x = x(ok); y = y(ok);

rho = NaN; p = NaN;
n = numel(x);
if n < 4, return; end

rx = tied_ranks(x) - (n + 1) / 2;
ry = tied_ranks(y) - (n + 1) / 2;
den = sqrt(sum(rx.^2) * sum(ry.^2));
if den == 0, return; end

rho = sum(rx .* ry) / den;

[~, ord] = sort(rand(n, n_perm), 1);       % one shuffle per column
null = (rx' * ry(ord)) / den;
cnt  = sum(abs(null) >= abs(rho) - 1e-12);
p    = (cnt + 1) / (n_perm + 1);
end


function r = tied_ranks(a)
% Ranks with tied values sharing their mean rank.
    n = numel(a);
    [sa, ord] = sort(a);
    rs = (1:n)';
    i = 1;
    while i <= n
        j = i;
        while j < n && sa(j + 1) == sa(i)
            j = j + 1;
        end
        if j > i, rs(i:j) = mean(i:j); end
        i = j + 1;
    end
    r = zeros(n, 1);
    r(ord) = rs;
end

function [rho, ci] = pt_rank_corr(x, y, n_boot, level)
% pt_rank_corr - Spearman rank correlation with a bootstrap interval
%
% Used for dose-response questions: does a metric grow monotonically with
% the size of the perturbation? Ranks are tie-averaged and the correlation
% is Pearson's r on the ranks. The interval is a percentile bootstrap over
% pairs (one pair = one perturbation realisation), so it says how much the
% correlation would move for a different set of random draws of the same
% sizes. No Statistics Toolbox is needed.
%
% USAGE:
%   [rho, ci] = pt_rank_corr(x, y)
%   [rho, ci] = pt_rank_corr(x, y, n_boot, level)
%
% INPUT:
%   x, y    - paired observations (pairs with any NaN are removed)
%   n_boot  - bootstrap draws (default 10000)
%   level   - interval coverage (default 0.95)
%
% OUTPUT:
%   rho     - Spearman rank correlation (NaN with fewer than 4 pairs)
%   ci      - [lo; hi] percentile bootstrap interval. Resamples in which
%             either variable is constant are dropped, since the correlation
%             is undefined there.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 3 || isempty(n_boot), n_boot = 10000; end
if nargin < 4 || isempty(level),  level  = 0.95;  end

x = x(:); y = y(:);
ok = ~isnan(x) & ~isnan(y);
x = x(ok); y = y(ok);

rho = NaN; ci = [NaN; NaN];
n = numel(x);
if n < 4, return; end

rho = spearman(x, y);
if isnan(rho), return; end

b = nan(n_boot, 1);
for k = 1:n_boot
    i    = randi(n, n, 1);
    b(k) = spearman(x(i), y(i));
end
a  = (1 - level) / 2;
ci = [pctl(b, a * 100); pctl(b, (1 - a) * 100)];
end


function r = spearman(x, y)
    n  = numel(x);
    rx = tied_ranks(x) - (n + 1) / 2;
    ry = tied_ranks(y) - (n + 1) / 2;
    den = sqrt(sum(rx.^2) * sum(ry.^2));
    if den == 0, r = NaN; else, r = sum(rx .* ry) / den; end
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

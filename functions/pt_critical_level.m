function [crit, txt] = pt_critical_level(levels, rate, target)
% pt_critical_level - Highest noise level at which a detection rate holds
%
% The rule stages 5 and 7 share: the largest noise level (a multiple of the
% system's own floor) at which the detection rate still reaches target,
% log-interpolated between the simulated levels. Level 0 is excluded.
%
% USAGE:
%   [crit, txt] = pt_critical_level(levels, rate, target)
%
% INPUT:
%   levels  - [1 x n] noise levels, x floor (0 allowed, ignored)
%   rate    - [1 x n] detection rate at each level
%   target  - required rate, e.g. 0.8
%
% OUTPUT:
%   crit    - the critical level; Inf if still detected at the largest
%             level, 0 if not detected even at the smallest non-zero one
%   txt     - the same as text ('> 8', '< 0.125' or the value)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

nz = levels > 0;
lv = levels(nz); rt = rate(nz);
if all(rt >= target)
    crit = Inf; txt = sprintf('> %g', lv(end));
elseif rt(1) < target
    crit = 0;   txt = sprintf('< %g', lv(1));
else
    i = find(rt >= target, 1, 'last');
    if i == numel(lv)
        crit = lv(end);
    else
        f = (rt(i) - target) / (rt(i) - rt(i+1));
        crit = exp(log(lv(i)) + f * (log(lv(i+1)) - log(lv(i))));
    end
    txt = sprintf('%.3g', crit);
end
end

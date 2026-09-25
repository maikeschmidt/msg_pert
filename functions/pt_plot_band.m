function h = pt_plot_band(ax, x, Y, color, style, label)
% pt_plot_band - Median line with an interquartile band
%
% USAGE:
%   h = pt_plot_band(ax, x, Y, color)
%   h = pt_plot_band(ax, x, Y, color, style, label)
%
% INPUT:
%   ax     - target axes
%   x      - [1 x n] x values
%   Y      - [n_rep x n] replicates (rows) at each x; the band is the IQR
%            across rows, the line the median. A single row draws the line
%            only.
%   color  - [1 x 3] RGB
%   style  - line style (default '-')
%   label  - legend text (default '')
%
% OUTPUT:
%   h      - line handle (the band is excluded from legends)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 5 || isempty(style), style = '-'; end
if nargin < 6, label = ''; end

hold(ax, 'on');
x = x(:)';

if size(Y, 1) > 1
    q  = pt_quantile(Y, [25 50 75], 1);
    lo = q(1, :); md = q(2, :); hi = q(3, :);
    ok = ~isnan(lo) & ~isnan(hi);
    if any(ok)
        % Draw each contiguous run separately so NaN gaps stay gaps
        runs = find_runs(ok);
        for k = 1:size(runs, 1)
            i = runs(k, 1):runs(k, 2);
            fill(ax, [x(i), fliplr(x(i))], [lo(i), fliplr(hi(i))], color, ...
                'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
else
    md = Y;
end

h = plot(ax, x, md, style, 'Color', color, 'LineWidth', 1.8, ...
    'DisplayName', label);
end


function runs = find_runs(ok)
    d = diff([false, ok(:)', false]);
    runs = [find(d == 1)', find(d == -1)' - 1];
end

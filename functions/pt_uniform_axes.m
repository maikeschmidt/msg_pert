function pt_uniform_axes(fig)
% pt_uniform_axes - The same axis limits across the panels of a figure
%
% So that panels can be compared by eye: in every tiled layout of the
% figure,
%   - every row shares one y range (rows hold one quantity: one metric, one
%     model),
%   - every column shares one x range (columns hold one condition, e.g. a
%     perturbation type whose size is in mm or in %),
%   - and when every labelled panel carries the same y label (or x label),
%     that axis is shared across the whole figure.
% Each shared range is the union of the panels' own ranges, so nothing is
% cut off.
%
% Called by pt_save_figure on every figure before it is written. Left alone:
%   - hidden axes (label tiles)
%   - axes showing an image, surface or contour, or with a fixed aspect
%     ratio (topoplots, sensor maps, heat maps), whose limits are geometry
%   - axes whose Tag is 'pt_own_limits', for a panel that must keep its own
%
% USAGE:
%   pt_uniform_axes(fig)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

for tl = findobj(fig, 'Type', 'tiledlayout')'
    gs = tl.GridSize;
    A = findobj(tl, '-depth', 1, 'Type', 'axes');
    A = A(arrayfun(@eligible, A));
    if numel(A) < 2, continue; end

    tile = arrayfun(@(a) tile_of(a), A);
    keep = ~isnan(tile);
    A = A(keep); tile = tile(keep);
    if numel(A) < 2, continue; end
    row = ceil(tile / gs(2));
    col = mod(tile - 1, gs(2)) + 1;

    ylab = arrayfun(@(a) strjoin(cellstr(a.YLabel.String), ' '), A, 'UniformOutput', false);
    xlab = arrayfun(@(a) strjoin(cellstr(a.XLabel.String), ' '), A, 'UniformOutput', false);
    same_y = all_same(ylab);
    same_x = all_same(xlab);

    % y: per row, or the whole figure when every labelled panel agrees
    if same_y, groups = ones(size(row)); else, groups = row; end
    for g = unique(groups)'
        share(A(groups == g), 'YLim', 'YScale');
    end
    % x: per column, or the whole figure
    if same_x, groups = ones(size(col)); else, groups = col; end
    for g = unique(groups)'
        share(A(groups == g), 'XLim', 'XScale');
    end
end
end


function tf = eligible(ax)
    tf = strcmp(ax.Visible, 'on') && ~strcmp(ax.Tag, 'pt_own_limits') && ...
         strcmp(ax.DataAspectRatioMode, 'auto') && ...
         isempty(findobj(ax, '-depth', 1, 'Type', 'image', '-or', 'Type', 'surface', ...
                         '-or', 'Type', 'contour')) && ...
         isa(ax.YAxis(1), 'matlab.graphics.axis.decorator.NumericRuler') && ...
         isa(ax.XAxis(1), 'matlab.graphics.axis.decorator.NumericRuler');
end

function t = tile_of(ax)
    t = NaN;
    try
        if isnumeric(ax.Layout.Tile) && all(ax.Layout.TileSpan == 1)
            t = ax.Layout.Tile;
        end
    catch
    end
end

function tf = all_same(lbl)
% Every non-empty label identical (and at least two of them)
    lbl = lbl(~cellfun(@isempty, lbl));
    tf = numel(lbl) >= 2 && numel(unique(lbl)) == 1;
end

function share(A, lim, scale)
% Union of the axes' current limits, applied to all of them
    if numel(A) < 2, return; end
    if numel(unique(arrayfun(@(a) string(a.(scale)), A))) > 1, return; end
    L = cell2mat(arrayfun(@(a) a.(lim), A, 'UniformOutput', false));
    L = [min(L(:, 1)), max(L(:, 2))];
    if all(isfinite(L)) && L(2) > L(1)
        set(A, lim, L);
    end
end

function pt_distinct_markers(fig)
% pt_distinct_markers - A different marker shape for every series in a plot
%
% So that every figure can be read in black and white or with colour-vision
% deficiency: within each axes, every named series (a line or a set of
% points with a legend entry) gets its own marker shape, in the order the
% series were drawn, from pt_marker. Line style and colour are kept.
%
% Called by pt_save_figure on every figure before it is written, so no plot
% has to set markers itself. What it leaves alone:
%   - series with no DisplayName or with HandleVisibility 'off' (reference
%     lines, fitted lines, interval bars, marks that carry their own
%     meaning such as filled = interval excludes zero)
%   - axes with a single named series, or with more than the number of
%     distinct markers (e.g. one line per sensor channel)
%
% On a line with many points the marker is drawn at about 10 positions,
% offset between series so markers of neighbouring lines do not stack.
%
% USAGE:
%   pt_distinct_markers(fig)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

n_mk = numel(pt_marker());
for ax = findobj(fig, 'Type', 'axes')'
    kids = flipud(ax.Children);                  % creation order
    is_series = arrayfun(@(h) any(strcmp(h.Type, {'line', 'scatter'})) && ...
        strcmp(h.HandleVisibility, 'on') && ~isempty(h.DisplayName), kids);
    S = kids(is_series);
    if numel(S) < 2 || numel(S) > n_mk, continue; end
    for k = 1:numel(S)
        h  = S(k);
        mk = pt_marker(k);
        h.Marker = mk;
        if strcmp(h.Type, 'line')
            n = numel(h.XData);
            if n > 15
                % ~10 markers per line, staggered between series
                step = max(1, floor(n / 10));
                off  = 1 + mod((k - 1) * max(1, floor(step / numel(S))), step);
                h.MarkerIndices = off:step:n;
            end
            if isempty(h.MarkerSize) || h.MarkerSize < 5, h.MarkerSize = 6; end
            if ismember(mk, {'x', '+', '*', '.'})
                h.MarkerFaceColor = 'none';
            elseif ischar(h.MarkerFaceColor) && strcmp(h.MarkerFaceColor, 'none')
                h.MarkerFaceColor = 'w';     % hollow, readable over other lines
            end
        end
    end
end
end

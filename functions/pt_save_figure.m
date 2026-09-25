function pt_save_figure(fig, save_dir, name, resolution)
% pt_save_figure - Export a figure as .fig and .png, then close it
%
% USAGE:
%   pt_save_figure(fig, save_dir, name)
%   pt_save_figure(fig, save_dir, name, resolution)
%
% INPUT:
%   fig         - figure handle
%   save_dir    - output folder (created if missing)
%   name        - file stem, no extension
%   resolution  - PNG resolution in dpi (default 600)
%
% NOTES:
%   A failed export is reported as a warning rather than an error, so one
%   figure that cannot be written does not abandon the rest of a stage and
%   its tables. The .fig is written first because it needs the figure
%   handle to still be valid.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 4 || isempty(resolution), resolution = 600; end
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

ok = true;
try
    savefig(fig, fullfile(save_dir, [name '.fig']));
catch err
    ok = false;
    warning('pt_save_figure:fig', 'Could not write %s.fig: %s', name, err.message);
end
try
    exportgraphics(fig, fullfile(save_dir, [name '.png']), 'Resolution', resolution);
catch err
    ok = false;
    warning('pt_save_figure:png', 'Could not write %s.png: %s', name, err.message);
end
if isgraphics(fig), close(fig); end
if ok, fprintf('    saved %s\n', name); end
end

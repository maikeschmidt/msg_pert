function pos_by_axis = sim_sensor_positions(geom_file, array_side, is_meg)
% sim_sensor_positions - Sensor positions for one array, split by sensor axis
%
% Mirrors the field lookup used by msg_fwd/plot_topoplots so the sim topoplots
% place sensors identically to the forward-model figures.
%
% USAGE:
%   pos_by_axis = sim_sensor_positions(geom_file, array_side, is_meg)
%
% INPUT:
%   geom_file  - full path to the geometry .mat
%   array_side - 'front' or 'back'
%   is_meg     - true  -> reads <side>_coils_3axis.chanpos  (3 axes)
%                false -> reads <side>_coils_2axis.elecpos  (2 axes)
%
% OUTPUT:
%   pos_by_axis - {1 x n_axes} cell, each [n_sens x 3] positions
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if ~isfile(geom_file)
    error('sim_sensor_positions: geometry file not found:\n  %s', geom_file);
end
geom = load(geom_file);

if is_meg
    fld = [array_side '_coils_3axis'];
    if ~isfield(geom, fld)
        error('sim_sensor_positions: field "%s" not found in:\n  %s', fld, geom_file);
    end
    pos     = geom.(fld).chanpos;
    n_axes  = 3;
else
    fld = [array_side '_coils_2axis'];
    if ~isfield(geom, fld)
        error('sim_sensor_positions: field "%s" not found in:\n  %s', fld, geom_file);
    end
    pos     = geom.(fld).elecpos;
    n_axes  = 2;
end

n_per_axis  = size(pos, 1) / n_axes;
pos_by_axis = cell(1, n_axes);
for ax = 1:n_axes
    idx1 = (ax-1)*n_per_axis + 1;
    idx2 =  ax   *n_per_axis;
    pos_by_axis{ax} = pos(idx1:idx2, :);
end
end

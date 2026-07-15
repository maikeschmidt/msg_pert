function pos_by_axis = sim_sensor_positions(geom_file, array_side, is_meg, n_axes)
% sim_sensor_positions - Sensor positions for one array, split by sensor axis
%
% MSG and ESG geometries differ ONLY in their sensor definitions — the source
% grid and meshes are identical — which is why they sit in separate folders.
% This function pulls the sensor positions out of whichever geometry file it is
% given, without assuming a single naming convention.
%
% USAGE:
%   pos_by_axis = sim_sensor_positions(geom_file, array_side, is_meg, n_axes)
%
% INPUT:
%   geom_file  - full path to the geometry .mat
%   array_side - 'front' or 'back'
%   is_meg     - true = magnetic array, false = electrode array
%   n_axes     - number of sensor axes to split into. Pass the value from the
%                matching leadfield (lf.n_sensor_axes) so the positions are
%                guaranteed to line up with the leadfield channel order rather
%                than relying on the field name to imply it.
%
% OUTPUT:
%   pos_by_axis - {1 x n_axes} cell, each [n_sens x 3] positions
%
% FIELD LOOKUP:
%   Tries, in order: <side>_coils_3axis, <side>_coils_2axis, <side>_sensors —
%   the same set of names run_bem_leadfields and run_conductivity_perturbation
%   accept. Within the chosen struct, reads chanpos, then coilpos, then elecpos,
%   whichever is present.
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

% ── Find the sensor struct for this array ────────────────────────────────
candidates = { ...
    [array_side '_coils_3axis'], ...
    [array_side '_coils_2axis'], ...
    [array_side '_sensors'] };

sens = [];
for c = 1:numel(candidates)
    if isfield(geom, candidates{c})
        sens = geom.(candidates{c});
        break
    end
end

if isempty(sens)
    error(['sim_sensor_positions: no sensor struct for the "%s" array in\n  %s\n' ...
           'Looked for: %s'], array_side, geom_file, strjoin(candidates, ', '));
end

% ── Find the position field ──────────────────────────────────────────────
% chanpos is the FieldTrip channel position; coilpos / elecpos are the raw
% MEG / EEG equivalents. Prefer chanpos where present so MEG positions match
% what the leadfield was computed on.
if isfield(sens, 'chanpos')
    pos = sens.chanpos;
elseif isfield(sens, 'coilpos')
    pos = sens.coilpos;
elseif isfield(sens, 'elecpos')
    pos = sens.elecpos;
else
    error(['sim_sensor_positions: sensor struct for the "%s" array has no ' ...
           'chanpos, coilpos, or elecpos field in\n  %s'], array_side, geom_file);
end

% ── Split into sensor axes ───────────────────────────────────────────────
% n_axes comes from the leadfield, so a mismatch here means the geometry and
% the leadfield disagree about the array — which would misalign every topoplot.
n_total = size(pos, 1);
if mod(n_total, n_axes) ~= 0
    error(['sim_sensor_positions: %d sensor positions in the "%s" array do not ' ...
           'divide into %d sensor axes.\nGeometry and leadfield disagree about ' ...
           'this array:\n  %s'], n_total, array_side, n_axes, geom_file);
end

n_per_axis  = n_total / n_axes;
pos_by_axis = cell(1, n_axes);
for ax = 1:n_axes
    idx1 = (ax-1)*n_per_axis + 1;
    idx2 =  ax   *n_per_axis;
    pos_by_axis{ax} = pos(idx1:idx2, :);
end
end

function lf = sim_load_leadfield(file, varname, scale, is_meg)
% sim_load_leadfield - Load one leadfield .mat and split it by sensor axis
%
% Unlike organise_leadfield (msg_fwd), this returns a STANDALONE struct for a
% single array rather than writing into a shared container. That matters here:
% organise_leadfield overwrites its target key on each call, so loading front
% and back into one key would silently discard the first. The sim scripts keep
% the two arrays separate on purpose.
%
% USAGE:
%   lf = sim_load_leadfield(file, varname, scale, is_meg)
%
% INPUT:
%   file    - full path to the leadfield .mat
%   varname - variable name inside it ('leadfield_cord', 'leadfield_bs', ...)
%   scale   - unit scale factor applied to every leadfield entry
%             (1e15 for T/nAm -> fT/nAm; 1e6 for V/nAm -> uV/nAm; 1 if already
%              in the target unit)
%   is_meg  - true = magnetic, false = electric
%
% OUTPUT:
%   lf.VD / lf.RC / lf.LR   - {n_axes x n_sources} cell, each [n_sens x 1]
%   lf.L                    - {n_axes x n_sources} cell, each [n_sens x 3]
%                             columns ordered [LR, RC, VD] (FieldTrip xyz)
%   lf.n_sources            - number of source positions
%   lf.n_sensor_axes        - 3 for triaxial MSG, 2 for the ESG electrode split
%   lf.n_sensors_per_axis   - sensors per axis
%   lf.is_meg               - passthrough
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if ~isfile(file)
    error('sim_load_leadfield: file not found:\n  %s', file);
end

tmp = load(file, varname);
if ~isfield(tmp, varname)
    error('sim_load_leadfield: variable "%s" not found in:\n  %s', varname, file);
end
lf_struct = tmp.(varname);

n_sources        = numel(lf_struct.leadfield);
n_channels_total = size(lf_struct.leadfield{1}, 1);

% Detect number of sensor axes the same way organise_leadfield does, so the
% sim module and the perturbation module agree on channel ordering.
if mod(n_channels_total, 3) == 0
    n_sensor_axes = 3;
elseif mod(n_channels_total, 2) == 0
    n_sensor_axes = 2;
else
    n_sensor_axes = 1;
end
n_sensors_per_axis = n_channels_total / n_sensor_axes;

lf                    = struct();
lf.VD                 = cell(n_sensor_axes, n_sources);
lf.RC                 = cell(n_sensor_axes, n_sources);
lf.LR                 = cell(n_sensor_axes, n_sources);
lf.L                  = cell(n_sensor_axes, n_sources);
lf.n_sources          = n_sources;
lf.n_sensor_axes      = n_sensor_axes;
lf.n_sensors_per_axis = n_sensors_per_axis;
lf.is_meg             = is_meg;

for s = 1:n_sources
    lf_matrix = lf_struct.leadfield{s} * scale;
    for ax = 1:n_sensor_axes
        idx1    = (ax-1)*n_sensors_per_axis + 1;
        idx2    =  ax   *n_sensors_per_axis;
        lf_axis = lf_matrix(idx1:idx2, :);

        lf.L{ax, s}  = lf_axis;
        lf.LR{ax, s} = lf_axis(:, 1);
        lf.RC{ax, s} = lf_axis(:, 2);
        lf.VD{ax, s} = lf_axis(:, 3);
    end
end
end

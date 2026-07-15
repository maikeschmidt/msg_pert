function f = sim_geom_file(model, geo)
% sim_geom_file - Path to the geometry .mat holding sensor positions
%
% Sensor positions come from the geometry file that the leadfield was built on.
% For a 'cond' variant the sensors are those of the geometry the conductivity
% perturbation was applied to (geo.short), which is the unperturbed original —
% conductivity changes do not move sensors.
%
% USAGE:
%   f = sim_geom_file(model, geo)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

f = fullfile(model.geoms_path, ['geometries_' geo.short '.mat']);
end

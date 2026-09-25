function s = pt_method_label(method)
% pt_method_label - Display name of a forward model method key
%
% USAGE:
%   s = pt_method_label('bslaw')   % -> 'Biot-Savart'
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

switch lower(method)
    case 'bslaw',  s = 'Biot-Savart';
    case 'bem',    s = 'BEM';
    case 'fem',    s = 'FEM';
    case 'sphere', s = 'Sphere';
    otherwise,     s = method;
end
end

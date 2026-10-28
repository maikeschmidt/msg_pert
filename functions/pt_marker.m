function m = pt_marker(k)
% pt_marker - The k-th marker shape of msg_pert's fixed sequence
%
% One sequence for every figure, so the same position in a legend always
% has the same shape. pt_distinct_markers assigns these automatically; call
% this directly only where a marker must match a series drawn separately.
%
% USAGE:
%   m = pt_marker(k)     % one marker, k = 1, 2, ... (cycles after the last)
%   M = pt_marker()      % the whole sequence
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

M = {'o', 's', '^', 'd', 'v', 'p', '>', 'h', '<', 'x', '+', '*'};
if nargin == 0
    m = M;
else
    m = M{mod(k - 1, numel(M)) + 1};
end
end

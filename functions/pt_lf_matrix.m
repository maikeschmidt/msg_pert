function L = pt_lf_matrix(E, ori, ax, n_trunc, src_idx)
% pt_lf_matrix - Leadfield matrix for one model, orientation and axis set
%
% Builds the comparison vectors used throughout msg_pert from one entry of
% the organised leadfields struct. Every metric in the pipeline is computed
% on matrices from this function, so all stages compare identical data.
%
% USAGE:
%   L = pt_lf_matrix(E, ori, ax, n_trunc)
%   L = pt_lf_matrix(E, ori, ax, n_trunc, src_idx)
%
% INPUT:
%   E        - one model from leadfields_organised.mat, e.g.
%              leadfields.bem_original_source_original
%   ori      - 'VD' | 'RC' | 'LR' for one dipole orientation, or 'ALL' to
%              stack [LR; RC; VD] (msg_fwd's concatenated convention)
%   ax       - sensor axis index, or 0 for the whole array (every sensor
%              axis stacked into one vector)
%   n_trunc  - sensors per axis to keep (the common count across models)
%   src_idx  - source indices to return (default: all)
%
% OUTPUT:
%   L        - [n_rows x numel(src_idx)], column s = comparison vector for
%              source src_idx(s). n_rows = n_ori * n_axes_used * n_trunc.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 5 || isempty(src_idx), src_idx = 1:E.n_sources; end

if ax == 0
    axes_used = 1:E.n_sensor_axes;
else
    if ax > E.n_sensor_axes
        error('pt_lf_matrix:badAxis', 'Axis %d requested, model has %d.', ...
            ax, E.n_sensor_axes);
    end
    axes_used = ax;
end

if strcmpi(ori, 'ALL')
    oris = {'LR', 'RC', 'VD'};
else
    oris = {upper(ori)};
end

n_blk = numel(oris) * numel(axes_used);
L = zeros(n_blk * n_trunc, numel(src_idx));

for si = 1:numel(src_idx)
    s   = src_idx(si);
    row = 0;
    for o = 1:numel(oris)
        for a = axes_used
            v = E.(oris{o}){a, s};
            L(row + (1:n_trunc), si) = v(1:n_trunc);
            row = row + n_trunc;
        end
    end
end
end

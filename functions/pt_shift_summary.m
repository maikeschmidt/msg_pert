function [v, X] = pt_shift_summary(S, metric, c, o)
% pt_shift_summary - One value per perturbation realisation, from a metric store
%
% Summarises each perturbation realisation by the median of a metric over
% cord positions. These per-realisation values are the unit of observation
% for every statistical test in msg_pert.
%
% USAGE:
%   [v, X] = pt_shift_summary(S, metric, c, o)
%
% INPUT:
%   S       - metric store PM.pert.<type>.M.<method> from pt_compute_metrics
%             (fields re, rsq, rdm, lnmag, each [n x n_src x n_cfg x n_ori])
%   metric  - 're' | 'rsq' | 'rdm' | 'gain' | 'absgain'
%   c       - axis-set index (position in PM.axcfg_ids)
%   o       - orientation index (position in PM.ori_labels)
%
% OUTPUT:
%   v       - [n x 1] median over cord of the metric, per realisation
%   X       - [n x n_src] the metric at every cord position
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

sub = struct('re', S.re(:, :, c, o), 'rsq', S.rsq(:, :, c, o), ...
             'rdm', S.rdm(:, :, c, o), 'lnmag', S.lnmag(:, :, c, o));
[~, X] = pt_metric_info(sub, metric);
v = median(X, 2, 'omitnan');
v(all(isnan(X), 2)) = NaN;
end

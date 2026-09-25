function M = pt_metrics_block(LA, LB, opts)
% pt_metrics_block - Vectorised four-metric comparison of matched leadfields
%
% Column-by-column RE, r2, RDM and lnMAG between a reference and a comparison
% leadfield matrix. Numerically identical to msg_fwd's lf_metrics_series under
% the toolbox default settings, but computed for every column at once, which
% matters in the noise simulation where the metrics are evaluated millions of
% times. If metric_defaults() selects a non-default mode the call is passed
% straight to lf_metrics_series, so definitions can never drift apart.
%
% USAGE:
%   M = pt_metrics_block(LA, LB)
%   M = pt_metrics_block(LA, LB, opts)
%
% INPUT:
%   LA    - [n x n_src] REFERENCE leadfields (the RE denominator)
%   LB    - [n x n_src] COMPARISON leadfields
%   opts  - metric options (default: metric_defaults())
%
% OUTPUT:
%   M     - struct of [1 x n_src] row vectors:
%             .re      relative error (%), ||LB-LA|| / ||LA|| * 100
%             .rsq     squared Pearson correlation
%             .rdm     ||LB/||LB|| - LA/||LA|| ||, topography only, [0, 2]
%             .lnmag   log(||LB|| / ||LA||), magnitude only
%             .re_ref  RE, reference-normalised convention
%             .re_sym  RE, symmetric L1 convention
%           Columns with a vanishing norm return NaN for every metric, and
%           constant columns return NaN for rsq, exactly as lf_metrics does.
%
% SEE ALSO:
%   lf_metrics, lf_metrics_series, metric_defaults
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 3 || isempty(opts), opts = metric_defaults(); end

re_mode  = 'reference';
rsq_mode = 'pearson';
if isfield(opts, 're_mode'),  re_mode  = lower(opts.re_mode);  end
if isfield(opts, 'rsq_mode'), rsq_mode = lower(opts.rsq_mode); end

if ~any(strcmp(re_mode, {'reference', 'eq13', 'symmetric'})) || ...
        ~strcmp(rsq_mode, 'pearson')
    M = lf_metrics_series(LA, LB, opts);
    return
end

if ~isequal(size(LA), size(LB))
    error('pt_metrics_block:sizeMismatch', ...
        'LA %s and LB %s must be the same size. Truncate first.', ...
        mat2str(size(LA)), mat2str(size(LB)));
end

nA = sqrt(sum(LA.^2, 1));
nB = sqrt(sum(LB.^2, 1));
D  = LB - LA;

M.re_ref = sqrt(sum(D.^2, 1)) ./ nA * 100;
M.re_sym = sum(abs(D), 1) ./ (sum(abs(LA), 1) + sum(abs(LB), 1)) * 100;
if strcmp(re_mode, 'symmetric')
    M.re = M.re_sym;
else
    M.re = M.re_ref;
end

M.rdm   = sqrt(sum((LB ./ nB - LA ./ nA).^2, 1));
M.lnmag = log(nB ./ nA);

Ac = LA - mean(LA, 1);
Bc = LB - mean(LB, 1);
r  = sum(Ac .* Bc, 1) ./ sqrt(sum(Ac.^2, 1) .* sum(Bc.^2, 1));
M.rsq = r.^2;
M.rsq(std(LA, 0, 1) < eps | std(LB, 0, 1) < eps) = NaN;

bad = nA < 1e-30 | nB < 1e-30;
for f = {'re', 'rsq', 'rdm', 'lnmag', 're_ref', 're_sym'}
    M.(f{1})(bad) = NaN;
end
end

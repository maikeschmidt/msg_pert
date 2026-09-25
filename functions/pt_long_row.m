function r = pt_long_row(meta, x, n_boot)
% pt_long_row - One row of the common long-format results table
%
% Every stage appends rows in this one layout to its *_long.csv, and
% pt_summary_table reads them all back, so any factor from any stage can be
% placed on the same scale as any other.
%
% USAGE:
%   r = pt_long_row(meta, x)
%   r = pt_long_row(meta, x, n_boot)
%
% INPUT:
%   meta   - struct: .stage .modality .system .method .factor .level
%            .orientation .axis_set .unit   (text; unit names what one
%            observation is, e.g. 'source' or 'realisation')
%   x      - struct of observation vectors: .re .rsq .rdm .gain .absgain
%   n_boot - bootstrap draws for the RE CI (default 10000)
%
% OUTPUT:
%   r      - 1-row table; pt_long_row([]) returns an empty table with the
%            same variables, for initialisation
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

vars = {'stage', 'modality', 'system', 'method', 'factor', 'level', ...
        'orientation', 'axis_set', 'unit', 'n', 're_median', 're_ci_lo', ...
        're_ci_hi', 're_q1', 're_q3', 'rsq_median', 'rdm_median', ...
        'gain_median', 'absgain_median'};

if isempty(meta)
    r = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
    return
end
if nargin < 3 || isempty(n_boot), n_boot = 10000; end

D  = pt_describe(x.re, n_boot);
md = @(v) median(v(~isnan(v)));
r = cell2table({meta.stage, meta.modality, meta.system, meta.method, ...
    meta.factor, meta.level, meta.orientation, meta.axis_set, meta.unit, ...
    D.n, D.median, D.ci_lo, D.ci_hi, D.q1, D.q3, md_or_nan(md, x.rsq), ...
    md_or_nan(md, x.rdm), md_or_nan(md, x.gain), md_or_nan(md, x.absgain)}, ...
    'VariableNames', vars);
end


function v = md_or_nan(md, x)
    if all(isnan(x(:))), v = NaN; else, v = md(x(:)); end
end

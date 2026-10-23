function r = pt_compare_row(R, family, subset, comparison)
% pt_compare_row - One row of a comparison table, from pt_compare_two
%
% Every comparison table in msg_pert has the same columns, so tables from
% different stages can be stacked and read the same way.
%
% USAGE:
%   vars = pt_compare_row()                          % the column names
%   r    = pt_compare_row(R, family, subset, comparison)
%
% INPUT:
%   R           - output of pt_compare_two
%   family      - comparison family, e.g. 'type_pairs' or 'between'
%   subset      - bundle name or 'pooled'
%   comparison  - label, 'A vs B'
%
% OUTPUT:
%   r           - 1 x 13 cell row matching the column names:
%                 family, subset, comparison, design, n_a, n_b, median_a,
%                 median_b, estimate, estimate_name, ci_lo, ci_hi, separated
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin == 0
    r = {'family', 'subset', 'comparison', 'design', 'n_a', 'n_b', ...
         'median_a', 'median_b', 'estimate', 'estimate_name', 'ci_lo', ...
         'ci_hi', 'separated'};
    return
end

if strcmp(R.design, 'paired'), nm = 'median of A-B'; else, nm = 'median A - median B'; end
r = {family, subset, comparison, R.design, R.n_a, R.n_b, R.median_a, ...
     R.median_b, R.median_diff, nm, R.ci_lo, R.ci_hi, R.separated};
end

function T = pt_fdr_by_group(T, group_vars, q)
% pt_fdr_by_group - Benjamini-Hochberg adjustment within groups of a table
%
% Adds p_fdr and sig columns. Rows are corrected together when they share
% the values of every variable in group_vars — typically the test family and
% the metric, so that, for example, every MSG-vs-ESG test on RE is one family
% across perturbation types, bundles, orientations and sensor axes.
%
% USAGE:
%   T = pt_fdr_by_group(T, group_vars)
%   T = pt_fdr_by_group(T, group_vars, q)
%
% INPUT:
%   T           - table with a numeric column p
%   group_vars  - cellstr of variable names defining a correction family
%   q           - FDR level for the sig column (default 0.05)
%
% OUTPUT:
%   T           - input table with p_fdr (adjusted) and sig (p_fdr < q)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 3 || isempty(q), q = 0.05; end

T.p_fdr = nan(height(T), 1);
T.sig   = false(height(T), 1);
if height(T) == 0, return; end

key = strings(height(T), 1);
for g = 1:numel(group_vars)
    key = key + "|" + string(T.(group_vars{g}));
end

[~, ~, gid] = unique(key);
for k = 1:max(gid)
    sel = gid == k;
    [padj, sig] = st_bh_fdr(T.p(sel), q);
    T.p_fdr(sel) = padj;
    T.sig(sel)   = sig;
end
end

function [info, V] = pt_metric_info(S, name)
% pt_metric_info - The metrics msg_pert reports, and how to read them from a store
%
% USAGE:
%   info = pt_metric_info()            % list of reported metrics
%   [~, V] = pt_metric_info(S, name)   % values of one metric from a store
%
% INPUT:
%   S     - struct with fields re, rsq, rdm, lnmag (arrays of equal size)
%   name  - 're' | 'rsq' | 'rdm' | 'gain' | 'absgain' | 'lnmag'
%
% OUTPUT:
%   info  - struct array, one entry per reported metric:
%             .name    field name
%             .label   axis / table label (TeX)
%             .txt     plain-text label
%             .tested  true if the metric enters the statistical tests
%             .worse   +1 if larger values mean a larger model error,
%                      -1 if smaller values do (r2)
%   V     - the requested metric, same size as S.re
%
% THE METRICS
%   re       relative error (%): magnitude and shape together
%   rsq      squared Pearson correlation: shape only
%   rdm      relative difference measure: shape only, on unit-norm fields
%   gain     (exp(lnMAG) - 1) * 100: signed amplitude change (%)
%   absgain  |gain|: size of the amplitude change whatever its sign
%
%   Signed gain is descriptive only. Across random shifts the sign varies
%   (moving the cord towards the array raises the field, away lowers it),
%   so a signed mean can sit near zero while every shift changes amplitude
%   substantially. Tests on amplitude therefore use absgain.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

info = struct( ...
    'name',   {'re',     'rsq',   'rdm',  'gain',     'absgain'}, ...
    'label',  {'RE (%)', 'r^2',   'RDM',  'Gain (%)', '|Gain| (%)'}, ...
    'txt',    {'RE (%)', 'r2',    'RDM',  'gain (%)', '|gain| (%)'}, ...
    'tested', {true,     true,    true,   false,      true}, ...
    'worse',  {1,        -1,      1,      1,          1});

if nargin < 2
    V = [];
    return
end

switch name
    case {'re', 'rsq', 'rdm', 'lnmag'}
        V = S.(name);
    case 'gain'
        V = (exp(S.lnmag) - 1) * 100;
    case 'absgain'
        V = abs(exp(S.lnmag) - 1) * 100;
    otherwise
        error('pt_metric_info:unknown', 'Unknown metric "%s".', name);
end
end

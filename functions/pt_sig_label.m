function s = pt_sig_label(p)
% pt_sig_label - Star label for an (FDR-adjusted) p-value
%
% USAGE:
%   s = pt_sig_label(p)
%
% OUTPUT:
%   '***' p < 0.001, '**' p < 0.01, '*' p < 0.05, 'n.s.' otherwise,
%   '' for NaN (no test was possible).
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if isnan(p),      s = '';
elseif p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
else,             s = 'n.s.';
end
end

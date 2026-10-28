function s = pt_sep_label(separated)
% pt_sep_label - Figure mark for a comparison whose interval excludes zero
%
% Figures in msg_pert mark comparisons by their interval, not by a p-value:
% a diamond means the bootstrap interval of the difference (or of the
% correlation) excludes zero. Every figure that uses it says so in its title.
%
% The mark is the Unicode white diamond (U+25C7), not a TeX command, so it
% renders whatever text interpreter MATLAB is set to; draw it with
% 'Interpreter', 'none'.
%
% USAGE:
%   s = pt_sep_label(separated)
%
% OUTPUT:
%   char(9671), the white diamond, when separated is true, '' otherwise (including NaN or
%   an empty input, where no interval could be formed).
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if ~isempty(separated) && islogical(separated) && separated(1)
    s = char(9671);
else
    s = '';
end
end

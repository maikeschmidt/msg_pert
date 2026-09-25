function [s, why] = pt_unit_scale(lf_struct, is_meg, method)
% pt_unit_scale - Unit scale factor for a saved leadfield, from its magnitude
%
% Leadfields reach msg_pert from several writers that do not agree on units:
% raw FieldTrip BEM output per nA*m or per A*m (a factor of 1e9 apart,
% depending on whether the cfg.dipoleunit patch is active), Biot-Savart
% already in fT/nAm, and conductivity-perturbation files that carry an extra
% x1e15 baked in at save time. Applying one fixed factor to all of them puts
% a perturbed leadfield orders of magnitude away from its reference, which
% leaves r2 untouched (it is scale invariant) but makes RE and lnMAG
% meaningless.
%
% The candidate factors for each modality are at least six orders of
% magnitude apart, and a physical leadfield sits in a band narrower than
% that, so the factor that lands the median |L| inside the band is chosen.
%
% USAGE:
%   [s, why] = pt_unit_scale(lf_struct, is_meg)
%   [s, why] = pt_unit_scale(lf_struct, is_meg, method)
%
% INPUT:
%   lf_struct - loaded leadfield struct (with a .leadfield cell array)
%   is_meg    - true: target unit fT/nAm; false: target unit uV/nAm
%   method    - optional, used only for the fallback when no candidate fits
%
% OUTPUT:
%   s         - multiplicative factor to the target unit
%   why       - short explanation, for the load log
%
% CANDIDATES:
%   MSG (fT/nAm), plausible median |L| in [1e-3, 1e2]:
%     1e15  raw T per nA*m          1e6   raw T per A*m
%     1     already fT/nAm          1e-9  fT per A*m (x1e15 baked, per A*m)
%   ESG (uV/nAm), plausible median |L| in [1e-4, 30]:
%     1e6   raw V per nA*m          1e-3  raw V per A*m
%     1e-9  V per nA*m, x1e15 baked 1e-18 V per A*m, x1e15 baked
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 3, method = ''; end

if is_meg
    cand  = [1e15, 1e6, 1, 1e-9];
    names = {'raw T per nA*m', 'raw T per A*m', 'already fT/nAm', ...
             'fT per A*m (baked x1e15)'};
    band  = [1e-3, 1e2];
    unit  = 'fT/nAm';
    if startsWith(lower(method), 'bem'), fallback = 1e15; else, fallback = 1; end
else
    cand  = [1e6, 1e-3, 1e-9, 1e-18];
    names = {'raw V per nA*m', 'raw V per A*m', ...
             'V per nA*m (baked x1e15)', 'V per A*m (baked x1e15)'};
    band  = [1e-4, 30];
    unit  = 'uV/nAm';
    fallback = 1e6;
end

m = NaN;
if isstruct(lf_struct) && isfield(lf_struct, 'leadfield')
    L = lf_struct.leadfield;
    if iscell(L)
        L = L(~cellfun(@isempty, L));
        if ~isempty(L), L = cell2mat(L(:)); else, L = []; end
    end
    if isnumeric(L) && ~isempty(L)
        m = median(abs(L(:)), 'omitnan');
    end
end

if ~isfinite(m) || m <= 0
    s   = fallback;
    why = sprintf('magnitude unreadable, fallback x%g', s);
    warning('pt_unit_scale:unreadable', ...
        'Could not measure leadfield magnitude; using fallback x%g.', s);
    return
end

ok = find(m * cand >= band(1) & m * cand <= band(2));

if isscalar(ok)
    s   = cand(ok);
    why = sprintf('%s (median %.3g %s)', names{ok}, m * s, unit);
elseif isempty(ok)
    s   = fallback;
    why = sprintf('no candidate plausible, fallback x%g', s);
    warning('pt_unit_scale:noPlausibleScale', ...
        ['Median |L| = %.3g fits no expected unit convention for %s. ' ...
         'Using x%g. Check the file, or set unit_scale_mode = ''fixed'' in ' ...
         'config_pert.'], m, unit, s);
else
    % Cannot happen with the bands above; kept so a future band edit fails loudly
    s   = cand(ok(1));
    why = sprintf('ambiguous, chose %s', names{ok(1)});
    warning('pt_unit_scale:ambiguous', ...
        'Median |L| = %.3g matches more than one convention; chose x%g.', m, s);
end
end

function L = pt_dedup_electrodes(L, n_copies)
% pt_dedup_electrodes - Keep one copy of an electrode set stored n times
%
% msg_coreg builds ESG electrodes from the coil positions of a triaxial
% magnetometer grid, which hold every grid point three times; the electrode
% array keeps two of those copies. An electrode measures a potential, which
% has no orientation, so the copies are the same electrodes measuring the
% same thing: their lead-field rows are identical. Treating them as two
% sensor "axes" counts every electrode twice and, through noise, gives each
% copy independent noise it would never have. This keeps the first copy.
%
% The lead field of the first copy is exactly the lead field of the single
% set, including its common-average reference: averaging over every
% electrode n times gives the same mean as averaging over each once. So no
% lead field has to be recomputed.
%
% The copies are checked before anything is dropped; if they differ, the
% array is not a duplicated set and the function stops rather than discard
% real channels.
%
% USAGE:
%   L = pt_dedup_electrodes(L, n_copies)
%
% INPUT:
%   L         - FieldTrip lead-field struct (.leadfield = {1 x n_src} of
%               [n_channels x 3]); optional .label
%   n_copies  - how many times the geometry stores each electrode
%               (config_pert: electrode_copies). 1 = return L unchanged.
%
% OUTPUT:
%   L         - the same struct with n_channels / n_copies rows per source
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 2 || isempty(n_copies) || n_copies <= 1, return; end

for s = 1:numel(L.leadfield)
    X = L.leadfield{s};
    if isempty(X), continue; end
    n_tot = size(X, 1);
    if mod(n_tot, n_copies) ~= 0
        error('pt_dedup_electrodes:size', ...
            '%d channels cannot be %d copies of one electrode set.', n_tot, n_copies);
    end
    n   = n_tot / n_copies;
    ref = X(1:n, :);
    tol = 1e-9 * max(abs(ref(:)));
    for c = 2:n_copies
        d = max(abs(X((c-1)*n + (1:n), :) - ref), [], 'all');
        if d > tol
            error('pt_dedup_electrodes:notCopies', ...
                ['Copy %d of the electrode set differs from copy 1 (max |diff| %.3g, ' ...
                 'source %d): these are not duplicated electrodes. Set ' ...
                 'electrode_copies = 1 for this array.'], c, d, s);
        end
    end
    L.leadfield{s} = ref;
end

if isfield(L, 'label') && iscell(L.label) && mod(numel(L.label), n_copies) == 0
    L.label = L.label(1:numel(L.label) / n_copies);
end
end

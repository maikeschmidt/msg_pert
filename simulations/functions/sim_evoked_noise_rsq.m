function [rsq_mean, rsq_sd, snr] = sim_evoked_noise_rsq(lfvec, w, sigma_vec, n_real, ref_lfvec)
% sim_evoked_noise_rsq - r^2 of a noisy evoked field vs a noise-free reference
%
% The exact closed-form used by the single-geometry sim, factored out so the
% geometry loop reuses identical maths. For each source the measured sensor data
% is the rank-1 evoked field plus Gaussian noise,
%
%     Y = g * w' + sigma * Z ,   Z ~ N(0,1),
%
% scored by Pearson r^2 against a noise-free reference field a*w' over all
% sensor-by-time elements. The reference a defaults to the signal g itself
% (r^2 vs the field's own noise-free version — pure noise robustness). Pass a
% DIFFERENT reference (e.g. the unperturbed/original leadfield) to score the
% combined effect of model error AND noise: r^2 then starts below 1 even at zero
% noise, by the amount the perturbed field already differs from the reference.
%
% Because both signal and reference are rank-1 the correlation collapses to a few
% scalar sums per source, so no [n_ch x n_time] matrix is formed. One noise draw
% is shared across all noise levels within a realisation (common random numbers).
%
% USAGE:
%   [rsq_mean, rsq_sd, snr] = sim_evoked_noise_rsq(lfvec, w, sigma_vec, n_real)
%   [...] = sim_evoked_noise_rsq(lfvec, w, sigma_vec, n_real, ref_lfvec)
%
% INPUT:
%   lfvec      [n_ch x n_src]  measured (perturbed) leadfield per source, one
%              orientation, already in the same physical unit as sigma_vec
%   w          [n_time x 1]    source waveform (nA*m)
%   sigma_vec  [1 x n_lev]     effective noise s.d. per level (bandwidth + trials)
%   n_real     scalar          noise realisations
%   ref_lfvec  [n_ch x n_src]  OPTIONAL noise-free reference leadfield. Defaults
%              to lfvec (self). Must share channel ordering and source grid.
%
% OUTPUT:
%   rsq_mean [n_src x n_lev]   mean r^2 over realisations
%   rsq_sd   [n_src x n_lev]   s.d. of r^2 across realisations
%   snr      [n_src x n_lev]   amplitude SNR, rms(signal)/sigma
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

if nargin < 5 || isempty(ref_lfvec)
    ref_lfvec = lfvec;                 % self reference (noise robustness only)
end
if ~isequal(size(ref_lfvec), size(lfvec))
    error('sim_evoked_noise_rsq: reference size %s does not match signal size %s.', ...
        mat2str(size(ref_lfvec)), mat2str(size(lfvec)));
end

[n_ch, n_src] = size(lfvec);
n_time = numel(w);
n_el   = n_ch * n_time;
n_lev  = numel(sigma_vec);

Sw  = sum(w);
Sww = w' * w;

% Signal (measured field) sums drive the noisy data moments
Sg  = sum(lfvec, 1);       % [1 x n_src]
Sgg = sum(lfvec.^2, 1);

% Reference (noise-free) sums drive the reference moments
Sa  = sum(ref_lfvec, 1);
Saa = sum(ref_lfvec.^2, 1);

rsq_all = nan(n_src, n_lev, n_real);
snr     = nan(n_src, n_lev);

for s = 1:n_src
    g   = lfvec(:, s);
    a   = ref_lfvec(:, s);
    Sx  = Sa(s)  * Sw;         % reference sum
    Sxx = Saa(s) * Sww;        % reference energy
    if Sxx <= 0
        continue               % degenerate (all-zero) reference field
    end
    Sag  = a' * g;             % reference-signal cross term
    denX = n_el*Sxx - Sx^2;

    for it = 1:n_real
        Z   = randn(n_ch, n_time);
        u   = Z * w;
        SZ  = sum(Z(:));
        SZZ = sum(Z(:).^2);
        gu  = g' * u;
        au  = a' * u;

        for L = 1:n_lev
            sg  = sigma_vec(L);
            Sy  = Sg(s)*Sw + sg*SZ;                 % noisy data sum
            Sxy = Sag*Sww  + sg*au;                 % reference . noisy
            Syy = Sgg(s)*Sww + 2*sg*gu + sg^2*SZZ;  % noisy energy

            denY = n_el*Syy - Sy^2;
            den  = denX * denY;
            if den <= 0
                continue
            end
            rr = (n_el*Sxy - Sx*Sy) / sqrt(den);
            rsq_all(s, L, it) = rr^2;

            if it == 1
                snr(s, L) = sqrt(Sgg(s)*Sww / n_el) / sg;   % rms(signal)/sigma
            end
        end
    end
end

rsq_mean = mean(rsq_all, 3, 'omitnan');
rsq_sd   = std(rsq_all, 0, 3, 'omitnan');
end

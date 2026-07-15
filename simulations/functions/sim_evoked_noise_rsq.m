function [rsq_mean, rsq_sd, snr] = sim_evoked_noise_rsq(lfvec, w, sigma_vec, n_real)
% sim_evoked_noise_rsq - r^2 of a noisy evoked field vs its noise-free self
%
% The exact closed-form used by the working single-geometry sim, factored out so
% the geometry loop reuses identical maths. For each source the sensor data is
% the rank-1 evoked field plus Gaussian noise,
%
%     Y = g * w' + sigma * Z ,   Z ~ N(0,1),
%
% scored by Pearson r^2 against the noise-free field g*w' over all sensor-by-time
% elements. Because the signal is rank-1 the correlation collapses to a few
% scalar sums per source, so no [n_ch x n_time] matrix is formed. One noise draw
% is shared across all noise levels within a realisation (common random numbers).
%
% USAGE:
%   [rsq_mean, rsq_sd, snr] = sim_evoked_noise_rsq(lfvec, w, sigma_vec, n_real)
%
% INPUT:
%   lfvec      [n_ch x n_src]  full-channel leadfield per source (one orientation;
%              already in the same physical unit as sigma_vec)
%   w          [n_time x 1]    source waveform (nA*m)
%   sigma_vec  [1 x n_lev]     effective noise s.d. per level (bandwidth + trials)
%   n_real     scalar          noise realisations
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

[n_ch, n_src] = size(lfvec);
n_time = numel(w);
n_el   = n_ch * n_time;
n_lev  = numel(sigma_vec);

Sw  = sum(w);
Sww = w' * w;

Sg  = sum(lfvec, 1);      % [1 x n_src]
Sgg = sum(lfvec.^2, 1);   % [1 x n_src]

rsq_all = nan(n_src, n_lev, n_real);
snr     = nan(n_src, n_lev);

for s = 1:n_src
    g   = lfvec(:, s);
    Sx  = Sg(s)  * Sw;
    Sxx = Sgg(s) * Sww;
    if Sxx <= 0
        continue        % degenerate (all-zero) field
    end
    denX = n_el*Sxx - Sx^2;

    for it = 1:n_real
        Z   = randn(n_ch, n_time);
        u   = Z * w;
        SZ  = sum(Z(:));
        SZZ = sum(Z(:).^2);
        gu  = g' * u;

        for L = 1:n_lev
            sg  = sigma_vec(L);
            Sy  = Sx  + sg*SZ;
            Sxy = Sxx + sg*gu;
            Syy = Sxx + 2*sg*gu + sg^2*SZZ;

            denY = n_el*Syy - Sy^2;
            den  = denX * denY;
            if den <= 0
                continue
            end
            rr = (n_el*Sxy - Sx*Sy) / sqrt(den);
            rsq_all(s, L, it) = rr^2;

            if it == 1
                snr(s, L) = sqrt(Sxx / n_el) / sg;   % rms(signal)/sigma
            end
        end
    end
end

rsq_mean = mean(rsq_all, 3, 'omitnan');
rsq_sd   = std(rsq_all, 0, 3, 'omitnan');
end

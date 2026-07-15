function rsq = sim_noise_rsq_refs(g, refs, w, sigma_mat, n_real)
% sim_noise_rsq_refs - Monte-Carlo r^2 of a noisy rank-1 signal vs several
%                      noise-free references, in closed form
%
% The simulated sensor data is the perturbed forward field driven by the source
% waveform, plus Gaussian sensor noise:
%
%     Y = g * w' + sigma * Z          Z ~ N(0,1) elementwise
%
% and we score it against one or more noise-free reference fields
%
%     X_r = a_r * w'
%
% by the Pearson r^2 over all sensor-by-time elements. Two references are used
% by sim_perturbation_noise:
%   * a = g            -> "self": how robust this perturbed field is to noise.
%   * a = h (perfect)  -> "combined": model error (perturbation) AND noise
%                         together, scored against the unperturbed field.
%
% CLOSED FORM
%   Both X and Y are rank-1 (an outer product) plus, for Y, noise. The Pearson
%   correlation over the n = n_ch*n_time elements collapses to a few scalar
%   sums, so the whole thing is O(n_ch) per reference after one O(n_ch*n_time)
%   noise contraction — no [n_ch x n_time] matrix is ever formed for scoring.
%   With
%       Sw=sum(w)  Sww=w'w                                   (waveform)
%       Sg=sum(g)  Sgg=g'g                                   (signal field)
%       Sa=sum(a)  Saa=a'a  Sag=a'g                          (reference field)
%       u = Z*w    SZ=sum(Z(:))  SZZ=sum(Z(:).^2)            (noise, per draw)
%   the correlation is
%       Sx=Sa*Sw           Sxx=Saa*Sww
%       Sy=Sg*Sw+sigma*SZ  Syy=Sgg*Sww+2*sigma*(g'u)+sigma^2*SZZ
%       Sxy=Sag*Sww+sigma*(a'u)
%       r = (n*Sxy - Sx*Sy) / sqrt((n*Sxx - Sx^2)*(n*Syy - Sy^2))
%   which is exact, not an approximation.
%
% COMMON RANDOM NUMBERS
%   One noise draw Z per realisation is shared across every reference, every
%   system, and every noise level (each just rescales it by sigma). Systems that
%   measure the same field therefore see the same noise pattern at different
%   amplitudes, so differences between their curves are due to the noise floor,
%   not to luck in the draw.
%
% USAGE:
%   rsq = sim_noise_rsq_refs(g, refs, w, sigma_mat, n_real)
%
% INPUT:
%   g          [n_ch x 1]  perturbed signal leadfield (the field being measured)
%   refs       {1 x n_ref} cell of [n_ch x 1] reference leadfields to score
%              against (e.g. {g, h_perfect})
%   w          [n_time x 1] source waveform (nA*m)
%   sigma_mat  [n_sys x n_lev] effective noise s.d. per system per level,
%              already including bandwidth and trial averaging
%   n_real     number of noise realisations
%
% OUTPUT:
%   rsq        [n_ref x n_sys x n_lev x n_real]  r^2 per reference/system/level/
%              realisation. NaN where the correlation is undefined (a degenerate
%              all-zero field or a zero-variance case).
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

    n_ref = numel(refs);
    [n_sys, n_lev] = size(sigma_mat);
    n_ch   = numel(g);
    n_time = numel(w);
    n_el   = n_ch * n_time;

    rsq = nan(n_ref, n_sys, n_lev, n_real);

    Sw  = sum(w);
    Sww = w' * w;

    Sg  = sum(g);
    Sgg = g' * g;

    % Per-reference constants
    Sa  = zeros(1, n_ref);
    Saa = zeros(1, n_ref);
    Sag = zeros(1, n_ref);
    Sx  = zeros(1, n_ref);
    Sxx = zeros(1, n_ref);
    for r = 1:n_ref
        a       = refs{r};
        Sa(r)   = sum(a);
        Saa(r)  = a' * a;
        Sag(r)  = a' * g;
        Sx(r)   = Sa(r)  * Sw;
        Sxx(r)  = Saa(r) * Sww;
    end

    for it = 1:n_real
        Z   = randn(n_ch, n_time);
        u   = Z * w;               % [n_ch x 1] — the only O(n_ch*n_time) step
        SZ  = sum(Z(:));
        SZZ = sum(Z(:).^2);
        gu  = g' * u;

        for r = 1:n_ref
            if Sxx(r) <= 0
                continue           % degenerate reference field
            end
            au   = refs{r}' * u;
            denX = n_el*Sxx(r) - Sx(r)^2;

            for s = 1:n_sys
                for L = 1:n_lev
                    sg = sigma_mat(s, L);

                    Sy  = Sg*Sw + sg*SZ;
                    Syy = Sgg*Sww + 2*sg*gu + sg^2*SZZ;
                    Sxy = Sag(r)*Sww + sg*au;

                    denY = n_el*Syy - Sy^2;
                    den  = denX * denY;
                    if den <= 0
                        continue
                    end
                    rr = (n_el*Sxy - Sx(r)*Sy) / sqrt(den);
                    rsq(r, s, L, it) = rr^2;
                end
            end
        end
    end
end

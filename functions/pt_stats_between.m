function T = pt_stats_between(va, vb, bund_a, bund_b, paired, opts, comparison)
% pt_stats_between - Compare two conditions at every bundle and pooled
%
% Used wherever the same perturbation is measured two ways — MSG against ESG,
% Biot-Savart against BEM, one noise system against another — and the
% question is whether the two respond differently. Each observation is one
% perturbation realisation summarised over the cord; each row is a median
% difference with a bootstrap interval (pt_compare_two).
%
% USAGE:
%   T = pt_stats_between(va, vb, bund_a, bund_b, paired, opts, comparison)
%
% INPUT:
%   va, vb          - [n x 1] per-realisation values of the two conditions
%   bund_a, bund_b  - [n x 1] bundle index of each value
%   paired          - true when va(i) and vb(i) are the same realisation
%                     (the same geometric shift) measured two ways
%   opts            - struct: .bundle_names, .n_boot, .level
%   comparison      - label, e.g. 'MSG vs ESG'
%
% OUTPUT:
%   T               - table with the columns of pt_compare_row, family
%                     'between', one row per bundle plus a pooled row
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

nb   = numel(opts.bundle_names);
vars = pt_compare_row();
rows = cell(nb + 1, numel(vars));

for b = 0:nb
    if b == 0
        ia = true(size(va)); ib = true(size(vb)); sub = 'pooled';
    else
        ia = bund_a == b;    ib = bund_b == b;    sub = opts.bundle_names{b};
    end
    pr = paired && sum(ia) == sum(ib) && isequal(ia, ib);
    R  = pt_compare_two(va(ia), vb(ib), pr, opts.n_boot, opts.level);
    rows(b + 1, :) = pt_compare_row(R, 'between', sub, comparison);
end

% Bundles first, pooled last
rows = rows([2:end, 1], :);

T = cell2table(rows, 'VariableNames', vars);
end

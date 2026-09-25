% test_pt_core - Regression test for msg_pert's metric and statistics core
%
% Verifies, against known answers:
%   - pt_metrics_block reproduces msg_fwd's lf_metrics_series exactly,
%     including degenerate columns
%   - the noise decomposition E||g + n - a||^2 = ||g - a||^2 + E||n||^2 used
%     by the stage-5 error budget
%   - pt_perm_test2: exact enumeration, symmetric null, p never 0 by Monte Carlo
%   - pt_cliffs_delta, pt_spearman_perm and pt_quantile against hand values
%   - pt_compare_two chooses the paired / unpaired route correctly
%   - pt_unit_scale picks the right factor for every writer convention
%   - pt_fdr_by_group corrects within groups only
%
% USAGE:
%   Run from anywhere with msg_pert and msg_fwd as sibling folders:
%     test_pt_core
%
% Every line should report OK. Any *** FAIL *** means a metric or statistic
% has changed behaviour and earlier results will not reproduce.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'functions'), ...
    fullfile(fileparts(root), 'msg_fwd', 'functions'));

n_fail = 0;
check = @(name, got, want, tol) report(name, got, want, tol);

rng(1);

% 1. pt_metrics_block == lf_metrics_series ------------------------------
LA = randn(60, 40); LB = LA + 0.3*randn(60, 40); LB(:, 5) = 3*LA(:, 5);
LA(:, 7) = 0;                                     % degenerate column
M1 = pt_metrics_block(LA, LB);
M2 = lf_metrics_series(LA, LB, metric_defaults());
for f = {'re', 'rsq', 'rdm', 'lnmag', 're_sym'}
    a = M1.(f{1}); b = M2.(f{1});
    same = isequal(isnan(a), isnan(b)) && max(abs(a(~isnan(a)) - b(~isnan(b)))) < 1e-10;
    n_fail = n_fail + check(sprintf('block == series: %s', f{1}), double(same), 1, 0);
end
n_fail = n_fail + check('pure gain x3: RE (%)',  M1.re(5),  200, 1e-9);
n_fail = n_fail + check('pure gain x3: RDM',     M1.rdm(5), 0,   1e-12);
n_fail = n_fail + check('degenerate column: NaN', double(isnan(M1.re(7))), 1, 0);

% 2. Noise decomposition of squared RE ----------------------------------
a  = randn(300, 1); g = a + 0.2*randn(300, 1); sd = 0.15; nr = 20000;
N  = sd * randn(300, nr);
re2_emp  = mean(sum((g + N - a).^2, 1)) / sum(a.^2);
re2_pred = (sum((g - a).^2) + 300 * sd^2) / sum(a.^2);
n_fail = n_fail + check('E RE^2 = RE_pert^2 + RE_noise^2', re2_emp, re2_pred, 0.01*re2_pred);

% 3. pt_perm_test2 ------------------------------------------------------
x = [5 6 7 8]; y = [1 2 3 4];
p = pt_perm_test2(x, y, 1e5, 'right');           % 70 relabellings, exact
n_fail = n_fail + check('perm2 exact, complete separation', p, 1/70, 1e-12);
p = pt_perm_test2(x, y, 1e5, 'both');
n_fail = n_fail + check('perm2 exact two-sided', p, 2/70, 1e-12);
p = pt_perm_test2(randn(30,1) + 5, randn(30,1), 2000, 'both');
n_fail = n_fail + check('perm2 Monte Carlo p never 0', p, 1/2001, 1e-12);
p = pt_perm_test2(randn(8,1), randn(8,1), 12870, 'both');
n_fail = n_fail + check('perm2 null p in (0,1]', double(p > 0 && p <= 1), 1, 0);

% 4. Effect sizes and correlation ---------------------------------------
n_fail = n_fail + check('Cliff delta, complete separation', pt_cliffs_delta(x, y), 1, 0);
n_fail = n_fail + check('Cliff delta, identical groups',   pt_cliffs_delta(y, y), 0, 0);
[rho, p] = pt_spearman_perm(1:10, (1:10).^3, 5000);
n_fail = n_fail + check('Spearman monotone = 1', rho, 1, 1e-12);
n_fail = n_fail + check('Spearman monotone p small', double(p < 0.01), 1, 0);
[rho, ~] = pt_spearman_perm(1:10, 10:-1:1, 100);
n_fail = n_fail + check('Spearman reversed = -1', rho, -1, 1e-12);

% 5. Quantiles match MATLAB's prctile convention -------------------------
v = [1 2 3 4 10];
q = pt_quantile(v, [25 50 75], 2);
n_fail = n_fail + check('quantile 25', q(1), 1.75, 1e-12);
n_fail = n_fail + check('quantile 50', q(2), 3,    1e-12);
n_fail = n_fail + check('quantile 75', q(3), 5.5,  1e-12);

% 6. pt_compare_two routes ---------------------------------------------
R = pt_compare_two([1 2 3 4 5 6]', [0 1 2 3 4 5]', true);
n_fail = n_fail + check('paired: all diffs +1 -> rank-biserial 1', R.effect, 1, 1e-12);
n_fail = n_fail + check('paired: median diff', R.median_diff, 1, 1e-12);
R = pt_compare_two([5 6 7 8]', [1 2 3 4]', false);
n_fail = n_fail + check('unpaired: Cliff delta 1', R.effect, 1, 1e-12);

% 7. Unit scale detection ----------------------------------------------
base_T = 1e-15 * (0.5 + rand(50, 3));            % ~1 fT/nAm expressed in T per nA*m
cases = { 'MSG raw T/nAm',          base_T,              true,  1e15; ...
          'MSG raw T/(A*m)',        base_T * 1e9,        true,  1e6;  ...
          'MSG already fT/nAm',     base_T * 1e15,       true,  1;    ...
          'MSG baked, per A*m',     base_T * 1e9 * 1e15, true,  1e-9; ...
          'ESG raw V/nAm',          base_T * 1e8,        false, 1e6;  ...
          'ESG raw V/(A*m)',        base_T * 1e17,       false, 1e-3; ...
          'ESG baked, per nA*m',    base_T * 1e23,       false, 1e-9; ...
          'ESG baked, per A*m',     base_T * 1e32,       false, 1e-18 };
for k = 1:size(cases, 1)
    lf = struct('leadfield', {{cases{k, 2}}});
    s  = pt_unit_scale(lf, cases{k, 3}, 'bem');
    n_fail = n_fail + check(['unit scale: ' cases{k, 1}], log10(s), log10(cases{k, 4}), 1e-9);
end

% 8. FDR within groups -------------------------------------------------
T = table({'a'; 'a'; 'b'; 'b'}, [0.01; 0.04; 0.01; 0.04], 'VariableNames', {'g', 'p'});
T = pt_fdr_by_group(T, {'g'});
n_fail = n_fail + check('FDR per group (not pooled)', T.p_fdr(3), 0.02, 1e-12);

fprintf('\n%d failure(s)\n', n_fail);


function f = report(name, got, want, tol)
    ok = abs(got - want) <= tol;
    if ok, s = 'OK'; else, s = '*** FAIL ***'; end
    fprintf('%-42s got %12.6g  want %12.6g  %s\n', name, got, want, s);
    f = double(~ok);
end

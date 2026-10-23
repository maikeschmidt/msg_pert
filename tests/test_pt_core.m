% test_pt_core - Regression test for msg_pert's metric and statistics core
%
% Verifies, against known answers:
%   - pt_metrics_block reproduces msg_fwd's lf_metrics_series exactly,
%     including degenerate columns
%   - the noise decomposition E||g + n - a||^2 = ||g - a||^2 + E||n||^2 used
%     by the stage-5 error budget
%   - pt_rank_corr and pt_quantile against hand values; pt_quantile agrees
%     with msg_fwd's pctl
%   - pt_compare_two chooses the paired / unpaired route and its interval
%     separates clear differences and covers zero for none
%   - pt_describe's coverage percentile and its vectorised bootstrap agree
%     with pctl
%   - pt_unit_scale picks the right factor for every writer convention
%   - pt_compute_leadfields writes msg_fwd's file layout and skips existing
%     files (Biot-Savart route, which needs no FieldTrip)
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

% 3-4. Rank correlation with interval ----------------------------------
[rho, ci] = pt_rank_corr(1:10, (1:10).^3, 2000);
n_fail = n_fail + check('Spearman monotone = 1', rho, 1, 1e-12);
n_fail = n_fail + check('Spearman monotone CI excludes 0', double(ci(1) > 0), 1, 0);
[rho, ~] = pt_rank_corr(1:10, 10:-1:1, 100);
n_fail = n_fail + check('Spearman reversed = -1', rho, -1, 1e-12);
[~, ci] = pt_rank_corr(1:3, 1:3, 100);
n_fail = n_fail + check('Spearman n < 4 -> NaN CI', double(all(isnan(ci))), 1, 0);

% 5. Quantiles match MATLAB's prctile convention -------------------------
v = [1 2 3 4 10];
q = pt_quantile(v, [25 50 75], 2);
n_fail = n_fail + check('quantile 25', q(1), 1.75, 1e-12);
n_fail = n_fail + check('quantile 50', q(2), 3,    1e-12);
n_fail = n_fail + check('quantile 75', q(3), 5.5,  1e-12);
w = randn(37, 1);
n_fail = n_fail + check('pt_quantile == pctl (p90)', pt_quantile(w, 90, 1), pctl(w, 90), 1e-12);

% 6. pt_compare_two routes and intervals --------------------------------
R = pt_compare_two([1 2 3 4 5 6]', [0 1 2 3 4 5]', true, 2000);
n_fail = n_fail + check('paired: median diff', R.median_diff, 1, 1e-12);
n_fail = n_fail + check('paired: constant diff -> CI [1, 1]', R.ci_hi - R.ci_lo, 0, 1e-12);
n_fail = n_fail + check('paired: separated', double(R.separated), 1, 0);
R = pt_compare_two([5 6 7 8]', [1 2 3 4]', false, 2000);
n_fail = n_fail + check('unpaired: difference of medians', R.median_diff, 4, 1e-12);
n_fail = n_fail + check('unpaired: complete separation', double(R.separated), 1, 0);
R = pt_compare_two(randn(8, 1), randn(8, 1) + 20, false, 2000);
n_fail = n_fail + check('unpaired: CI below 0 when B >> A', double(R.ci_hi < 0), 1, 0);
z = randn(24, 1);
R = pt_compare_two(z, z, true, 2000);
n_fail = n_fail + check('paired: identical -> not separated', double(R.separated), 0, 0);
R = pt_compare_two([1 2]', [3 4]', false, 200);
n_fail = n_fail + check('n < 3 -> NaN CI, not separated', ...
    double(isnan(R.ci_lo) && ~R.separated), 1, 0);

% 6b. Coverage percentile ----------------------------------------------
v = (1:30)';
D = pt_describe(v, 2000);
n_fail = n_fail + check('describe: p95 == pctl', D.cover, pctl(v, 95), 1e-12);
n_fail = n_fail + check('describe: p95 CI brackets p95', ...
    double(D.cover_ci_lo <= D.cover && D.cover <= D.cover_ci_hi), 1, 0);

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

% 8. Lead fields written in msg_fwd's layout, and resumed --------------
gd   = tempname; mkdir(gd);
grad = struct('coilpos', [0 0 50; 10 0 50; 0 10 50], 'coilori', repmat([0 0 1], 3, 1), ...
              'tra', eye(3), 'label', {{'c1'; 'c2'; 'c3'}}, 'unit', 'mm');
g = struct('back_coils_3axis', grad, 'sources_cent', struct('pos', [0 0 0; 0 5 0]));
save(fullfile(gd, 'geometries_k1.mat'), '-struct', 'g');
Rl = pt_compute_leadfields('bslaw', gd, {'k1'}, fullfile(gd, 'bs'));
f  = fullfile(gd, 'bs', 'leadfield_geometries_k1_bslaw_back.mat');
n_fail = n_fail + check('BS file in msg_fwd layout', double(isfile(f) && Rl.written == 1), 1, 0);
if isfile(f)
    Lf = load(f, 'leadfield_bs');
    Lr = bs_leadfield(grad, g.sources_cent.pos);
    n_fail = n_fail + check('BS file == bs_leadfield', ...
        max(abs(Lf.leadfield_bs.leadfield{2}(:) - Lr{2}(:))), 0, 1e-12);
end
Rl = pt_compute_leadfields('bslaw', gd, {'k1'}, fullfile(gd, 'bs'));
n_fail = n_fail + check('existing file skipped on rerun', ...
    double(Rl.written == 0 && Rl.skipped == 1), 1, 0);

fprintf('\n%d failure(s)\n', n_fail);


function f = report(name, got, want, tol)
    ok = abs(got - want) <= tol;
    if ok, s = 'OK'; else, s = '*** FAIL ***'; end
    fprintf('%-42s got %12.6g  want %12.6g  %s\n', name, got, want, s);
    f = double(~ok);
end

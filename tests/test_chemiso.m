% ====================================================================== %
%  test_chemiso.m
%  Regression tests for core/chemisorption_eq.m against known-good
%  output of this exact shipped code (not against Rothschild et al.'s
%  independently published numbers -- that comparison lives in the
%  SoftwareX paper's validation table). A failure here means the core
%  solver's behavior changed, not that it disagrees with the reference
%  paper.
%
%  Run with:  matlab -batch "run('tests/test_chemiso.m')"
%  from the repository root, or cd into tests/ and run directly.
% ====================================================================== %

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','core'), fullfile(here,'..','presets'));
base = load_CdS_O2();

n_pass = 0; n_fail = 0;

% ---------------------------------------------------------------------
% Test 1: coverage conservation, Theta_tot = Theta_minus + Theta_zero
% Tolerance: this is an algebraic identity by construction
% (theta_minus = theta*fA_minus, theta_zero = theta*(1-fA_minus)), so
% any deviation above floating-point round-off (~1e-10) indicates a
% real bug, not numerical noise.
% ---------------------------------------------------------------------
Pset = logspace(-10,0,21);
[~, th_tot, th_m, th_0] = chemisorption_eq(base, Pset);
dev = max(abs(th_tot - (th_m+th_0)));
[n_pass,n_fail] = check('Conservation Theta=Theta^-+Theta^0', dev, 0, 1e-10, n_pass, n_fail);

% ---------------------------------------------------------------------
% Test 2: Vs_eq at the Fig. 3a reference point (P=1 atm, ND=1e16, T=300K)
% Expected value verified directly against this shipped code
% (v1.0.3, 400-point grid) on 2026-09-15: 0.688722 eV.
% Tolerance: 1e-4 eV -- checks reproducibility of THIS computation
% across MATLAB versions/platforms, not absolute grid accuracy (the
% grid itself only resolves Vs to ~1.2/399 = 0.003 eV, see
% test_grid_convergence.m for that separate question).
% ---------------------------------------------------------------------
par = base; par.ND = 1e16; par.T = 300;
Pset3a = logspace(-13,0,13);
[Vs_eq,~,~,~] = chemisorption_eq(par, Pset3a);
[n_pass,n_fail] = check('Fig.3a Vs_eq(P=1atm,ND=1e16)', Vs_eq(end), 0.688722, 1e-4, n_pass, n_fail);

% ---------------------------------------------------------------------
% Test 3: Fig. 4 slopes d(eV_s)/d(log10 P) at kT=300/400/500K.
% Expected values verified directly against this shipped code on
% 2026-09-15: 0.057744 / 0.075188 / 0.093233.
% Tolerance 1e-4: same reproducibility rationale as Test 2. Separately
% verified (test_grid_convergence.m) that these slopes are stable to
% <0.4% under 10x grid refinement, so this is also a reasonable
% physical-accuracy check, not just a regression one.
% ---------------------------------------------------------------------
T_list = [300 400 500];
expected_slopes = [0.057744, 0.075188, 0.093233];
Pset4 = logspace(-10,0,11); logP4 = log10(Pset4);
for k = 1:numel(T_list)
    par = base; par.T = T_list(k);
    [Vs_eq,~,~,~] = chemisorption_eq(par,Pset4);
    switch T_list(k)
        case 300, idx = (logP4>=-7)&(logP4<=-3);
        case 400, idx = (logP4>=-6)&(logP4<=-2);
        case 500, idx = (logP4>=-3)&(logP4<=-1);
    end
    p = polyfit(logP4(idx), Vs_eq(idx), 1);
    [n_pass,n_fail] = check(sprintf('Fig.4 slope T=%dK', T_list(k)), ...
        p(1), expected_slopes(k), 1e-4, n_pass, n_fail);
end

% ---------------------------------------------------------------------
% Test 4: Fig. 5b Theta^- saturation (P=1atm, T=300K), ND=1e14/1e16/1e18.
% Expected values verified directly against this shipped code
% (v1.0.3, 400-point grid) on 2026-09-15: 2.044647e-05 / 1.995532e-04 /
% 2.183927e-03. NOTE: these do NOT match the values previously printed
% in the paper's draft (2.05e-5/2.05e-4/2.18e-3) -- that discrepancy was
% found and the paper corrected; see 15-09-26-c105 in the collaboration
% bitacora. Theta^- depends exponentially on Vs_eq, so it is far more
% sensitive to grid resolution than Vs_eq itself (verified: -5.5%/+1.6%/
% -3.5% shift between 400 and 40000 grid points -- see
% test_grid_convergence.m). Tolerance here is loose (1%) precisely
% because this quantity is known to be numerically sensitive; it is a
% regression check against this exact shipped grid, not a claim of
% high absolute accuracy.
% ---------------------------------------------------------------------
Pset5 = logspace(-10,0,11);
ND_list = [1e14 1e16 1e18];
expected_thm = [2.044647e-05, 1.995532e-04, 2.183927e-03];
for i = 1:numel(ND_list)
    par = base; par.ND = ND_list(i); par.T = 300;
    [~,~,thm,~] = chemisorption_eq(par, Pset5);
    rel_tol = 0.01 * expected_thm(i);
    [n_pass,n_fail] = check(sprintf('Fig.5b Theta^- saturation ND=%.0e', ND_list(i)), ...
        thm(end), expected_thm(i), rel_tol, n_pass, n_fail);
end

fprintf('\n%d passed, %d failed.\n', n_pass, n_fail);
if n_fail > 0
    error('test_chemiso:failed', '%d regression test(s) failed.', n_fail);
end

function [n_pass,n_fail] = check(name, actual, expected, tol, n_pass, n_fail)
d = abs(actual - expected);
if d <= tol
    fprintf('[PASS] %s (actual=%.6g, expected=%.6g, |diff|=%.3g <= tol=%.3g)\n', ...
        name, actual, expected, d, tol);
    n_pass = n_pass + 1;
else
    fprintf('[FAIL] %s (actual=%.6g, expected=%.6g, |diff|=%.3g > tol=%.3g)\n', ...
        name, actual, expected, d, tol);
    n_fail = n_fail + 1;
end
end

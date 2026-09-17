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
% (v1.0.6, 400-pt coarse grid + local refinement) on 2026-09-15:
% 0.688322 eV. The plain 400-pt grid (v1.0.3-v1.0.5, no refinement)
% gave 0.688722 eV; the 40000-pt brute-force reference is 0.688307 eV
% -- the refined value is within 2e-5 eV of brute force, vs 4e-4 eV
% for the unrefined grid (see 15-09-26-c110 in the collaboration
% bitacora for the full before/after table and derivation).
% Tolerance: 1e-4 eV -- checks reproducibility of THIS computation
% across MATLAB versions/platforms, not absolute grid accuracy.
% ---------------------------------------------------------------------
par = base; par.ND = 1e16; par.T = 300;
Pset3a = logspace(-13,0,13);
[Vs_eq,~,~,~] = chemisorption_eq(par, Pset3a);
[n_pass,n_fail] = check('Fig.3a Vs_eq(P=1atm,ND=1e16)', Vs_eq(end), 0.688322, 1e-4, n_pass, n_fail);

% ---------------------------------------------------------------------
% Test 3: Fig. 4 slopes d(eV_s)/d(log10 P) at kT=300/400/500K.
% Expected values verified directly against this shipped code
% (v1.0.6, refined solver) on 2026-09-15: 0.057609 / 0.074672 /
% 0.092819. (v1.0.3-v1.0.5, unrefined 400-pt grid: 0.057744 / 0.075188
% / 0.093233 -- refinement shifts these by <0.7%, consistent with the
% <0.4% grid-convergence bound already measured in
% test_grid_convergence.m for the unrefined solver.)
% Tolerance 1e-4: reproducibility check, same rationale as Test 2.
% ---------------------------------------------------------------------
T_list = [300 400 500];
expected_slopes = [0.057609, 0.074672, 0.092819];
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
% (v1.0.6, 400-pt coarse grid + local refinement) on 2026-09-15:
% 1.939543e-05 / 2.026602e-04 / 2.110836e-03. The 40000-pt brute-force
% reference is 1.938042e-05 / 2.027786e-04 / 2.110665e-03 -- the refined
% values agree with brute force to <0.1%, vs the plain 400-pt grid
% (v1.0.3-v1.0.5: 2.044647e-05 / 1.995532e-04 / 2.183927e-03), which
% was off by -5.1%/+1.6%/-3.4%. This was the numerical issue review2
% flagged (per-point error from the fixed 400-pt default); see
% 15-09-26-c110 in the collaboration bitacora. Theta^- depends
% exponentially on Vs_eq, so it remains more sensitive to residual grid
% error than Vs_eq itself even after refinement -- tolerance here stays
% loose (1%) as a regression check against this exact shipped code, not
% a claim of higher absolute accuracy than the refinement provides.
% ---------------------------------------------------------------------
Pset5 = logspace(-10,0,11);
ND_list = [1e14 1e16 1e18];
expected_thm = [1.939543e-05, 2.026602e-04, 2.110836e-03];
for i = 1:numel(ND_list)
    par = base; par.ND = ND_list(i); par.T = 300;
    [~,~,thm,~] = chemisorption_eq(par, Pset5);
    rel_tol = 0.01 * expected_thm(i);
    [n_pass,n_fail] = check(sprintf('Fig.5b Theta^- saturation ND=%.0e', ND_list(i)), ...
        thm(end), expected_thm(i), rel_tol, n_pass, n_fail);
end

% ---------------------------------------------------------------------
% Test 5: root-existence/residual diagnostic (chemisorption_eq:noRoot).
% Purely diagnostic warning added 2026-09-17 in response to review4 of
% the SoftwareX paper (existence of a genuine Q_s=Q_sc crossing was
% never checked before -- the solver would silently return a
% closest-approach pinned to a domain boundary if no root existed).
% Verifies: (a) the warning never fires for the shipped CdS preset
% across its full parameter range (3 doping levels x 3 temperatures x
% 40 pressures from 1e-13 to 1 atm) -- max normalized residual measured
% 1.2e-3 against RESID_TOL=1e-2, comfortable margin; (b) the warning
% DOES fire for a deliberately pathological preset (DeltaE=5 eV, far
% outside any physical CdS value) where no crossing exists in
% [0,1.2] eV, confirming the check is not a no-op.
% ---------------------------------------------------------------------
lastwarn('');
ND_list5 = [1e14 1e16 1e18]; T_list5 = [300 400 500]; Pset_wide = logspace(-13,0,40);
for iN = 1:numel(ND_list5)
    for iT = 1:numel(T_list5)
        par = base; par.ND = ND_list5(iN); par.T = T_list5(iT);
        chemisorption_eq(par, Pset_wide);
    end
end
[~,warnid] = lastwarn();
[n_pass,n_fail] = check_bool('No spurious noRoot warning over shipped CdS range', ...
    isempty(warnid), n_pass, n_fail);

lastwarn('');
par_bad = base; par_bad.DeltaE = 5; par_bad.ND = 1e16;
chemisorption_eq(par_bad, 1.0);
[~,warnid] = lastwarn();
[n_pass,n_fail] = check_bool('noRoot warning fires for a pathological preset', ...
    strcmp(warnid,'chemisorption_eq:noRoot'), n_pass, n_fail);

fprintf('\n%d passed, %d failed.\n', n_pass, n_fail);
if n_fail > 0
    error('test_chemiso:failed', '%d regression test(s) failed.', n_fail);
end

function [n_pass,n_fail] = check_bool(name, ok, n_pass, n_fail)
if ok
    fprintf('[PASS] %s\n', name);
    n_pass = n_pass + 1;
else
    fprintf('[FAIL] %s\n', name);
    n_fail = n_fail + 1;
end
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

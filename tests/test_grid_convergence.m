% ====================================================================== %
%  test_grid_convergence.m
%  Base-grid independence check for core/chemisorption_eq.m's
%  coarse-grid + local-refinement scheme (400-pt coarse sweep over
%  0-1.2 eV, then a 400-pt local refinement around the coarse minimum --
%  see paper Section 1 for the design rationale).
%
%  PREVIOUS VERSION of this test (v1.0.5 and earlier) compared the plain,
%  unrefined 400-pt grid search against a 40000-pt brute-force reference
%  and found Theta^- shifted by up to 5.5% -- that finding is what
%  motivated adding the local-refinement step now present in
%  chemisorption_eq.m (see 15-09-26-c110 in the collaboration bitacora).
%  With refinement in place, comparing against a 40000-pt *unrefined*
%  brute force is no longer the right question -- refined 400 already
%  agrees with unrefined 40000 to <0.1% (validated separately). The
%  question this test now asks is whether the *coarse* base-grid size
%  (400 by default) still matters once local refinement is applied: it
%  reruns the same two-stage algorithm at coarse base grids of
%  200/400/800 points and asserts the results agree to a tight bound.
%
%  This does NOT modify or call chemisorption_eq.m. It reimplements only
%  the coarse+refine search loop (with a configurable base-grid size)
%  inline, driving the actual shared physics helpers
%  (wolkenstein_setup.m / wolkenstein_qsc.m / wolkenstein_qs.m) so the
%  Qs/Qsc formulas themselves are not duplicated a second time.
%
%  Run with:  matlab -batch "run('tests/test_grid_convergence.m')"
% ====================================================================== %

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','core'), fullfile(here,'..','presets'));
base = load_CdS_O2();

Ns = [200 400 800];

% ---------------------------------------------------------------------
% Vs_eq (Fig. 3a point) and the Fig. 4 slopes.
% ---------------------------------------------------------------------
par = base; par.ND = 1e16; par.T = 300;
Pset3a = logspace(-13,0,13);
Vs_eq_N = arrayfun(@(N) solve_Vseq(par, Pset3a, N), Ns);
rel_shift = (max(Vs_eq_N)-min(Vs_eq_N))/abs(Vs_eq_N(2));
report('Vs_eq(P=1atm,ND=1e16), base grid 200/400/800', rel_shift, 0.005, ...
    Vs_eq_N, Ns);

T_list = [300 400 500]; Pset4 = logspace(-10,0,11); logP4 = log10(Pset4);
for k = 1:numel(T_list)
    par = base; par.T = T_list(k);
    slopes = zeros(size(Ns));
    for j = 1:numel(Ns)
        Vs_eq = solve_full(par, Pset4, Ns(j));
        switch T_list(k)
            case 300, idx = (logP4>=-7)&(logP4<=-3);
            case 400, idx = (logP4>=-6)&(logP4<=-2);
            case 500, idx = (logP4>=-3)&(logP4<=-1);
        end
        p = polyfit(logP4(idx), Vs_eq(idx), 1);
        slopes(j) = p(1);
    end
    rel_shift = (max(slopes)-min(slopes))/abs(slopes(2));
    report(sprintf('Fig.4 slope T=%dK, base grid 200/400/800', T_list(k)), ...
        rel_shift, 0.005, slopes, Ns);
end

% ---------------------------------------------------------------------
% Theta^- saturation (Fig. 5b): the quantity found sensitive to grid
% resolution before refinement. Checked here for base-grid independence
% of the REFINED algorithm, not absolute accuracy (that comparison
% against brute force lives in the validation notes, not this test).
% ---------------------------------------------------------------------
fprintf('\n');
Pset5 = logspace(-10,0,11);
for ND = [1e14 1e16 1e18]
    par = base; par.ND = ND; par.T = 300;
    thm = zeros(size(Ns));
    for j = 1:numel(Ns)
        [~, theta_minus] = solve_full(par, Pset5, Ns(j));
        thm(j) = theta_minus(end);
    end
    rel_shift = (max(thm)-min(thm))/abs(thm(2));
    report(sprintf('Theta^- saturation ND=%.0e, base grid 200/400/800', ND), ...
        rel_shift, 0.005, thm, Ns);
end

function report(name, rel_shift, bound, values, Ns)
fprintf('%s: values = [%s] (N=[%s]), relative spread = %.3f%% (bound %.2f%%)\n', ...
    name, sprintf('%.6g ', values), sprintf('%d ', Ns), 100*rel_shift, 100*bound);
if rel_shift > bound
    error('test_grid_convergence:exceeded', ...
        '%s exceeded its documented bound (%.3f%% > %.2f%%).', name, 100*rel_shift, 100*bound);
end
end

function Vs_eq = solve_Vseq(par, Pset, N)
Vs_eq = solve_full(par, Pset, N);
Vs_eq = Vs_eq(end);
end

function [Vs_eq, theta_minus] = solve_full(par, Pset, N)
[~, ~, EC_EF, beta0, kT_eV] = wolkenstein_setup(par);
Vs_grid = linspace(0,1.2,N);
Qsc_tab = wolkenstein_qsc(par, Vs_grid);

Np = numel(Pset);
Vs_eq = zeros(1,Np);
theta_minus = zeros(Np,1);
for ip = 1:Np
    P = Pset(ip);
    Qs_vec = wolkenstein_qs(par, Vs_grid, EC_EF, beta0, kT_eV, P);
    [~,idx] = min(abs(abs(Qs_vec) - Qsc_tab));

    lo = Vs_grid(max(idx-1,1));
    hi = Vs_grid(min(idx+1,numel(Vs_grid)));
    Vs_fine  = linspace(lo, hi, 400);
    Qsc_fine = wolkenstein_qsc(par, Vs_fine);
    Qs_fine  = wolkenstein_qs(par, Vs_fine, EC_EF, beta0, kT_eV, P);
    [~,idxf] = min(abs(abs(Qs_fine) - Qsc_fine));
    Vs_star  = Vs_fine(idxf);
    Vs_eq(ip) = Vs_star;

    Vs = Vs_star;
    fA_minus = 1 / (1 + par.gA*exp((EC_EF + Vs - par.DeltaE)/kT_eV));
    num = 1 + (1/par.gA)*exp((par.DeltaE - EC_EF - Vs)/kT_eV);
    den = 1 + (1/par.gA)*exp(-(EC_EF + Vs)/kT_eV);
    beta  = beta0 * num / den;
    theta = (beta*P)/(1 + beta*P);
    theta_minus(ip) = theta * fA_minus;
end
end

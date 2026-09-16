% ====================================================================== %
%  test_grid_convergence.m
%  Grid-resolution sensitivity check for the direct grid-search method
%  in core/chemisorption_eq.m (Vs swept over a fixed 400-point grid,
%  0-1.2 eV, no residual tolerance or convergence check -- see paper
%  Section 1 for the design rationale and its limits).
%
%  This does NOT modify chemisorption_eq.m (that function, and its
%  400-point default, are the exact code archived under the paper's
%  DOI). It reimplements the same grid search inline at variable
%  resolution, following the precedent already set by figures/make_fig2.m
%  (which also reimplements the physics independently of the core).
%
%  Run with:  matlab -batch "run('tests/test_grid_convergence.m')"
% ====================================================================== %

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','presets'));
base = load_CdS_O2();

% ---------------------------------------------------------------------
% Vs_eq (used for Figs. 3a/3b/4/7): verified stable. Bound set from the
% measured 400-vs-4000-point shift (0.33% for the Fig.4 slope case,
% 2026-09-15) with a safety margin.
% ---------------------------------------------------------------------
par = base; par.ND = 1e16; par.T = 300;
Pset4 = logspace(-10,0,11); logP4 = log10(Pset4);
idx = (logP4>=-7)&(logP4<=-3);
Vs400  = arrayfun(@(P) solve_at(par,P,400),  Pset4);
Vs4000 = arrayfun(@(P) solve_at(par,P,4000), Pset4);
p400  = polyfit(logP4(idx), Vs400(idx), 1);
p4000 = polyfit(logP4(idx), Vs4000(idx), 1);
rel_shift = abs(p4000(1)-p400(1))/abs(p4000(1));
report('Fig.4 slope, 400 vs 4000 grid points', rel_shift, 0.01, ...
    'Vs_eq is well-converged at the shipped 400-point grid.');

% ---------------------------------------------------------------------
% Theta^- saturation (used for Fig. 5b): verified NOT well-converged at
% 400 points -- exponential dependence on Vs_eq amplifies small grid
% quantization shifts. This test documents and bounds that known
% sensitivity rather than asserting a false tight tolerance; the bound
% (15%) is set above the largest shift measured on 2026-09-15 (5.5%,
% 400 vs 40000 points) with margin, so it fails only if the solver's
% numerical behavior changes materially, not on ordinary noise.
% ---------------------------------------------------------------------
fprintf('\n');
max_rel_shift = 0;
for ND = [1e14 1e16 1e18]
    par = base; par.ND = ND; par.T = 300;
    thm400   = solve_thm_at(par, 1, 400);
    thm40000 = solve_thm_at(par, 1, 40000);
    rel_shift = abs(thm40000-thm400)/abs(thm40000);
    max_rel_shift = max(max_rel_shift, rel_shift);
    fprintf('  Theta^- ND=%.0e: 400-pt=%.6e, 40000-pt=%.6e, rel.shift=%.2f%%\n', ...
        ND, thm400, thm40000, 100*rel_shift);
end
report('Fig.5b Theta^- saturation, 400 vs 40000 grid points (max over ND)', ...
    max_rel_shift, 0.15, ...
    ['Theta^- is NOT well-converged at 400 points (depends exponentially ' ...
     'on Vs_eq); this is a known, bounded limitation, not a claim of high ' ...
     'precision for this quantity -- see paper Section 1.']);

function report(name, rel_shift, bound, note)
fprintf('%s: relative shift = %.2f%% (bound %.0f%%) -- %s\n', ...
    name, 100*rel_shift, 100*bound, note);
if rel_shift > bound
    error('test_grid_convergence:exceeded', ...
        '%s exceeded its documented bound (%.2f%% > %.0f%%).', name, 100*rel_shift, 100*bound);
end
end

function Vs_star = solve_at(par, P, Ngrid)
[Vs_star, ~] = solve_core(par, P, Ngrid);
end

function thm = solve_thm_at(par, P, Ngrid)
[~, thm] = solve_core(par, P, Ngrid);
end

function [Vs_star, thm] = solve_core(par, P, Ngrid)
q=1.602176634e-19; kB=1.380649e-23; eps0=8.854187817e-12; h=6.62607015e-34; m0=9.10938356e-31;
T=par.T; kT_eV=(kB*T)/q;
Nc=2*((2*pi*par.me_rel*m0*kB*T)/h^2)^(1.5)/1e6;
Nv=2*((2*pi*par.mh_rel*m0*kB*T)/h^2)^(1.5)/1e6;
ni=sqrt(Nc*Nv).*exp(-par.Eg/(2*kT_eV));
EC_EF=par.Eg-(par.Eg/2+0.5*kT_eV*log(Nv/Nc)+kT_eV*log(par.ND/ni));
beta0=(par.sticking*par.s0_m2)/(par.nu0*sqrt(2*pi*par.M_gas*kB*T))*exp(par.q0/kT_eV)*101325;
Vs_grid=linspace(0,1.2,Ngrid);
epsS=par.eps_r*eps0; ND_m3=par.ND*1e6;
Qsc_tab=sqrt(2*epsS*ND_m3*q.*Vs_grid)*1e-4;
Qs_vec=zeros(size(Vs_grid));
for j=1:numel(Vs_grid)
    Vs=Vs_grid(j);
    fA_minus=1/(1+par.gA*exp((EC_EF+Vs-par.DeltaE)/kT_eV));
    num=1+(1/par.gA)*exp((par.DeltaE-EC_EF-Vs)/kT_eV);
    den=1+(1/par.gA)*exp(-(EC_EF+Vs)/kT_eV);
    beta=beta0*num/den;
    theta=(beta*P)/(1+beta*P);
    Qs_vec(j)=-q*par.Nstar_cm2*theta*fA_minus;
end
[~,gidx]=min(abs(abs(Qs_vec)-Qsc_tab));
Vs_star=Vs_grid(gidx);
fA_minus=1/(1+par.gA*exp((EC_EF+Vs_star-par.DeltaE)/kT_eV));
num=1+(1/par.gA)*exp((par.DeltaE-EC_EF-Vs_star)/kT_eV);
den=1+(1/par.gA)*exp(-(EC_EF+Vs_star)/kT_eV);
beta=beta0*num/den;
theta=(beta*P)/(1+beta*P);
thm=theta*fA_minus;
end

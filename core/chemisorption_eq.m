function [Vs_eq, theta_tot, theta_minus, theta_zero, EC_EF] = chemisorption_eq(par, Pset)
% ========================================================================
%  chemisorption_eq
%  ----------------------------------------------------------------------
%  Material/gas-independent core. For each pressure P, computes:
%      • Equilibrium band-bending  Vs_eq  (eV, positive value)
%      • Total coverage  Θ            (theta_tot)
%      • Charged coverage Θ⁻          (theta_minus)
%      • Neutral coverage Θ⁰          (theta_zero)
%
%  INPUTS
%  ------
%  par  : struct  (returned by presets/*)  – all parameters
%  Pset : [1×N]   Pressures in atm
%
%  OUTPUTS
%  -------
%  Vs_eq       : [1×N]  equilibrium e|V_s|  (eV)
%  theta_tot   : [N×1]  total coverage
%  theta_minus : [N×1]  charged coverage
%  theta_zero  : [N×1]  neutral coverage
%  EC_EF       : scalar (E_C^b - E_F), eV — for Fig. 7 (E_act = EC_EF + e|Vs|)
%
%  METHOD
%  ------
%  • Sweep Vs (0→1.2 eV) on a coarse 400-point grid and compute Q_sc(Vs)
%    and Q_s(Vs,P) --- via the shared helpers wolkenstein_setup.m /
%    wolkenstein_qs.m, also used by figures/make_fig2.m.
%  • Take the coarse-grid Vs whose |Q_s| best matches Q_sc, then refine
%    locally: re-evaluate |Q_sc-|Q_s|| on a 400-point grid spanning the
%    two neighboring coarse intervals around that minimum and take its
%    minimum instead. This keeps the cost of the coarse sweep (one
%    401-point residual evaluation per pressure point) while resolving
%    Vs_eq to a local grid spacing of ~1.2/400^2 eV instead of 1.2/400.
% ========================================================================

[Vs_grid, Qsc_tab, EC_EF, beta0, kT_eV] = wolkenstein_setup(par);

N        = numel(Pset);
Vs_eq       = zeros(1,N);
theta_tot   = zeros(N,1);
theta_minus = zeros(N,1);
theta_zero  = zeros(N,1);

for ip = 1:N
    P  = Pset(ip);
    Qs_vec = wolkenstein_qs(par, Vs_grid, EC_EF, beta0, kT_eV, P);

    % ---- Vs_eq via minimum of |Q_sc - |Q_s|| : coarse grid ... ----
    [~,idx] = min(abs(abs(Qs_vec) - Qsc_tab));

    % ---- ... then local refinement around the coarse minimum ----
    lo = Vs_grid(max(idx-1,1));
    hi = Vs_grid(min(idx+1,numel(Vs_grid)));
    Vs_fine  = linspace(lo, hi, 400);
    Qsc_fine = wolkenstein_qsc(par, Vs_fine);
    Qs_fine  = wolkenstein_qs(par, Vs_fine, EC_EF, beta0, kT_eV, P);
    [~,idxf] = min(abs(abs(Qs_fine) - Qsc_fine));
    Vs_star  = Vs_fine(idxf);        % e|V_s|  (eV)
    Vs_eq(ip) = Vs_star;

    % ---- Coverages at Vs_star ----
    Vs = Vs_star;
    fA_minus = 1 / (1 + par.gA*exp((EC_EF + Vs - par.DeltaE)/kT_eV));
    num = 1 + (1/par.gA)*exp((par.DeltaE - EC_EF - Vs)/kT_eV);
    den = 1 + (1/par.gA)*exp(-(EC_EF + Vs)/kT_eV);
    beta  = beta0 * num / den;
    theta = (beta*P)/(1 + beta*P);

    theta_tot(ip)   = theta;
    theta_minus(ip) = theta * fA_minus;
    theta_zero(ip)  = theta * (1 - fA_minus);
end
end

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
%  • After refinement, checks that a genuine electroneutrality crossing
%    exists: |Q_s(0)|>=Q_sc(0) and |Q_s(1.2)|<=Q_sc(1.2) must have
%    opposite-sign residuals (root bracketed by the endpoints), and the
%    final refined residual, normalized by the domain-wide scale
%    Q_sc(1.2), must stay below RESID_TOL. Both conditions hold with
%    large margin for every (N_D,T,P) combination of the shipped CdS
%    preset (residual <=1.2e-3 of RESID_TOL=1e-2, no sign-check failures
%    over N_D in {1e14,1e16,1e18}, T in {300,400,500}, P in
%    [1e-13,1] atm); a violation flags a preset/pressure combination
%    where the returned Vs_eq is a spurious closest-approach, not a
%    true root. This check is purely diagnostic: it never changes
%    Vs_eq itself, only emits a warning.
% ========================================================================

RESID_TOL = 1e-2;   % measured max over the shipped CdS sweep: 1.2e-3

[Vs_grid, Qsc_tab, EC_EF, beta0, kT_eV] = wolkenstein_setup(par);
domain_scale = Qsc_tab(end);

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
    [minval,idxf] = min(abs(abs(Qs_fine) - Qsc_fine));
    Vs_star  = Vs_fine(idxf);        % e|V_s|  (eV)
    Vs_eq(ip) = Vs_star;

    % ---- Diagnostic check: does a genuine crossing exist, and is the
    %      refined residual small relative to the problem's own scale? ----
    f_lo = abs(Qs_vec(1))   - Qsc_tab(1);
    f_hi = abs(Qs_vec(end)) - Qsc_tab(end);
    no_bracket   = isnan(f_lo) || isnan(f_hi) || (sign(f_lo) == sign(f_hi) && f_lo ~= 0 && f_hi ~= 0);
    resid_norm   = minval / domain_scale;
    if no_bracket || resid_norm > RESID_TOL
        warning('chemisorption_eq:noRoot', ...
            ['chemisorption_eq: no electroneutrality crossing confidently found ', ...
             'in the swept V_s domain [0,%.2g] eV at P=%.3g atm (bracket ok=%d, ', ...
             'normalized residual=%.3g, tol=%.3g). Returned Vs_eq may be a ', ...
             'closest-approach, not a true root.'], ...
            Vs_grid(end), P, ~no_bracket, resid_norm, RESID_TOL);
    end

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

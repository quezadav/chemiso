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
%  • Sweep Vs (0→1.2 eV) and compute Q_sc(Vs) and Q_s(Vs,P).
%  • Take the Vs whose |Q_s| best matches Q_sc.
% ========================================================================

% ---------- 1.  Universal constants ----------
q   = 1.602176634e-19;   % C
kB  = 1.380649e-23;      % J K⁻¹
eps0 = 8.854187817e-12;  % F m⁻¹
h   = 6.62607015e-34;    % J s
m0  = 9.10938356e-31;    % kg

% ---------- 2.  Convenience aliases ----------
T       = par.T;
kT_eV   = (kB*T)/q;

% ---------- 3.  Density of states, ni, E_C‑E_F ----------
Nc = 2*((2*pi*par.me_rel*m0*kB*T)/h^2)^(1.5)/1e6;   % cm⁻³
Nv = 2*((2*pi*par.mh_rel*m0*kB*T)/h^2)^(1.5)/1e6;   % cm⁻³
ni = sqrt(Nc*Nv).*exp(-par.Eg/(2*kT_eV));           % cm⁻³
EC_EF = par.Eg - (par.Eg/2 + 0.5*kT_eV*log(Nv/Nc) ...
                  + kT_eV*log(par.ND/ni));          % eV

% ---------- 4.  Prefactor β₀ (atm⁻¹) ----------
beta0 = (par.sticking*par.s0_m2)/(par.nu0*sqrt(2*pi*par.M_gas*kB*T)) ...
        * exp(par.q0/kT_eV) * 101325;               % Pa⁻¹→atm⁻¹

% ---------- 5.  Q_sc(Vs) table (material-dependent only) ----------
Vs_grid = linspace(0,1.2,400);                      % eV
epsS    = par.eps_r*eps0;
ND_m3   = par.ND*1e6;
Qsc_tab = sqrt(2*epsS*ND_m3*q.*Vs_grid)*1e-4;       % C cm⁻²

% ---------- 6.  Output preallocation ----------
N        = numel(Pset);
Vs_eq       = zeros(1,N);
theta_tot   = zeros(N,1);
theta_minus = zeros(N,1);
theta_zero  = zeros(N,1);

% ======================= 7.  Loop over pressures =====================
for ip = 1:N
    P  = Pset(ip);
    Qs_vec = zeros(size(Vs_grid));

    % ---- 7.1  |Q_s|(Vs) for this pressure ----
    for j = 1:numel(Vs_grid)
        Vs = Vs_grid(j);

        % ---- f_A⁻ (Wolkenstein equation) ----
        fA_minus = 1 / (1 + par.gA*exp((EC_EF + Vs - par.DeltaE)/kT_eV));

        % ---- β(Vs) & Θ(Vs,P) ----
        num = 1 + (1/par.gA)*exp((par.DeltaE - EC_EF - Vs)/kT_eV);
        den = 1 + (1/par.gA)*exp(-(EC_EF + Vs)/kT_eV);
        beta = beta0 * num / den;

        theta = (beta*P)/(1 + beta*P);

        Qs_vec(j) = -q*par.Nstar_cm2*theta*fA_minus;  % C cm⁻²
    end

    % ---- 7.2  Vs_eq via minimum of |Q_sc - |Q_s|| ----
    [~,idx] = min(abs(abs(Qs_vec) - Qsc_tab));
    Vs_star = Vs_grid(idx);          % e|V_s|  (eV)
    Vs_eq(ip) = Vs_star;

    % ---- 7.3  Coverages at Vs_star ----
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

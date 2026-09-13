function par = load_CdS_O2()
% ===============================================================
%  Parameters for CdS / O2 (based on Rothschild et al., 2002)
%  -> Returns struct 'par' that feeds chemisorption_eq.
% ===============================================================

% -- Temperature & doping --
par.T   = 300;      % K
par.ND  = 1e16;     % cm⁻³

% -- Band gap & effective masses --
par.Eg      = 2.42; % eV
par.me_rel  = 0.21; % m*/m0 (electrons)
par.mh_rel  = 0.80; % m*/m0 (holes)
par.eps_r   = 5.4;

% -- Surface energies --
par.DeltaE = 0.80;  % eV (E_C^s - E_A^s)
par.gA     = 2;     % degeneracy factor

% -- Neutral adsorption --
par.q0       = 0.10;   % eV
par.s0_m2    = 1e-19;  % m²
par.sticking = 1;      % σ
par.nu0      = 1e13;   % s⁻¹
par.nu_ion   = 1e13;   % s⁻¹
par.Nstar_cm2= 1e15;   % cm⁻²

% -- Gas --
par.M_gas = 32*1.66054e-27; % kg
end

function [Vs_grid, Qsc_tab, EC_EF, beta0, kT_eV] = wolkenstein_setup(par)
% ========================================================================
%  wolkenstein_setup
%  ----------------------------------------------------------------------
%  Material/gas-dependent but pressure-independent quantities shared by
%  chemisorption_eq.m (core solver) and figures/make_fig2.m: the surface
%  band-bending sweep Vs_grid, the space-charge curve Qsc(Vs), the bulk
%  Fermi-level position (E_C^b-E_F), and the isotherm prefactor beta0.
%
%  INPUTS
%  ------
%  par : struct  (returned by presets/*) -- all parameters
%
%  OUTPUTS
%  -------
%  Vs_grid : [1x400]  surface band-bending sweep, e|Vs| (eV), 0 to 1.2
%  Qsc_tab : [1x400]  space-charge curve Qsc(Vs) (C/cm^2)
%  EC_EF   : scalar   bulk (E_C^b - E_F) (eV)
%  beta0   : scalar   isotherm prefactor beta0 (atm^-1)
%  kT_eV   : scalar   kB*T in eV
% ========================================================================

q    = 1.602176634e-19;   % C
kB   = 1.380649e-23;      % J K^-1
eps0 = 8.854187817e-12;   % F m^-1
h    = 6.62607015e-34;    % J s
m0   = 9.10938356e-31;    % kg

T     = par.T;
kT_eV = (kB*T)/q;

Nc = 2*((2*pi*par.me_rel*m0*kB*T)/h^2)^(1.5)/1e6;   % cm^-3
Nv = 2*((2*pi*par.mh_rel*m0*kB*T)/h^2)^(1.5)/1e6;   % cm^-3
ni = sqrt(Nc*Nv).*exp(-par.Eg/(2*kT_eV));           % cm^-3
EC_EF = par.Eg - (par.Eg/2 + 0.5*kT_eV*log(Nv/Nc) ...
                  + kT_eV*log(par.ND/ni));          % eV

beta0 = (par.sticking*par.s0_m2)/(par.nu0*sqrt(2*pi*par.M_gas*kB*T)) ...
        * exp(par.q0/kT_eV) * 101325;               % Pa^-1 -> atm^-1

Vs_grid = linspace(0,1.2,400);                      % eV
epsS    = par.eps_r*eps0;
ND_m3   = par.ND*1e6;
Qsc_tab = sqrt(2*epsS*ND_m3*q.*Vs_grid)*1e-4;       % C cm^-2
end

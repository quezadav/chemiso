function Qsc = wolkenstein_qsc(par, Vs)
% ========================================================================
%  wolkenstein_qsc
%  ----------------------------------------------------------------------
%  Space-charge curve Qsc(Vs) (Poisson/depletion-region term). Shared by
%  wolkenstein_setup.m (coarse grid tabulation) and chemisorption_eq.m's
%  local refinement step (evaluated again on a fine local Vs window),
%  so the closed-form Qsc(Vs) formula lives in exactly one place.
%
%  INPUTS
%  ------
%  par : struct   (returned by presets/*)
%  Vs  : [1xN]    surface band-bending values (eV)
%
%  OUTPUT
%  ------
%  Qsc : [1xN]  space-charge curve (C/cm^2)
% ========================================================================

q    = 1.602176634e-19;   % C
eps0 = 8.854187817e-12;   % F m^-1

epsS  = par.eps_r*eps0;
ND_m3 = par.ND*1e6;
Qsc   = sqrt(2*epsS*ND_m3*q.*Vs)*1e-4;   % C cm^-2
end

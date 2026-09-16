function Qs_curve = wolkenstein_qs(par, Vs_grid, EC_EF, beta0, kT_eV, P)
% ========================================================================
%  wolkenstein_qs
%  ----------------------------------------------------------------------
%  Surface charge Qs(Vs,P) over a Vs sweep, at a single pressure P.
%  Shared by chemisorption_eq.m (core solver) and figures/make_fig2.m --
%  both need the full Qs(Vs) curve (chemisorption_eq.m to search it for
%  the electroneutrality crossing, make_fig2.m to plot it directly).
%
%  INPUTS
%  ------
%  par     : struct   (returned by presets/*)
%  Vs_grid : [1xN]    surface band-bending sweep (eV), from wolkenstein_setup
%  EC_EF   : scalar    (E_C^b - E_F) (eV), from wolkenstein_setup
%  beta0   : scalar    isotherm prefactor (atm^-1), from wolkenstein_setup
%  kT_eV   : scalar    kB*T (eV), from wolkenstein_setup
%  P       : scalar    pressure (atm)
%
%  OUTPUT
%  ------
%  Qs_curve : [1xN]  surface charge Qs(Vs,P) over Vs_grid (C/cm^2)
% ========================================================================

q = 1.602176634e-19;   % C

Qs_curve = zeros(size(Vs_grid));
for j = 1:numel(Vs_grid)
    Vs = Vs_grid(j);

    fA_minus = 1 / (1 + par.gA*exp((EC_EF + Vs - par.DeltaE)/kT_eV));

    num = 1 + (1/par.gA)*exp((par.DeltaE - EC_EF - Vs)/kT_eV);
    den = 1 + (1/par.gA)*exp(-(EC_EF + Vs)/kT_eV);
    beta = beta0 * num / den;

    theta = (beta*P)/(1 + beta*P);

    Qs_curve(j) = -q*par.Nstar_cm2*theta*fA_minus;
end
end

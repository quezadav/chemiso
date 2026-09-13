% ====================================================================== %
%  make_fig7.m  -  Fig. 7 of Rothschild et al.                           %
%  Comparison between the computed curve e|Vs^0|+0.21 eV (for n-type    %
%  CdS with ND=10^16 cm^-3 at T=400K) and the experimental data from     %
%  Ref. 19 cited in the paper.                                           %
%  NOTE: the previous version of this file plotted Vs vs ND at fixed P   %
%  (P=1e-6, T=300K) - that does NOT correspond to the real Fig. 7 of the %
%  paper; that real figure is Activation energy vs PRESSURE at fixed     %
%  ND and T.                                                              %
%  Fixed 2026-09-13 after comparing against page 7095 of the PDF.        %
% ====================================================================== %
function make_fig7(base)
if nargin == 0
    addpath(fullfile(pwd,'core'), fullfile(pwd,'presets'));
    base = load_CdS_O2();  % CdS / O2 preset
end

par = base; par.ND = 1e16; par.T = 400;   % same as Fig. 7 of the paper
Pset = logspace(-13,-1,25);               % atm
[Vs_eq, ~, ~, ~, EC_EF] = chemisorption_eq(par, Pset);
% E_act = e|Vs^0| + (E_C^b-E_F); EC_EF is computed in chemisorption_eq.m
% with this preset's ND/T (~0.204 eV, matching the ~0.21 eV cited in the
% paper's text for T=400K) - a fixed constant is no longer used, so the
% figure stays correct if ND or T are changed above.
Eact = Vs_eq + EC_EF;

figure; hold on; box on;
set(gca,'XScale','log');
plot(Pset, Eact, '--o', 'LineWidth', 1.6, 'MarkerSize', 5);
xlabel('P (atm)');
ylabel('Activation Energy (eV)');
title('Fig. 7 - e|V_s^0|+0.21 eV vs P  (N_D=10^{16} cm^{-3}, T=400K)');
ylim([0 0.9]); grid on;
end

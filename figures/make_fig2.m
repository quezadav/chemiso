% ====================================================================== %
%  make_fig2.m                                                           %
%  ------------------------------------------------------------------   %
%  Reproduces Fig. 2 of Rothschild et al. (2002):                        %
%      |Q_s| and Q_sc vs e|V_s| for O2/CdS (ND = 1·10^16 cm‑3, T = 300 K) %
%                                                                       %
%  ▸ Requires:                                                          %
%        • load_CdS_O2.m       (preset with all parameters)             %
%        • wolkenstein_setup.m (shared: Vs sweep, Qsc(Vs), EC_EF, beta0) %
%        • wolkenstein_qs.m    (shared: Qs(Vs,P) curve)                 %
%                                                                       %
%  ▸ What it does:                                                      %
%        1) Computes Q_sc (Poisson) over a sweep of Vs                  %
%        2) Computes Q_s(P,Vs) for P = 10^{‑10} … 1 atm                 %
%        3) Finds the crossings |Q_s| = Q_sc and prints them to screen  %
%        4) Plots all curves in the style of the original article       %
% ====================================================================== %

function make_fig2(par)
% -----------------------------------------------------------------------
%  0) SHARED SETUP (Vs sweep, Qsc(Vs), EC_EF, beta0) -- same core physics
%     used by chemisorption_eq.m for every other figure.
% -----------------------------------------------------------------------
[Vs_eV, Qsc, EC_EF, beta0, kT_eV] = wolkenstein_setup(par);

% -----------------------------------------------------------------------
%  1) Q_s(P,Vs) CALCULATION FOR SEVERAL PRESSURES
% -----------------------------------------------------------------------
Pset = [1e-10 1e-8 1e-6 1e-4 1e-2 1];              % atm
Qs   = zeros(numel(Pset), numel(Vs_eV));

for ip = 1:numel(Pset)
    Qs(ip,:) = wolkenstein_qs(par, Vs_eV, EC_EF, beta0, kT_eV, Pset(ip));
end

% -----------------------------------------------------------------------
%  2) PRINT CROSSINGS  |Q_s| = Q_sc  (sensor-relevant info)
% -----------------------------------------------------------------------
fprintf('\n--- Crossings |Q_s| = Q_sc  (Fig. 2)  -----------------------------\n');
for ip = 1:numel(Pset)
    [~, idx] = min(abs(abs(Qs(ip,:)) - Qsc));   % crossing index
    fprintf('P = %.1e atm  ->  e|V_s| ~ %.4f eV\n', Pset(ip), Vs_eV(idx));
end
fprintf('----------------------------------------------------------------\n\n');

% -----------------------------------------------------------------------
%  3) PLOT:  |Q_s| and Q_sc vs e|V_s|
% -----------------------------------------------------------------------

set(gca,'YScale','log','FontSize',10)
hold on; box on
plot(Vs_eV, Qsc, 'k', 'LineWidth',2)               % Q_sc (black)

cmap = lines(numel(Pset));
for ip = 1:numel(Pset)
    plot(Vs_eV, abs(Qs(ip,:)), 'Color',cmap(ip,:), 'LineWidth',1.6);
end

xlabel('Surface band-bending  e|V_s|  (eV)')
ylabel('Charge density  |Q_s|; Q_{sc}  (C/cm^2)')

    legend( ['Q_{sc}', arrayfun(@(p)sprintf('P = 10^{%d} atm',log10(p)), Pset,'uni',0)], ...
            'Location','best')
    title(sprintf('Fig. 2  -  O_2/CdS  (T=%d K,  N_D=10^{%d} cm^{-3})', ...
          par.T,round(log10(par.ND))));
grid on
xlim([0 1.2]); ylim([1e-10 1e-2]);
end

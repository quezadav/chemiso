% ====================================================================== %
%  make_fig2.m                                                           %
%  ------------------------------------------------------------------   %
%  Reproduces Fig. 2 of Rothschild et al. (2002):                        %
%      |Q_s| and Q_sc vs e|V_s| for O2/CdS (ND = 1·10^16 cm‑3, T = 300 K) %
%                                                                       %
%  ▸ Requires:                                                          %
%        • load_CdS_O2.m   (preset with all parameters)                 %
%        • chemisorption_eq.m  (core that computes Qs, Vs, coverages)   %
%                                                                       %
%  ▸ What it does:                                                      %
%        1) Computes Q_sc (Poisson) over a sweep of Vs                  %
%        2) Computes Q_s(P,Vs) for P = 10^{‑10} … 1 atm                 %
%        3) Finds the crossings |Q_s| = Q_sc and prints them to screen  %
%        4) Plots all curves in the style of the original article       %
% ====================================================================== %

function make_fig2(par)
% -----------------------------------------------------------------------
%  0) CONSTANTS AND BASIC PARAMETERS
% -----------------------------------------------------------------------
q     = 1.602176634e-19;         % C
kB    = 1.380649e-23;            % J·K⁻¹
T     = par.T;                   % K  (taken from preset)
kT_eV = kB*T/q;                  % eV
eps0  = 8.854187817e-12;         % F·m⁻¹
h     = 6.62607015e-34;          % J·s
m0    = 9.10938356e-31;          % kg
M_O2  = 32*1.6605390666e-27;     % kg  O2 molecule mass

% -- density of states (cm⁻³) and ni -------------------------------------
Nc = 2*((2*pi*par.me_rel*m0*kB*T)/h^2)^(1.5)/1e6;
Nv = 2*((2*pi*par.mh_rel*m0*kB*T)/h^2)^(1.5)/1e6;
ni = sqrt(Nc*Nv).*exp(-par.Eg/(2*kT_eV));

% -- EF position relative to EC (in eV) ----------------------------------
EC_EF = par.Eg - (par.Eg/2 + 0.5*kT_eV*log(Nv/Nc) + ...
                  kT_eV*log(par.ND/ni));

% -- prefactor β₀  (Pa⁻¹ → atm⁻¹) ---------------------------------------
beta0 = (par.sticking*par.s0_m2)/(par.nu0*sqrt(2*pi*M_O2*kB*T)) ...
       * exp(par.q0/kT_eV)*101325;   % ← 101325 converts Pa→atm

% -----------------------------------------------------------------------
%  1) Vs SWEEP  and  Q_sc
% -----------------------------------------------------------------------
Vs_eV = linspace(0,1.2,400);                       % e|V_s| (eV)
Qsc   = sqrt(2*par.eps_r*eps0*par.ND*1e6*q).*sqrt(Vs_eV)*1e-4;  % C·cm⁻²

% -----------------------------------------------------------------------
%  2) Q_s(P,Vs) CALCULATION FOR SEVERAL PRESSURES
% -----------------------------------------------------------------------
Pset = [1e-10 1e-8 1e-6 1e-4 1e-2 1];              % atm
Qs   = zeros(numel(Pset), numel(Vs_eV));

for ip = 1:numel(Pset)
    P = Pset(ip);
    for j = 1:numel(Vs_eV)
        Vs = Vs_eV(j);

        % -- fraction of charged acceptors (f_A⁻) ----------------------
        expo = (EC_EF + Vs - par.DeltaE)/kT_eV;
        fAminus = 1 / (1 + par.gA*exp(expo));

        % -- coefficient β(Vs) (Wolkenstein) ---------------------------
        num = 1 + (1/par.gA)*exp((par.DeltaE - EC_EF - Vs)/kT_eV);
        den = 1 + (1/par.gA)*exp(-(EC_EF + Vs)/kT_eV);
        beta = beta0 * num / den;

        % -- total coverage θ -------------------------------------------
        theta = (beta*P)/(1 + beta*P);

        % -- surface charge density Q_s ---------------------------------
        Qs(ip,j) = -q*par.Nstar_cm2*theta*fAminus;   % C·cm⁻²
    end
end

% -----------------------------------------------------------------------
%  3) PRINT CROSSINGS  |Q_s| = Q_sc  (sensor-relevant info)
% -----------------------------------------------------------------------
fprintf('\n--- Crossings |Q_s| = Q_sc  (Fig. 2)  -----------------------------\n');
for ip = 1:numel(Pset)
    [~, idx] = min(abs(abs(Qs(ip,:)) - Qsc));   % crossing index
    fprintf('P = %.1e atm  ->  e|V_s| ~ %.4f eV\n', Pset(ip), Vs_eV(idx));
end
fprintf('----------------------------------------------------------------\n\n');

% -----------------------------------------------------------------------
%  4) PLOT:  |Q_s| and Q_sc vs e|V_s|
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

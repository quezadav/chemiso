% ====================================================================== %
%  make_fig4.m  -  Fig. 4 of Rothschild et al.                           %
% ====================================================================== %
function make_fig4(base)
if nargin==0
    addpath(fullfile(pwd,'core'), fullfile(pwd,'presets'));
    base = load_CdS_O2();           % CdS / O2 preset
end

T_list = [300 400 500];
q      = 1.602176634e-19;           % C
kB     = 1.380649e-23;              % J K^-1
kT_eV  = kB.*T_list./q;

Pset   = logspace(-10,0,11);        % 1e‑10 … 1 atm
logP   = log10(Pset);

S = zeros(size(T_list));

for k = 1:numel(T_list)
    par = base; par.T = T_list(k);
    [Vs_eq,~,~,~] = chemisorption_eq(par,Pset);   % |Vs^0|(eV)

    % --- material-specific linear segment ------------------------------
    switch par.T
        case 300, idx = (logP>=-7)&(logP<=-3);   % 10^-7-10^-3
        case 400, idx = (logP>=-6)&(logP<=-2);   % 10^-6-10^-2
        case 500, idx = (logP>=-3)&(logP<=-1);   % 10^-3-10^-1  (fit)
        otherwise, error('T=%g K not handled',par.T)
    end

    p   = polyfit(logP(idx), Vs_eq(idx),1);       % slope
    S(k)= p(1);
end

% --------------------------- plot -----------------------------------

plot(kT_eV,S,'-o','LineWidth',1.8,'MarkerSize',6),hold on; box on
xlabel('kT  (eV)')
ylabel('S(T) = d(e|V_s^0|)/d(log_{10}P)  (eV)')
title('Fig. 4  -  Slope d(e|V_s^0|)/d(log_{10}P) vs kT')
xlim([0.024 0.045]); ylim([0.05 0.10]); grid on
end

function make_fig3a(base)
% -------------------------------------------------------------
%  Figure 3a – Vs_eq vs P for ND = 1e14, 1e16, 1e18 cm⁻³
%  (T = 300 K) – Rothschild Fig. 3(a)-style curve
% -------------------------------------------------------------
Pset   = logspace(-13,0,13);
NDlist = [1e14 1e16 1e18];
clr    = lines(numel(NDlist));

hold on; box on;
for i = 1:numel(NDlist)
    par = base;
    par.ND = NDlist(i);
    par.T  = 300;
    [Vs_eq,~,~,~] = chemisorption_eq(par,Pset);
    semilogx(Pset, Vs_eq, 'LineWidth',1.6,'Color',clr(i,:))
end
xlabel('P_{O_2}  (atm)')
ylabel('e|V_s|_{eq}  (eV)')
title('Fig. 3a - e|V_s|_{eq} vs P  (T = 300 K)')
legend(arrayfun(@(n)sprintf('N_D=10^{%d}',log10(n)),NDlist,'uni',0), ...
       'Location','best')
set(gca,'XScale','log','XMinorGrid','on'); xlim([1e-13 1]); grid on
end

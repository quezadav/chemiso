function make_fig3b(base)
% -------------------------------------------------------------
%  Figure 3b – Vs_eq vs P for T = 300, 400, 500 K
%  (N_D = 1e16 cm⁻³) – Rothschild Fig. 3(b)-style curve
% -------------------------------------------------------------
Pset  = logspace(-13,0,13);
Tlist = [300 400 500];
clr   = lines(numel(Tlist));

hold on; box on;
for i = 1:numel(Tlist)
    par = base;
    par.ND = 1e16;
    par.T  = Tlist(i);
    [Vs_eq,~,~,~] = chemisorption_eq(par,Pset);
    semilogx(Pset, Vs_eq, 'LineWidth',1.6,'Color',clr(i,:))
end
xlabel('P_{O_2}  (atm)')
ylabel('e|V_s|_{eq}  (eV)')
title('Fig. 3b - e|V_s|_{eq} vs P  (N_D = 10^{16} cm^{-3})')
legend(arrayfun(@(T)sprintf('%d K',T),Tlist,'uni',0),'Location','best')
set(gca,'XScale','log','XMinorGrid','on'); xlim([1e-13 1]); grid on
end

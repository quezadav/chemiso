% ====================================================================== %
%  fit_ZnAl2O4_titulacion.m - Tercer material: ZnAl2O4 pellets, prueba
%  ESTATICA (Guillen-Bonilla et al. 2021, Sensors 21, 2362, Fig. 10a).
%  Contexto: se uso este paper primero para intentar un discriminador A/B
%  (estatico vs dinamico, mismo material) - quedo INCONCLUSO (ver bitacora,
%  demasiados confundentes: concentracion real post-dilucion desconocida
%  en el protocolo dinamico, lectura dinamica no necesariamente en estado
%  estacionario). Pero la serie ESTATICA sola sirve como tercer punto de
%  comparacion para la pregunta "converge K o degenera" (CuSb2O6 vs GdCoO3).
%
%  OJO: salvo los dos anclajes exactos en texto (S=2.8 a 200C/500ppm,
%  S=7.8 a 300C/500ppm), el resto de la serie se leyo VISUALMENTE de la
%  Fig. 10a (no hay tabla numerica) - menor confianza que CuSb2O6/GdCoO3.
% ====================================================================== %

C = [1,5,50,100,200,300,400,500];
S_300 = [0.3,0.8,7.0,7.1,7.3,7.5,7.6,7.8];  % lectura visual + anclas exactas
S_200 = [0.1,0.3,1.7,2.0,2.3,2.4,2.6,2.8];  % lectura visual + anclas exactas

opts = optimset('Display','off','TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',20000,'MaxIter',20000);
modelWolk = @(x,P) exp(-x(1).*(x(2).*P./(1+x(2).*P))); % rho_rel(0)=1

results = struct();
figure; hold on; box on; set(gca,'YScale','log');
for k=1:2
    if k==1, S=S_300; T=573; lbl='ZnAl2O4_300C'; col='b'; else, S=S_200; T=473; lbl='ZnAl2O4_200C'; col='r'; end
    rho_rel = 1./(1+S);
    logdata = log(rho_rel);
    cost = @(x) sum((log(modelWolk(x,C))-logdata).^2);
    [xfit,~] = fminsearch(cost, [3,0.01], opts);
    rhofit = modelWolk(xfit,C);
    R2 = 1 - sum((log(rhofit)-logdata).^2)/sum((logdata-mean(logdata)).^2);
    kT_eV = 8.617e-5*T;
    results.(lbl) = struct('A',xfit(1),'K',xfit(2),'R2',R2,'DeltaVsMax',xfit(1)*kT_eV);
    fprintf('%-16s A=%.4f  K=%.5f/ppm  R2=%.4f  DeltaVsMax=%.4f eV\n', lbl, xfit(1), xfit(2), R2, xfit(1)*kT_eV);

    semilogy(C, rho_rel, 'o', 'Color', col, 'MarkerFaceColor', col);
    Cfine = linspace(1,500,200);
    semilogy(Cfine, modelWolk(xfit,Cfine), '-', 'Color', col);
end
legend({'300C (datos)','300C (ajuste)','200C (datos)','200C (ajuste)'},'Location','best');
xlabel('C_3H_8 (ppm)'); ylabel('\rho_G/\rho_0');
title('ZnAl2O4 (estatico): ajuste tipo Wolkenstein-titulacion NO degenera (3er material)');
grid on;

% CONCLUSION (ver bitacora): NO degenera (K finito, R2 alto), igual que
% GdCoO3 y a diferencia de CuSb2O6. K aqui es ~20-30x mayor que en GdCoO3
% (satura mucho antes, ~50ppm) pero la forma se sostiene igual.
% Con esto: 2 de 3 materiales (GdCoO3, ZnAl2O4) bien portados; CuSb2O6 es
% la anomalia, no el patron - la afirmacion mas solida disponible para
% un futuro paper, mas solida que cualquier conclusion sobre Escenario A/B.

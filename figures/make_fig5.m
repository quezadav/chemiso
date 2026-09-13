% ====================================================================== %
%  make_fig5.m  -  Fig. 5 of Rothschild et al.                           %
% ====================================================================== %
function make_fig5(base)
if nargin == 0
    addpath(fullfile(pwd,'core'), fullfile(pwd,'presets'));
    base = load_CdS_O2();  % CdS / O2 preset
end

Pset = logspace(-10, 0, 11);  % atm

%% ----------- (a) ND sweep at T = 300 K ---------------------------
figure; hold on; box on;
set(gca, 'XScale','log', 'YScale','log')
title('Fig. 5a - Coverages vs P   (T = 300 K)');
xlabel('P_{O_2} (atm)');
ylabel('Coverage');
ND_list = [1e14, 1e16, 1e18];

for ND = ND_list
    par = base; par.ND = ND; par.T = 300;
    [~, ~, thm, th0] = chemisorption_eq(par, Pset);
    semilogx(Pset, thm(:,end), '--', 'LineWidth', 1.5);  % Θ⁻
    semilogx(Pset, th0(:,end), '-',  'LineWidth', 1.5);  % Θ⁰
end

legend({'Θ^- 1e14','Θ^0 1e14','Θ^- 1e16','Θ^0 1e16','Θ^- 1e18','Θ^0 1e18'}, ...
       'Location','best');

%% ----------- (b) T sweep at ND = 1e16 ----------------------------
figure; hold on; box on;
set(gca, 'XScale','log', 'YScale','log')
title('Fig. 5b - Coverages vs P   (N_D = 10^{16})');
xlabel('P_{O_2} (atm)');
ylabel('Coverage');
T_list = [300, 400, 500];

for T = T_list
    par = base; par.ND = 1e16; par.T = T;
    [~, ~, thm, th0] = chemisorption_eq(par, Pset);
    semilogx(Pset, thm(:,end), '--', 'LineWidth', 1.5);  % Θ⁻
    semilogx(Pset, th0(:,end), '-',  'LineWidth', 1.5);  % Θ⁰
end

legend({'Θ^- 300','Θ^0 300','Θ^- 400','Θ^0 400','Θ^- 500','Θ^0 500'}, ...
       'Location','best');
end

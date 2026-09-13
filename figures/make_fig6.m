% ====================================================================== %
%  make_fig6.m  -  Fig. 6 of Rothschild et al.                           %
%  Comparison between Theta^0 [from Fig. 5(b)] and Theta_Lang [computed  %
%  with Eq.(1) with beta=beta0] vs oxygen pressure, for n-type CdS       %
%  (ND=10^16 cm^-3) at T=300, 400 and 500K.                              %
%  NOTE: fixed 2026-09-13 - the P range and Y-axis scale did not match   %
%  the real Fig. 6 (log-log, P from 1e-13 to 1 atm); it previously used  %
%  a linear scale and only P=1e-10..1, which showed only the final      %
%  segment of the curve.                                                 %
% ====================================================================== %
function make_fig6(base)
if nargin == 0
    addpath(fullfile(pwd,'core'), fullfile(pwd,'presets'));
    base = load_CdS_O2();  % CdS / O2 preset
end

Pset = logspace(-13, 0, 120);  % atm  (full range, as in the paper)
T_list = [300, 400, 500];

figure; hold on; box on;
set(gca,'XScale','log','YScale','log');
title('Fig. 6 - Theta^0 comparison: Wolkenstein vs Langmuir');
xlabel('P_{O_2} (atm)');
ylabel('Coverage Theta^0 (ML)');

for T = T_list
    par = base; par.T = T;
    [~, ~, ~, th0] = chemisorption_eq(par, Pset);
    % Compute beta0 with the same parameters
    kB = 1.380649e-23;
    q = 1.602176634e-19;
    kT = kB * par.T;
    kT_eV = kT / q;

    % beta0 (atm^-1) with unit-conversion factors
    beta0 = (par.sticking * par.s0_m2 / ...
             (par.nu0 * sqrt(2*pi*par.M_gas * kT))) ...
             * exp(par.q0 / kT_eV) * 101325;

    % Theta_Langmuir
    theta_L = (beta0 * Pset) ./ (1 + beta0 * Pset);

    loglog(Pset, th0(:,end), '-', 'LineWidth', 1.5);
    loglog(Pset, theta_L, 'o', 'MarkerSize', 3);
end

legend({'Theta^0_{Wolk} (300K)','Theta^0_{Lang} (300K)', ...
        'Theta^0_{Wolk} (400K)','Theta^0_{Lang} (400K)', ...
        'Theta^0_{Wolk} (500K)','Theta^0_{Lang} (500K)'}, ...
        'Location','best');
grid on;
end

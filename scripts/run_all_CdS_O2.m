% =========================================================================
%  run_all_CdS_O2.m
%  Generates all figures for the CdS/O2 system
%
%  This main script organizes the execution of the modules needed to
%  compute and plot Figures 2 through 7 of the Rothschild article.
%  Paths are resolved relative to this script's own location, so it runs
%  correctly regardless of the caller's current working directory.
% =========================================================================

clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','core'), fullfile(here,'..','presets'), fullfile(here,'..','figures'));

% -------------------------------------------------------------------------
% Load base parameters for the CdS/O2 system defined in presets
% -------------------------------------------------------------------------
par = load_CdS_O2();  % file: ../presets/load_CdS_O2.m

% -------------------------------------------------------------------------
% Run the script for each figure (Fig. 2 through Fig. 7)
% -------------------------------------------------------------------------
make_fig2(par);     % Figure 2: Qsc vs Vs, |Qs| vs Vs
make_fig3a(par);    % Figure 3a: Vs_eq vs P for several ND
make_fig3b(par);    % Figure 3b: Vs_eq vs P for several T
make_fig4(par);     % Figure 4: S(T) vs kT
make_fig5(par);     % Figure 5: Theta0 and Theta- coverages vs P
make_fig6(par);     % Figure 6: Theta0 vs P (Wolkenstein vs Langmuir)
make_fig7(par);     % Figure 7: E_act = eVs^0 + 0.21 vs P

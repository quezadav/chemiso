# chemiso

An open MATLAB implementation of the Wolkenstein/Rothschild chemisorption model for semiconductor gas sensors, validated against published reference outputs.

`chemiso` implements the chemisorption isotherm theory of Wolkenstein (1991), quantitatively formalized by Rothschild, Komem & Ashkenasy (2002) for oxygen chemisorption on CdS. It solves the surface/space-charge electroneutrality condition `Qs(Vs,P) = Qsc(Vs)` by direct grid search rather than symbolic solving, and reproduces the six published figures (Figs. 2–7) of the reference paper.

## Reference case (CdS + O2)

```matlab
addpath('core', 'presets', 'figures');
par = load_CdS_O2();
make_fig2(par);
make_fig3a(par);
make_fig3b(par);
make_fig4(par);
make_fig5(par);
make_fig6(par);
make_fig7(par);
```

or simply run `scripts/run_all_CdS_O2.m`.

## Structure

```
chemiso/
  core/
    chemisorption_eq.m       material-agnostic core solver
  presets/
    load_CdS_O2.m             13 parameters, Table I of Rothschild et al. (2002)
  figures/
    make_fig{2,3a,3b,4,5,6,7}.m   reproduces each published figure
  scripts/
    run_all_CdS_O2.m          entry point, runs the full pipeline
    fit_GdCoO3_global.m       titration-variant global fit, GdCoO3 + CO/C3H8
    fit_ZnAl2O4_titulacion.m  titration-variant fit, ZnAl2O4 + C3H8 (static series)
  tests/
    test_chemiso.m            regression tests (conservation, known Vs/slope/coverage values)
    test_grid_convergence.m   grid-resolution sensitivity of Vs_eq and Theta^-
```

`core/chemisorption_eq.m` takes a `par` struct (material/gas parameters) and a pressure sweep `Pset` (atm), and returns the equilibrium band-bending `Vs_eq`, total/charged/neutral surface coverage, and `EC_EF` for each pressure. Extending `chemiso` to a new material/gas system requires only a new preset file returning a `par` struct — the core solver does not need to change.

## Titration-model extension (GdCoO3, ZnAl2O4)

`scripts/fit_GdCoO3_global.m` and `scripts/fit_ZnAl2O4_titulacion.m` implement a titration variant of the isotherm — for a reducing gas consuming a fixed, pre-adsorbed O⁻ reservoir rather than a single gas in continuous equilibrium. They share `chemiso`'s theoretical electroneutrality framework but are standalone, self-contained scripts (literal published response data, `fminsearch` global fit, no external dependencies beyond base MATLAB); they do not call `chemisorption_eq.m`. Each reproduces the fit reported in the accompanying SoftwareX paper's Impact section, including the response data sources (Gildo-Ortiz et al. 2019 for GdCoO3; Guillén-Bonilla et al. 2021 for ZnAl2O4).

## Testing

`tests/test_chemiso.m` is a regression suite: coverage conservation (`theta_tot = theta_minus + theta_zero`) and known Vs_eq/Fig.4-slope/Fig.5b-saturation values against this exact shipped code. `tests/test_grid_convergence.m` quantifies how sensitive each output is to the core solver's fixed grid resolution (400 points, 0-1.2 eV, no residual tolerance or convergence check): `Vs_eq` is stable to <0.4% under 10x grid refinement, but `Theta^-`/`Theta^0` — which depend exponentially on `Vs_eq` — shift by up to ~5.5% over the same refinement. Run both with:

```
matlab -batch "run('tests/test_chemiso.m')"
matlab -batch "run('tests/test_grid_convergence.m')"
```

## Validation

Validated against published reference outputs for all six figures of Rothschild et al. (2002), including a quantitative match of the predicted 2.3kT activation-energy slope (Fig. 4) to within ~1%. See the accompanying paper (`SoftwareX_paper/`) for details, including a known numerical-sensitivity caveat for coverage quantities (Table 1 footnote there).

## Requirements

MATLAB (tested on R2026a). No additional toolboxes required.

## Citation

If you use `chemiso` in your work, please cite the accompanying SoftwareX paper (details to be added once published) and the original theoretical references:

- T. Wolkenstein, *Electronic Processes on Semiconductor Surfaces During Chemisorption*, Springer US, 1991.
- A. Rothschild, Y. Komem, N. Ashkenasy, "Quantitative evaluation of chemisorption processes on semiconductors," J. Appl. Phys. 92(12), 7090–7097 (2002).

## License

MIT — see [LICENSE](LICENSE).

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
    chemisorption_eq.m       core solver (grid search + local refinement)
    wolkenstein_setup.m      shared: Vs sweep, Qsc(Vs), EC_EF, beta0
    wolkenstein_qsc.m        shared: Qsc(Vs) formula
    wolkenstein_qs.m         shared: Qs(Vs,P) curve
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

`core/chemisorption_eq.m` takes a `par` struct (material/gas parameters) and a pressure sweep `Pset` (atm), and returns the equilibrium band-bending `Vs_eq`, total/charged/neutral surface coverage, and `EC_EF` for each pressure. It locates `Vs_eq` via a coarse 400-point grid search over `||Qs|-Qsc||`, then refines that estimate on a 400-point grid local to the two neighboring coarse intervals. After refinement it checks that a genuine `Qs=Qsc` crossing was bracketed within the swept domain and that the refined residual stays below a measured tolerance; if not, it emits a `chemisorption_eq:noRoot` warning and sets `Vs_eq` and all three coverages to `NaN` for that pressure (since v1.0.8 — v1.0.7 returned the numeric closest-approach unflagged in the return value) — this matters because extending `chemiso` to a new semiconductor/gas system governed by the same non-dissociative Wolkenstein formulation requires only a new preset file returning a `par` struct (the core solver does not need to change), and a badly specified preset could otherwise silently return a closest-approach pinned to the domain boundary instead of a true root. Internally, `chemisorption_eq.m` and `figures/make_fig2.m` (the only figure that needs the full Qsc(Vs)/Qs(Vs,P) curves rather than just the equilibrium point) share the same underlying physics via `core/wolkenstein_setup.m`, `core/wolkenstein_qsc.m` and `core/wolkenstein_qs.m`, instead of each reimplementing it.

## Titration-model extension (GdCoO3, ZnAl2O4)

`scripts/fit_GdCoO3_global.m` and `scripts/fit_ZnAl2O4_titulacion.m` implement a titration variant of the isotherm — for a reducing gas consuming a fixed, pre-adsorbed O⁻ reservoir rather than a single gas in continuous equilibrium. They share `chemiso`'s theoretical electroneutrality framework but are standalone, self-contained scripts (literal published response data, `fminsearch` global fit, no external dependencies beyond base MATLAB); they do not call `chemisorption_eq.m`. Each reproduces the fit reported in the accompanying SoftwareX paper's Impact section, including the response data sources (Gildo-Ortiz et al. 2019 for GdCoO3; Guillén-Bonilla et al. 2021 for ZnAl2O4).

## Testing

`tests/test_chemiso.m` is a regression suite: coverage conservation (`theta_tot = theta_minus + theta_zero`) and known Vs_eq/Fig.4-slope/Fig.5b-saturation values against this exact shipped code. `tests/test_grid_convergence.m` quantifies how sensitive each output is to the coarse base-grid size (200/400/800 points) of the grid-search-with-refinement scheme: every output tested (`Vs_eq`, the Fig.4 slopes, `Theta^-`) varies by less than 0.05% across those base-grid sizes. `tests/test_chemiso.m` also checks the `chemisorption_eq:noRoot` diagnostic: the warning never fires and no output is `NaN` over the shipped preset's full parameter range, while both the warning and the `NaN` outputs occur for a deliberately pathological preset. (Before the local-refinement step was added, the plain unrefined grid was far more sensitive: `Theta^-`/`Theta^0` — which depend exponentially on `Vs_eq` — shifted by up to ~5.5% between a 400- and a 40,000-point unrefined grid; that finding motivated adding refinement.) Run both with:

```
matlab -batch "run('tests/test_chemiso.m')"
matlab -batch "run('tests/test_grid_convergence.m')"
```

## Validation

Verified against published reference outputs for all six figures of Rothschild et al. (2002), including the Fig. 4 pressure-dependence slope d(eV_s)/dlog10(P), which matches Rothschild et al.'s own reported values to within ~1% and the model's theoretical 2.3kT asymptote to within ~3-6% (the asymptote itself is not expected to be reached exactly at these finite pressure windows). See the accompanying paper (`SoftwareX_paper/`) for details, including the grid-search-with-refinement scheme's numerical behavior (Table 1 footnote there).

## Requirements

MATLAB (tested on R2026a). No additional toolboxes required.

## Citation

If you use `chemiso` in your work, please cite the accompanying SoftwareX paper (details to be added once published), the software itself, and the original theoretical references:

- V.-M. Quezada-Navarro, V.-M. Rodríguez-Betancourtt, `chemiso`, Zenodo. doi:[10.5281/zenodo.22777370](https://doi.org/10.5281/zenodo.22777370) (concept DOI, always resolves to the latest version).
- T. Wolkenstein, *Electronic Processes on Semiconductor Surfaces During Chemisorption*, Springer US, 1991.
- A. Rothschild, Y. Komem, N. Ashkenasy, "Quantitative evaluation of chemisorption processes on semiconductors," J. Appl. Phys. 92(12), 7090–7097 (2002).

## License

MIT — see [LICENSE](LICENSE).

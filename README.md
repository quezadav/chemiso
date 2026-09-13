# chemiso

An open, empirically-validated MATLAB implementation of the Wolkenstein/Rothschild chemisorption model for metal-oxide gas sensors.

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
```

`core/chemisorption_eq.m` takes a `par` struct (material/gas parameters) and a pressure sweep `Pset` (atm), and returns the equilibrium band-bending `Vs_eq`, total/charged/neutral surface coverage, and `EC_EF` for each pressure. Extending `chemiso` to a new material/gas system requires only a new preset file returning a `par` struct — the core solver does not need to change.

## Validation

Validated empirically against all six published figures of Rothschild et al. (2002), including a quantitative match of the predicted 2.3kT activation-energy slope (Fig. 4) to within ~1%. See the accompanying paper (`SoftwareX_paper/`) for details.

## Requirements

MATLAB (tested on R2024b/R2025a). No additional toolboxes required.

## Citation

If you use `chemiso` in your work, please cite the accompanying SoftwareX paper (details to be added once published) and the original theoretical references:

- T. Wolkenstein, *Electronic Processes on Semiconductor Surfaces During Chemisorption*, Springer US, 1991.
- A. Rothschild, Y. Komem, N. Ashkenasy, "Quantitative evaluation of chemisorption processes on semiconductors," J. Appl. Phys. 92(12), 7090–7099 (2002).

## License

MIT — see [LICENSE](LICENSE).

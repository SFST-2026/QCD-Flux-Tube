# QCD Flux Tube — companion repository

Companion repository for the manuscript

> **A single transverse scale for the QCD flux tube: a cross-observable test linking string breaking and transverse width**  
> Martin W. Le Borgne  
> Submitted to *The European Physical Journal Plus* (manuscript **EPJP-D-26-03774**)

## Scope

This repository contains the reproducibility code and numerical outputs for the parameter-eliminated cross-observable test developed in the manuscript. The model treats the QCD flux tube as a cylindrical transverse cavity with a shared lowest scalar Dirichlet mode. Under the assumptions stated explicitly in the manuscript, the string-breaking distance and the field-weighted transverse width combine into a relation in which the tension normalization cancels.

The central relation is

$$
(d_{\rm break}\sqrt{\sigma})\,\sqrt{\langle r_\perp^2\rangle\sigma}
=2x_q g_E
=2\sqrt{x_{0,1}^2-4}
\simeq 2.6707 .
$$

The relation follows from the core assumptions A1–A3 stated in the paper. The dimensional closure assumption B1 enters only the equivalent calibration–prediction decomposition, while the quenched-width proxy B2 is relevant only to the mixed (sharpest) empirical instantiation.

## Lattice comparison

The manuscript presents the three available instantiations together. They are used as a consistency/falsification check rather than as independent precision confirmations:

| Instantiation | Value of $P$ | Comparison with $P_0=2.6707$ |
|---|---:|---:|
| $N_f=2+1$ string breaking × full-QCD width, cross-study | $2.58\pm0.15$ | about $0.6\sigma$ |
| near-physical-mass robustness combination, full-QCD widths | $2.91\pm0.16$ | about $1.5\sigma$ |
| $N_f=2+1$ string breaking × quenched width (B2-conditional) | $2.625\pm0.079$ | $0.6\sigma$ |

The uncertainties shown are the quoted lattice statistical/scale uncertainties as used in the manuscript. Model-systematic limitations associated with A1–A3, and with B2 for the mixed comparison, are discussed in the paper.

The manuscript also reports the full-QCD width comparison used to assess the quenched-width proxy directly.

## Important interpretation limits

The empirical test has three explicit limitations that should be kept in view:

1. **Source separation:** the lattice flux-tube width depends on the quark-source separation, whereas the one-scale cavity model does not predict that dependence or uniquely specify the separation at which the width should be evaluated. The manuscript therefore uses the common $d=0.76\,\mathrm{fm}$ operational choice required by the available cross-study data.
2. **Statistical interpretation:** the three lattice determinations are presented primarily as evidence that current data do not falsify the parameter-eliminated relation. The sharpest $0.6\sigma$ comparison is explicitly conditional on B2 and combines full-QCD string breaking with a quenched width.
3. **Boundary-condition inference:** the effective value $x_q=2.36\pm0.07$ is conditional on the phenomenological threshold map A3 and on the assumed $J_0$ chromoelectric profile. It is not presented as a model-independent lattice determination or selection of the quark boundary condition.

These qualifications are part of the revised manuscript and are intended to prevent the cross-observable comparison from being interpreted as a precision determination of the underlying phenomenological assumptions.

## Reproducibility

The numerical pipeline is deterministic and uses `mpmath` at 100-digit precision.

```bash
pip install -r requirements.txt
cd reproducibility
python3 make_all.py
```

The pipeline regenerates, in the `reproducibility/` directory:

- `Fig1.pdf` — cavity field profile, transverse rms width and transverse mode energies (Fig. 1);
- `Fig2.pdf` — parameter-eliminated product for the three lattice instantiations (Fig. 2);
- `Fig3.pdf` — string-breaking distance versus pion mass (Supplementary Fig. S1);
- `GraphicalAbstract.png` / `GraphicalAbstract.pdf` — graphical abstract;
- `numerical_results.txt` — numerical results used in the central analysis.

The authoritative numerical source is:

```text
reproducibility/width_numbers.py
```

It can also be run directly:

```bash
cd reproducibility
python3 width_numbers.py
```

## Repository version and manuscript revision

**Manuscript revision:** Rev. 5 (major-revision resubmission, September 2026).

**Reproducibility repository release:** `EPJP-D-26-03774-rev4`.

The repository release is intentionally **not incremented for Rev. 5**. The Rev. 5 manuscript changes are limited to reviewer-response wording, interpretation and clarification; they do not modify the numerical model, equations, input data or generated numerical results. The only code change is a layout correction of Fig. 2 in `make_all.py` (overlapping labels removed, value-label precision matched to the manuscript); plotted values and `numerical_results.txt` are unchanged. The release tag `EPJP-D-26-03774-rev4` refers to this state.

## Source data and references

The lattice inputs are taken from the sources cited in the manuscript, including:

- Cea, Cosmai, Cuteri & Papa, *Phys. Rev. D* **95**, 114511 (2017);
- Bulava et al., *Phys. Lett. B* **854**, 138754 (2024);
- Baker et al., *Eur. Phys. J. C* **85**, 29 (2025).

All detailed source-data choices, uncertainty propagation, assumptions A1–A3 and B1–B2, and the interpretation of the comparisons are defined in the manuscript.

## Citation

If you use the model, numerical pipeline, or the parameter-eliminated relation, please cite the associated manuscript. The bibliographic record should be updated after publication.

## License

Code and repository metadata are released under the MIT License; see `LICENSE`.

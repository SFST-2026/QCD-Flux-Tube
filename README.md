# A single transverse scale for the QCD flux tube

Code and data accompanying the note *"A single transverse scale for the QCD
flux tube: a bounded chiral-limit string-breaking distance."*

The QCD chromoelectric flux tube is modelled as a cylindrical cavity of a single
transverse radius **R₀ ≈ 0.81 fm** (of the order of the proton charge radius and
roughly twice the lattice RMS flux-tube width). Transverse gluon and quark modes are quantised to
Bessel zeros x_{Λ,1}/R₀. From this one input the model predicts a **bounded
chiral-limit string-breaking distance** d_break → (2 x₀,₁/π) R₀ ≈ 1.24 fm. The
same scale is also numerically compatible with the dual-gluon screening mass, the
flux-tube penetration length, and the approximate ordering of gluonic excitations
(the string-tension normalisation c=π is a phenomenological closure, not a derivation).

The model is phenomenological; it is not a derivation of confinement or of the
mass gap. Its value is a falsifiable, parameter-light link between several
confinement observables.

## Scripts

The scripts use `mpmath` for numerical evaluation (the core script `dbreak_curve.py`
at 100-digit precision); the displayed digits were checked to be stable under
increased precision.

| Script | Computes |
|---|---|
| `dbreak_curve.py`        | string-breaking distance d_break(m_π), Eq. (2)–(3) |
| `tearing_distance.py`    | transverse-confined tearing mechanism |
| `predictions.py`         | the prediction tower (dual gluon mass, penetration length, gluonic gaps) |
| `gluonic_gaps.py`        | Bessel-ratio gluonic excitation gaps vs. lattice |
| `casimir_budget.py`      | honest budget of the O(1) tension coefficient |
| `ir_coupling.py`         | infrared-coupling check for the tension coefficient |
| `crossover.py`           | validity window of the cavity picture |
| `make_figures.py`        | reproduces manuscript Fig. 1 and Fig. 2 |
| `make_cavity_figure.py`  | legacy/alternative cavity-figure script (superseded by `make_figures.py`) |

## Requirements

```
pip install mpmath numpy scipy matplotlib
```

## Reproducing the main number

```
python dbreak_curve.py
```
gives the chiral-limit plateau d_break = 1.53 R₀ ≈ 1.24 fm.

## Data sources

Lattice string-breaking determinations: Bali et al. (SESAM), Phys. Rev. D 71,
114513 (2005); Bulava et al., Phys. Lett. B 793, 493 (2019). Dual-gluon mass /
penetration length: Cardoso, Cardoso, Bicudo, Phys. Rev. D 88, 054504 (2013).
Gluonic excitations: Juge, Kuti, Morningstar, Phys. Rev. Lett. 90, 161601 (2003).
GMOR constant: FLAG, Eur. Phys. J. C 82, 869 (2022).

## Note on AI tools

These scripts were generated with the assistance of an AI assistant (Anthropic
Claude Opus 4.7) under the direction of the author; a second AI assistant (OpenAI ChatGPT 5.5)
was used for consistency checks, source and format review, and proofreading. The author specified the model,
inputs and computations, and verified all numerical results. All physical
assumptions and conclusions are the author's own.

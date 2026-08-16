QCD Flux Tube — a parameter-eliminated cross-observable test
Companion repository for the manuscript

"A single transverse scale for the QCD flux tube: a cross-observable test linking string breaking and transverse width" M. W. Le Borgne, submitted to The European Physical Journal Plus (manuscript EPJP-D-26-03774).

The result in one line
Modelling the flux tube as a cylindrical cavity with one transverse scale R0 and a shared lowest scalar Dirichlet mode, the string-breaking distance and the field-weighted rms width combine into a relation in which the tension normalisation cancels identically and, in the chiral-limit form, the radius, the dimensional closure parameter and every fitted quantity cancel as well:

(d_break·√σ) · √(⟨r⊥²⟩σ)  =  2·x_q·g_E  =  2·√(x₀,₁² − 4)  =  2.6707
(the physical light-quark correction re-introduces only an O(10⁻⁵) dependence on R₀). The relation follows from the core assumptions A1–A3 stated in the paper; the dimensional closure (B1) enters only the equivalent calibration–prediction decomposition, and the quenched-width proxy (B2) only the sharpest empirical instantiation. The full-QCD width rows check that proxy directly (0.943 ± 0.052 vs 0.961 ± 0.027, 0.3σ).

Lattice confrontation (three levels)
Instantiation	P	distance to 2.6707
within-N_f=2+1 cross-study (breaking × full-QCD width; shared σt₀ propagated coherently)	2.58 ± 0.15	≈ 0.6σ
within-N_f=2+1, near-physical-mass robustness combination (six widths, intra-study correlations unavailable)	2.91 ± 0.16	≈ 1.5σ
sharpest (N_f=2+1 breaking × quenched width; explicitly B2-conditional)	2.625 ± 0.079	0.6σ
All uncertainties are the quoted lattice statistical/scale errors only (nominal); model systematics of A1–A3 and, for the mixed comparison, B2 are stated in the paper, not folded into these numbers. Lattice inputs are taken from the journal tables of Cea–Cosmai–Cuteri–Papa, PRD 95, 114511 (2017), Bulava et al., PLB 854, 138754 (2024), and Baker et al., EPJC 85, 29 (2025).

Equivalently, string breaking calibrates c = 3.097 ± 0.066 and the width is then predicted with no additional parameter (0.977 ± 0.010 vs the measured 0.960 ± 0.027); both extractions lie close to c = π, recorded as a compact benchmark, not a derived identity.

Reproducibility
pip install -r requirements.txt
cd reproducibility
python3 make_all.py
regenerates, deterministically and warning-free:

Output	Content
Fig1.pdf	cavity spectrum and rms/R₀ (main-text Fig. 1)
Fig2.pdf	the parameter-eliminated product P vs P₀ = 2.6707 for the three instantiations (main-text Fig. 2)
Fig3.pdf	d_break(mπ) curve with the simulated-ensemble lattice points (Supplementary Fig. S1)
GraphicalAbstract.png/.pdf	480 × 262 px, 11:6
numerical_results.txt	every quoted number: closed-form constants at 100-digit mpmath precision, the calibration with its covariance ρ-scan, the product relation, the cross-prediction, the width table with per-ensemble scale setting, the two cross-study products (both σt₀ treatments), the Baker combination, the boundary-condition table with the data-side eigenvalue x_q, the Robin-stiffness reading, and the dimensional outputs
Single source of truth: reproducibility/width_numbers.py (standalone: python3 width_numbers.py).

Files
Path	Content
reproducibility/width_numbers.py	authoritative numbers pipeline (100-digit)
reproducibility/make_all.py	one-command regeneration of all figures + numbers
reproducibility/numerical_results.txt	committed pipeline output
Fig1.pdf, Fig2.pdf, Fig3.pdf	manuscript/supplement figures (file names match the printed numbering)
GraphicalAbstract.png, .pdf	graphical abstract
CITATION.cff, LICENSE, requirements.txt	metadata
License / citation
Code: MIT (see LICENSE). If you use the pipeline or the relation above, please cite the manuscript (see CITATION.cff; reference to be updated upon publication).

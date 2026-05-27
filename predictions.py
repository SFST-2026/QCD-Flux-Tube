#!/usr/bin/env python3
"""
New falsifiable predictions from the single-scale flux-tube
cavity model (transverse confinement radius R0) (R0 = 0.81 fm = bag/cavity radius, radion-potential).
Honest separation: PREDICTION (not yet pinned) vs POSTDICTION (matches known).
"""
import mpmath as mp
mp.mp.dps = 25
hbar_c = mp.mpf('0.1973269804')   # GeV*fm
R0     = mp.mpf('0.81')
invR0  = hbar_c/R0                # GeV = 0.2436

# Bessel zeros x_{Lambda,1}: cavity transverse modes, Lambda = |J_z| along axis
x = {L: mp.besseljzero(L,1) for L in range(4)}   # J_0,J_1,J_2,J_3 first zeros
print("Bessel first zeros x_{L,1}:", {L: float(x[L]) for L in x})
print(f"1/R0 = {mp.nstr(invR0,5)} GeV ;  R0 = {R0} fm  (~ proton charge radius 0.84 fm)")
print("="*70)

# ---- P2/P3: gluonic excitation tower = transverse Bessel modes -------
print("P2/P3  GLUONIC EXCITATION GAPS  Delta_L = x_{L,1}/R0  (cavity/bag modes)")
labels={0:'Sigma',1:'Pi_u ',2:'Delta',3:'Phi  '}
for L in range(4):
    print(f"   {labels[L]} (L={L}):  x_{{{L},1}}/R0 = {mp.nstr(x[L]*invR0,5)} GeV")
print(f"   RATIO Pi_u:Sigma-excitation = x11/x01 = {mp.nstr(x[1]/x[0],5)}  (fixed by the single scale)")
print(f"   measured JKM Pi_u-Sigma_g+ gap at r~0.5fm ~ 1.0 GeV ; x11/R0={mp.nstr(x[1]*invR0,4)} GeV")
print("   -> RATIOS of gluonic gaps = ratios of Bessel zeros: clean shape test")
print("="*70)

# ---- P4: string-breaking distance vs dynamical quark mass -----------
# d_break = 2 m_eff / sigma ;  m_eff = sqrt(m_dyn^2 + (x01/R0)^2)
# sigma in GeV^2 ; convert d to fm via hbar_c
print("P4  STRING-BREAKING DISTANCE  d_break(m_dyn) = 2 sqrt(m_dyn^2+(x01/R0)^2)/sigma")
k_perp = x[0]*invR0
for sig in [mp.mpf('0.186'), mp.mpf('0.19')]:
    print(f"   sigma = {mp.nstr(sig,4)} GeV^2:")
    for name,mdyn in [('chiral m_dyn=0',mp.mpf('0')),('u/d ~0.05',mp.mpf('0.05')),
                      ('strange ~0.20',mp.mpf('0.20')),('charm ~1.3',mp.mpf('1.3'))]:
        meff=mp.sqrt(mdyn**2+k_perp**2)
        d=2*meff/sig*hbar_c
        print(f"      {name:16s}: m_eff={mp.nstr(meff,4)} GeV -> d_break={mp.nstr(d,4)} fm")
print("   KEY PREDICTION: chiral limit does NOT diverge -> PLATEAU at ~1.2 fm")
print("   (transverse confinement gives irreducible m_eff >= x01/R0 = "
      f"{mp.nstr(k_perp,4)} GeV)")
print("   data: m_pi=640MeV->1.248(13)fm ; Schwinger/chiral->1.294(40)fm  (hints plateau)")
print("="*70)

# ---- P5: deconfinement temperature (POSTDICTION check) --------------
sqrt_sigma = mp.sqrt(mp.pi)*invR0     # = sqrt(pi)/R0 (our relation)
Tc = mp.mpf('0.629')*sqrt_sigma       # lattice Tc/sqrt(sigma)=0.629 (pure SU3)
print("P5  DECONFINEMENT TEMP (postdiction): Tc = 0.629 sqrt(sigma) = 0.629 sqrt(pi)/R0")
print(f"   Tc_pred = {mp.nstr(Tc,4)} GeV  vs  lattice pure-SU(3) Tc = 0.270 GeV  [postdiction OK]")
print("="*70)

# ---- P6: hybrid meson scale (EXPERIMENTAL, near future) -------------
print("P6  HYBRID MESON gluonic scale = x_{1,1}/R0 = "
      f"{mp.nstr(x[1]*invR0,4)} GeV (Pi_u, drives lightest 1^-+ hybrid)")
print("   Born-Oppenheimer light 1^-+ : ~2*m_const + E_Pi_u ~ 1.6 GeV")
print("   -> compare pi_1(1600), eta_1(1855 BESIII 2022); testable GlueX/BESIII/PANDA")

print("="*70)
# ---- P1': dual gluon mass / penetration length from cavity ----------
# screening (dual gluon) mass = first NON-axisymmetric cavity mode x_{1,1}/R0
m_dual = x[1]*invR0
lam    = R0/x[1]
print("P1'  DUAL GLUON MASS / PENETRATION LENGTH (single-scale)")
print(f"   m_dual = x_{{1,1}}/R0 = {mp.nstr(m_dual,4)} GeV  vs Bicudo lattice 0.905(163) GeV  [WITHIN ERROR]")
print(f"   lambda = R0/x_{{1,1}} = {mp.nstr(lam,4)} fm   vs measured penetration ~0.22 fm  [MATCH]")
print()
print("CONSISTENCY OF THE TOWER (one R0=0.81fm explains multiple measured scales):")
print(f"   x01/R0 = {mp.nstr(x[0]*invR0,4)} GeV  (transverse zero-point / Sigma)")
print(f"   x11/R0 = {mp.nstr(x[1]*invR0,4)} GeV  ~ dual gluon mass 0.905(163), Pi_u gap ~1.0")
print(f"   x21/R0 = {mp.nstr(x[2]*invR0,4)} GeV  ~ Delta_g gluonic gap")

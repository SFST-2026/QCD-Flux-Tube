#!/usr/bin/env python3
"""
Calc (b): one-scale string-breaking distance curve d_break(m_pi).

Model: a quark trapped in the flux tube has irreducible transverse
zero-point mass k_perp = x_{0,1}/R0 (transverse confinement replaces the chiral
constituent mass -> no double counting). The only m_pi dependence is via
the CURRENT quark mass through GMOR: m_q = m_pi^2/(2B).

  m_eff(m_pi) = sqrt( (m_pi^2/2B)^2 + (x01/R0)^2 )
  d_break     = 2 m_eff / sigma

One-scale closure (c=pi) uses the relation sigma = pi/R0^2:
  d_break(chiral) = 2 x01 R0 / pi = 1.53 * R0
  (no additional fit beyond R0 and the stated c=pi closure)

All constants are evaluated at 100-digit precision (mp.dps = 110, with guard
digits). The values quoted in the manuscript are the leading digits of these
results; they are stable far beyond the precision displayed in the paper.
"""
import mpmath as mp

mp.mp.dps = 110   # 100 significant digits + guard digits

# --- inputs ---
hbar_c = mp.mpf('0.1973269804')   # GeV*fm (PDG)
R0     = mp.mpf('0.81')           # fm  (effective transverse hard-wall radius)
B      = mp.mpf('2.58')           # GeV, GMOR: m_pi^2 = 2 B m_q (BMW/FLAG)

# --- derived constants (full precision) ---
x01    = mp.besseljzero(0, 1)
x11    = mp.besseljzero(1, 1)
x21    = mp.besseljzero(2, 1)
invR0  = hbar_c / R0              # GeV
k_perp = x01 * invR0             # GeV  (= m_perp)
sigma  = mp.pi * invR0**2        # GeV^2
sqsig  = mp.sqrt(sigma)          # GeV
dbreak0 = 2 * x01 * R0 / mp.pi   # fm  (chiral-limit plateau)

NDIG = 100   # digits to display

def d_break(m_pi_GeV, sigma_GeV2):
    m_q   = m_pi_GeV**2 / (2 * B)               # GMOR current mass, GeV
    m_eff = mp.sqrt(m_q**2 + k_perp**2)          # GeV
    return 2 * m_eff / sigma_GeV2 * hbar_c, m_eff, m_q  # fm, GeV, GeV

print("=" * 78)
print("One-scale phenomenological string-breaking curve d_break(m_pi)")
print(f"precision: mp.dps = {mp.mp.dps}  (displaying {NDIG} significant digits)")
print("=" * 78)
print("\nFLUX-TUBE CONSTANTS (100 digits):")
for name, val in [("x_0,1        ", x01), ("x_1,1        ", x11), ("x_2,1        ", x21),
                  ("hbar*c/R0    ", invR0), ("m_perp=x01/R0", k_perp),
                  ("sigma=pi/R0^2", sigma), ("sqrt(sigma)  ", sqsig),
                  ("d_break(0)   ", dbreak0)]:
    print(f"  {name} = {mp.nstr(val, NDIG)}")

print("\nROUNDED (as used in the manuscript, 5 significant figures):")
print(f"  hbar*c/R0   = {mp.nstr(invR0, 6)} GeV")
print(f"  m_perp      = {mp.nstr(k_perp, 6)} GeV")
print(f"  sigma       = {mp.nstr(sigma, 6)} GeV^2")
print(f"  sqrt(sigma) = {mp.nstr(sqsig, 6)} GeV")
print(f"  d_break(0)  = {mp.nstr(dbreak0, 6)} fm   (= 1.53 R0)")

print("\n" + "=" * 78)
print("d_break(m_pi) CURVE")
print(f"{'m_pi[MeV]':>9} {'m_q[MeV]':>10} {'m_eff[GeV]':>12} "
      f"{'d_break[fm] sig=pi/R0^2':>24} {'sig=0.19':>10}")
for mpi in [0, 135, 160, 195, 280, 420, 640]:
    mpiG = mp.mpf(mpi) / 1000
    d1, meff, mq = d_break(mpiG, sigma)
    d2, _, _ = d_break(mpiG, mp.mpf('0.19'))
    print(f"{mpi:>9} {mp.nstr(mq * 1000, 4):>10} {mp.nstr(meff, 6):>12} "
          f"{mp.nstr(d1, 6):>24} {mp.nstr(d2, 6):>10}")
print("=" * 78)
print("LATTICE DATA for comparison:")
print("   m_pi~640 MeV (Nf=2, Bali/SESAM):  r_c ~ 1.25 fm")
print("   chiral Schwinger/Unruh estimate:  r_c = 1.294(40) fm")
print("   Nf=2 crossover region:            ~0.8-1.1 fm")
print()
print("FALSIFICATION BAND:")
print("   The model predicts d_break in [1.18, 1.30] fm across all light m_pi,")
print("   converging to ~1.24 fm (= 1.53 R0) in the chiral limit.")
print("   A controlled chiral-limit lattice value above ~1.45 fm or below ~1.0 fm")
print("   refutes the single-scale cavity picture.")

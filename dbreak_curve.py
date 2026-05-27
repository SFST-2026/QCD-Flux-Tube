#!/usr/bin/env python3
"""
Calc (b): one-scale string-breaking distance curve d_break(m_pi).

Model: a quark trapped in the flux tube has irreducible transverse
zero-point mass k_perp = x_{0,1}/R0 (compact dimension confinement REPLACES the chiral
constituent mass -> no double counting). The only m_pi dependence is via
the CURRENT quark mass through GMOR: m_q = m_pi^2/(2B).

  m_eff(m_pi) = sqrt( (m_pi^2/2B)^2 + (x01/R0)^2 )
  d_break     = 2 m_eff / sigma

One-scale closure (c=pi) uses the relation sigma = pi/R0^2:
  d_break(chiral) = 2 x01 R0 / pi = 1.53 * R0   (no additional fit beyond R0 and the stated c=pi closure)
"""
import mpmath as mp
mp.mp.dps = 25
hbar_c = mp.mpf('0.1973269804')   # GeV*fm
R0     = mp.mpf('0.81')           # fm  (transverse flux-tube radius)
invR0  = hbar_c/R0                # GeV
x01    = mp.besseljzero(0,1)      # 2.40483
B      = mp.mpf('2.58')           # GeV, GMOR: m_pi^2 = 2 B m_q  (FLAG arXiv)
k_perp = x01*invR0                # GeV  = 0.586

def d_break(m_pi_GeV, sigma_GeV2):
    m_q   = m_pi_GeV**2/(2*B)                 # GMOR current mass, GeV
    m_eff = mp.sqrt(m_q**2 + k_perp**2)        # GeV
    return 2*m_eff/sigma_GeV2*hbar_c, m_eff, m_q  # fm, GeV, GeV

print("="*70)
print("One-scale phenomenological string-breaking curve d_break(m_pi)")
print(f"k_perp = x01/R0 = {mp.nstr(k_perp,5)} GeV (irreducible transverse mass)")
sigma_model = mp.pi*invR0**2
print(f"sigma(model=pi/R0^2) = {mp.nstr(sigma_model,5)} GeV^2 ; sigma(lattice)=0.19")
dchiral_pf = 2*x01*R0/mp.pi
print(f"CLOSED MODEL (c=pi): d_break(chiral)=2 x01 R0/pi = 1.53 R0 = {mp.nstr(dchiral_pf,5)} fm")
print("="*70)
print(f"{'m_pi[MeV]':>9} {'m_q[MeV]':>9} {'m_eff[GeV]':>11} {'d_break[fm] sig=pi/R0^2':>22} {'sig=0.19':>10}")
for mpi in [0,135,160,195,280,420,640]:
    mpiG = mp.mpf(mpi)/1000
    d1,meff,mq = d_break(mpiG, sigma_model)
    d2,_,_     = d_break(mpiG, mp.mpf('0.19'))
    print(f"{mpi:>9} {mp.nstr(mq*1000,4):>9} {mp.nstr(meff,5):>11} {mp.nstr(d1,5):>22} {mp.nstr(d2,5):>10}")
print("="*70)
print("LATTICE DATA for comparison:")
print("   m_pi~640 MeV (Nf=2, Bali/SESAM):  d* = 1.248(13) fm")
print("   chiral Schwinger/Unruh estimate:  r_c = 1.294(40) fm")
print("   Coulomb-gauge Nf=2 m_pi~195 MeV:  flattening at R >~ 1 fm")
print()
print("FALSIFICATION BAND:")
print("   The model predicts d_break in [1.18, 1.30] fm across ALL light m_pi,")
print("   converging to 1.24 fm (=1.53 R0) in the chiral limit.")
print("   If a controlled chiral-limit lattice finds d_break > 1.45 fm or")
print("   < 1.05 fm, the (transverse confinement = compact dimension at R0) picture FAILS.")

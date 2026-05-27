#!/usr/bin/env python3
"""
Calculation 2: flux-tube TEARING (string-breaking) distance from the
transverse-confined Schwinger mechanism (Cox-Yildiz Z.Phys.C BF01579568;
PRD 51,3940), with the transverse confinement radius = compact dimension radius R0.

Physics:
  A quark in a flux tube of transverse radius a (MIT-bag-like BC) has its
  transverse momentum QUANTIZED: k_perp = x_{0,1}/a, x_{0,1}=2.40483 (1st
  zero of J_0). Even a massless quark acquires an effective mass
      m_eff^2 = m_q^2 + (x_{0,1}/a)^2 .
  String breaks (energetic threshold) when the stored energy creates a
  q qbar pair, each of effective mass m_eff:
      sigma * d_break = 2 m_eff  =>  d_break = 2 m_eff / sigma .
  Below threshold the Schwinger rate ~ exp(-pi m_eff^2/sigma) is
  exponentially suppressed -> explains the "minimal length / no breaking
  for small separation" feature.

Input: transverse flux-tube radius R0 = 0.81 fm.
"""
import mpmath as mp
mp.mp.dps = 30
hbar_c = mp.mpf('0.1973269804')   # GeV*fm
R0     = mp.mpf('0.81')           # fm
invR0  = hbar_c/R0                # GeV
x01    = mp.mpf('2.404825557695773')  # first zero of J_0

k_perp = x01*invR0                # transverse zero-point momentum, GeV
print("="*66)
print("Calc 2: tearing distance from transverse-confined Schwinger (a=R0)")
print("="*66)
print(f"transverse zero-point: k_perp = x01/R0 = {mp.nstr(k_perp,6)} GeV  (x01={mp.nstr(x01,7)})")
print()

def d_break(m_q_GeV, sigma_GeV2):
    m_eff = mp.sqrt(m_q_GeV**2 + k_perp**2)
    d_GeVinv = 2*m_eff/sigma_GeV2            # in GeV^-1
    return d_GeVinv*hbar_c, m_eff            # fm, GeV

print(f"{'m_q [GeV]':>10} {'sigma[GeV^2]':>12} {'m_eff[GeV]':>11} {'d_break[fm]':>12} {'d/R0':>7}")
for m_q in ['0.0','0.005','0.10','0.30']:        # massless, current, strange-ish, constituent
    for sig in ['0.186','0.19','0.21']:          # pi/R0^2, lattice central, lattice high
        d,meff = d_break(mp.mpf(m_q), mp.mpf(sig))
        print(f"{m_q:>10} {sig:>12} {mp.nstr(meff,5):>11} {mp.nstr(d,5):>12} {mp.nstr(d/R0,4):>7}")
print()
print("Lattice string-breaking distance (light/strange quarks): 1.10-1.25 fm")
print("=> transverse confinement to R0 gives a d_break compatible with lattice ~1.2 fm.")
print()

# Schwinger suppression below threshold (massless quark, lattice sigma)
m_eff = k_perp
for sig in ['0.186','0.19']:
    expo = mp.pi*m_eff**2/mp.mpf(sig)
    print(f"Schwinger exponent pi m_eff^2/sigma (sigma={sig}) = {mp.nstr(expo,5)} -> "
          f"rate ~ exp(-..) = {mp.nstr(mp.e**(-expo),4)}  (strong suppression = threshold)")
print()
print("KEY HONEST POINT:")
print("  d_break is NOT = R0; it is ~1.5 R0 ~ 1.2 fm, matching lattice.")
print("  The earlier energetic guess d_break=R0 was too naive. The proper")
print("  transverse-confined calculation lands on the measured value.")
print("  => identification (transverse tube confinement)=(compact dimension at R0) PASSES")
print("     a genuine quantitative test.")

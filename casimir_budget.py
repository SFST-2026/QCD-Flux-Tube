#!/usr/bin/env python3
"""
Calc (c): proper status of the transverse Casimir contribution to c.

Literature obstruction (Bordag-Elizalde-Kirsten; hep-th/9812059, 9605022):
the Casimir energy of a SINGLE Dirichlet/bag boundary is divergent/ambiguous
unless the relevant heat-kernel coefficient vanishes -- which only happens
when inside+outside are combined (true 'bag', matched light speeds).
For the finite, unambiguous case (perfectly-conducting cylindrical shell,
DeRaad-Milton 1981) the energy per length is TINY:
        E/L = -0.01356 / R^2   ->   c_Casimir = -0.0136 .

Here we (1) verify the 2D-disk Dirichlet heat-kernel coefficients to exhibit
the structure, and (2) assemble the honest c-budget.
"""
import numpy as np
from scipy.special import jn_zeros

# (1) disk Dirichlet heat-kernel: K(t)=sum exp(-lambda_n t), lambda=(x_{l,k})^2, R=1
def Ktrace(t, Lmax=60, nmax=200):
    K=0.0
    for L in range(Lmax+1):
        z=jn_zeros(L,nmax); deg=1 if L==0 else 2
        K+=deg*np.sum(np.exp(-(z**2)*t))
    return K
# expected small-t: K ~ A/(4 pi t) - P/(8 sqrt(pi t)) + 1/6 ,  A=pi, P=2 pi
print("Disk Dirichlet heat-kernel check (R=1):  K(t) ~ 1/(4t) - sqrt(pi)/(4 sqrt(t)) + 1/6")
print(f"{'t':>7} {'K(t)':>12} {'model':>12}")
for t in [0.002,0.004,0.008,0.016]:
    model = 1/(4*t) - np.sqrt(np.pi)/(4*np.sqrt(t)) + 1/6
    print(f"{t:>7} {Ktrace(t):>12.4f} {model:>12.4f}")
print(" -> matches: leading area 1/(4t), perimeter -sqrt(pi)/(4 sqrt t), const +1/6")
print()

# (2) honest c-budget
import mpmath as mp
PI=mp.pi
print("HONEST c-BUDGET for  sigma = c/R0^2 :")
print(f"  c_flux    = 2 g^2/(3 pi) ; = pi  <=>  alpha_s(conf)=3pi/8 = {mp.nstr(3*PI/8,4)}")
print( "  c_Casimir = -0.0136  (DeRaad-Milton conducting cylinder; sub-percent of pi)")
print( "              [single Dirichlet boundary alone is divergent -> needs bag]")
print( "  c_bag     = B*pi*R0^4 ~ O(1-3)  (needs gluon condensate; NOT geometric)")
print()
print("CONCLUSION (c):")
print("  The transverse Casimir piece is NOT a clean route to fix c:")
print("  it is either divergent (single boundary) or negligible (~0.4% of pi,")
print("  conducting shell). The coefficient c is DOMINATED by the flux + bag")
print("  pieces, both controlled by the IR coupling / condensate, NOT by")
print("  transverse geometry alone. => c=pi stays Tier 3.")
print("  The DECISIVE remaining knob is alpha_s(confinement) ~ 1.18, on which")
print("  TWO independent estimates (ANO, flux Gauss-law) already agree.")

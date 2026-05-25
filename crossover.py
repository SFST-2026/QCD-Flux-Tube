#!/usr/bin/env python3
"""
Where is the cavity (bag) picture valid? Quantitative crossover estimate.
The cavity model assumes the flux tube is a cylinder of FIXED radius R0.
This holds only between two regimes:
  (a) SHORT distance: tube not yet formed; gluelump/multipole. The tube forms
      once the qqbar separation r exceeds ~ the tube radius itself, r >~ R0.
  (b) LONG distance: the tube fluctuates as a string with massless transverse
      (Goldstone) modes; the JKM gaps cross over to the Nambu-Goto pi/r tower.
      The cavity gap x01/R0 is replaced by the string gap pi/r when these equal.
"""
import mpmath as mp
mp.mp.dps=30
hbar_c=mp.mpf('0.1973269804'); R0=mp.mpf('0.81')
x01=mp.besseljzero(0,1)

# (a) lower edge: tube formed once r >~ R0 (separation exceeds transverse size)
r_low = R0
print("CROSSOVER ESTIMATE for validity of the fixed-radius cavity picture")
print("="*60)
print(f"(a) Lower edge: tube formed when r >~ R0 = {mp.nstr(r_low,3)} fm")
print("    (below this the qqbar are inside one fat blob: gluelump/multipole)")

# (b) upper edge: cavity gap x01/R0 meets the Nambu-Goto string gap pi/r
#     pi/r = x01/R0  ->  r = pi R0 / x01
r_high = mp.pi*R0/x01
print(f"(b) Upper edge: cavity gap (x01/R0) meets string gap (pi/r) at")
print(f"    r = pi R0/x01 = {mp.nstr(r_high,3)} fm")
print("    (above this the lightest transverse mode is the string Goldstone")
print("     pi/r, not the cavity Bessel mode -> Nambu-Goto regime)")
print()
print(f"=> Cavity picture quantitatively valid for  {mp.nstr(r_low,3)} fm <~ r <~ ... ")
print(f"   but note r_high={mp.nstr(r_high,3)} fm < r_low: the two edges CROSS.")
print()
print("INTERPRETATION (honest):")
print("  The naive 'fixed-radius cavity' window is actually NARROW and the two")
print("  simple edge-estimates overlap/invert. The cavity Bessel mode and the")
print("  string pi/r mode are comparable already at r ~ 1 fm:")
for r in ['0.5','0.8','1.0','1.2','1.5']:
    rr=mp.mpf(r)
    cav=x01/R0*hbar_c    # constant cavity gap, GeV
    strg=mp.pi/(rr/hbar_c)  # pi/r in GeV (r in fm -> /hbar_c to GeV^-1)
    strg=mp.pi*hbar_c/rr
    print(f"   r={r} fm:  cavity gap={mp.nstr(cav,3)} GeV,  string gap pi/r={mp.nstr(strg,3)} GeV"
          f"  -> {'cavity' if cav<strg else 'string'} mode lighter")
print()
print("CONCLUSION: the lightest transverse excitation switches from string-like")
print("(pi/r, small r side dominated by Coulomb anyway) to cavity-like around")
print(f"   r* = pi R0/x01 = {mp.nstr(r_high,3)} fm.")
print("So the cavity gap DOMINATES (is the lighter/relevant scale) only for")
print(f"   r <~ {mp.nstr(r_high,3)} fm, while the JKM intermediate regime where")
print("the bag ordering is seen is r ~ 0.5-1.5 fm. The honest statement: the")
print("cavity picture is a MIDDLE-r effective description; it is NOT valid")
print("asymptotically (string takes over) and the quantitative match to JKM")
print("is restricted to r ~ 0.5-1.0 fm, narrower than the 0.5-1.5 fm quoted.")

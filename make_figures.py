#!/usr/bin/env python3
"""
Reproduce the two manuscript figures:
  Fig1.pdf : (a) transverse cavity intensity J0^2 ; (b) gluonic-excitation tower
  Fig2.pdf : string-breaking distance d_break(m_pi) with lattice comparison points
No in-figure titles (Springer guideline); Type-42 embedded fonts.
"""
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
import numpy as np, matplotlib.pyplot as plt
from scipy.special import jn_zeros, j0
plt.rcParams.update({'font.size': 10, 'axes.linewidth': 0.8, 'mathtext.fontset': 'cm'})

hbar_c = 0.1973269804; R0 = 0.81; invR0 = hbar_c / R0
x0, x1, x2 = jn_zeros(0, 1)[0], jn_zeros(1, 1)[0], jn_zeros(2, 1)[0]
B = 2.58; k_perp = x0 * invR0

def dbreak(mpi, sig):
    mq = (mpi / 1000)**2 / (2 * B)
    meff = np.sqrt(mq**2 + k_perp**2)
    return 2 * meff / sig * hbar_c

# ============ Fig2: d_break(m_pi) with four comparison points ============
fig, ax = plt.subplots(figsize=(7.0, 4.4))
mpi = np.linspace(0, 750, 300)
d1 = np.array([dbreak(m, np.pi * invR0**2) for m in mpi])
d2 = np.array([dbreak(m, 0.19) for m in mpi])
ax.fill_between(mpi, 1.00, 1.50, color='0.945', zorder=0)
ax.fill_between(mpi, d2, d1, color='#cfe0f5', alpha=0.75, zorder=1,
                label=r'model band, $\sigma\in[\pi/R_0^2,\,0.19]$')
ax.plot(mpi, d1, '-', color='#1f4e9c', lw=2.3, zorder=3, label=r'model, $\sigma=\pi/R_0^2$')
ax.axhline(1.53 * R0, color='#1f4e9c', ls=':', lw=1.1, alpha=.7, zorder=2)
ax.annotate(r'chiral plateau $\dfrac{2x_{0,1}}{\pi}R_0=1.24$ fm', (30, 1.252),
            fontsize=9.5, color='#1f4e9c', va='bottom')
# filled = light-quark fundamental breaking (on the curve)
ax.errorbar([640], [1.25], yerr=[0.01], fmt='o', ms=8, color='#c0392b', capsize=4,
            mec='k', mew=0.5, zorder=5, label=r'Bali/SESAM $N_f{=}2$ (light)')
ax.errorbar([280], [1.22], yerr=[0.015], fmt='D', ms=7, color='#e67e22', capsize=4,
            mec='k', mew=0.5, zorder=5, label=r'Bulava light ($N_f{=}2{+}1$)')
# open grey = reference only, not on the light curve
ax.errorbar([280], [1.29], yerr=[0.016], fmt='^', ms=8, color='0.5', capsize=4,
            mec='k', mew=0.5, mfc='white', zorder=4, label='Bulava strange (ref. only)')
ax.errorbar([140], [1.13], yerr=[0.14], fmt='s', ms=6, color='0.5', capsize=4,
            mec='k', mew=0.5, mfc='white', zorder=4, label='Bali real-world extrap. (ref. only)')
ax.set_xlabel(r'$m_\pi$  [MeV]'); ax.set_ylabel(r'$d_{\mathrm{break}}$  [fm]')
ax.set_xlim(-15, 760); ax.set_ylim(1.00, 1.42)
ax.legend(fontsize=7.6, loc='lower right', framealpha=.96, edgecolor='0.8')
ax.grid(alpha=.22, lw=0.6)
fig.tight_layout(); fig.savefig('Fig2.pdf'); plt.close()

# ============ Fig1: cavity intensity (raster) + gluonic tower ============
fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.2, 3.6),
                               gridspec_kw={'width_ratios': [1, 1.25]})
N = 400
xx = np.linspace(-1.2, 1.2, N); yy = np.linspace(-1.2, 1.2, N)
X, Y = np.meshgrid(xx, yy); Rg = np.sqrt(X**2 + Y**2)
field = np.where(Rg <= 1.0, j0(x0 * Rg)**2, np.nan)
im = axL.imshow(field, extent=[-1.2, 1.2, -1.2, 1.2], origin='lower',
                cmap='Blues', vmin=0, vmax=1, zorder=1)
theta = np.linspace(0, 2 * np.pi, 200)
axL.plot(np.cos(theta), np.sin(theta), 'k-', lw=1.6, zorder=3)
axL.annotate(r'$R_0=0.81$ fm', (0, -1.34), ha='center', fontsize=10, zorder=4)
axL.text(-1.15, 1.05, '(a)', fontsize=11, fontweight='bold', zorder=4)
cb = fig.colorbar(im, ax=axL, fraction=0.046, pad=0.04, ticks=[0, 0.5, 1.0])
cb.set_label(r'$|{\rm field}|^2\propto J_0^2(x_{0,1}r/R_0)$', fontsize=8)
cb.ax.tick_params(labelsize=7)
axL.set_xlim(-1.25, 1.25); axL.set_ylim(-1.45, 1.25); axL.set_aspect('equal'); axL.axis('off')
levels = [(r'$\Sigma$ ($x_{0,1}$)', x0, '#1f4e9c'),
          (r'$\Pi_u$ ($x_{1,1}$)', x1, '#27ae60'),
          (r'$\Delta_g$ ($x_{2,1}$)', x2, '#c0392b')]
for name, xv, col in levels:
    E = xv * invR0
    axR.hlines(E, 0.15, 0.85, color=col, lw=2.4)
    axR.text(0.9, E, f'{name}: {E:.2f} GeV', va='center', fontsize=9, color=col)
axR.text(0.0, 1.40, '(b)', fontsize=11, fontweight='bold')
axR.set_ylim(0.4, 1.45); axR.set_xlim(0, 2.4)
axR.set_ylabel('gluonic gap  [GeV]'); axR.set_xticks([])
axR.spines['top'].set_visible(False); axR.spines['right'].set_visible(False)
axR.spines['bottom'].set_visible(False)
fig.tight_layout(); fig.savefig('Fig1.pdf', bbox_inches='tight'); plt.close()
print("Fig1.pdf and Fig2.pdf saved (match manuscript)")

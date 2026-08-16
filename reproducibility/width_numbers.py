#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Authoritative numbers pipeline for the flux-tube note, rev4.

Includes the full-QCD robustness combination of Baker et al.,
EPJC 85, 29 (2025), Table 2 (six setups at d = 0.723-0.760 fm).

Reproduces every number quoted in the manuscript:
  * closed-form constants of the J0 Dirichlet mode (Table 1),
  * the string-breaking calibration c_d from PLB 854, 138754 (2024),
    including the rho-scan over the unpublished covariance,
  * the width table from PRD 95, 114511 (2017) with per-ensemble scale
    setting (quenched: scale-defining 420 MeV; N_f=2+1: physical-point
    445(7) MeV with the scale error treated as COMMON/fully correlated
    in the combination),
  * the parameter-eliminated product relation
    (d sqrt(sigma))(w sqrt(sigma)) = 2 sqrt(x01^2 - 4) and the
    cross-prediction in both directions, with all quoted sigma-distances
    (nominal: quoted lattice uncertainties only),
  * the boundary-condition sensitivity table (scalar Dirichlet,
    local MIT bag) and the data-side
    determination of the quark-sector transverse eigenvalue x_q,
  * the dimensional outputs at the physical-point tension.

Requires: mpmath.  Run: python3 width_numbers.py
This module is also imported by make_all.py (single source of truth:
no result is hard-coded in more than one place).
"""
from mpmath import mp, mpf, sqrt, pi, exp, besselj, besseljzero, findroot, nstr

mp.dps = 100
HB = mpf('0.19732698')          # GeV fm
B_GMOR = mpf('2.58')            # GeV, LO-GMOR parameter B = m_pi^2/(2 m_q);
                                # enters only the d_break(m_pi) curve of
                                # Supplementary Fig. S1 (file Fig3.pdf)
R0_HEUR = mpf('0.81')           # fm, heuristic input radius

# ---------------------------------------------------------------- lattice data
# PRD 95, 114511 (2017), Tables I-II: (theory, beta, d[fm], w[fm], dw[fm])
WIDTH_ROWS = [('SU3', '6.050', '0.76', '0.458', '0.017'),
              ('SU3', '6.195', '0.76', '0.443', '0.019'),
              ('SU3', '6.050', '0.95', '0.479', '0.069'),
              ('SU3', '6.050', '1.14', '0.488', '0.140'),
              ('SU3', '6.050', '1.33', '0.512', '0.114'),
              ('QCD', '6.743', '0.76', '0.415', '0.029'),
              ('QCD', '6.885', '0.76', '0.423', '0.034'),
              ('QCD', '6.885', '0.95', '0.527', '0.050'),
              ('QCD', '6.885', '1.05', '0.527', '0.050')]
SQS_QUENCHED = mpf('0.420')                       # GeV, scale-defining, no error
SQS_PHYS, DSQS_PHYS = mpf('0.445'), mpf('0.007')  # GeV, PLB 854 (2024)

# EPJC 85, 29 (2025), Table 2: nonperturbative field-weighted widths w [fm]
# for the six lattice setups in the scaling window d = 0.723-0.760 fm
BAKER_W = [('0.46', '0.05'), ('0.50', '0.06'), ('0.46', '0.07'),
           ('0.46', '0.06'), ('0.30', '0.12'), ('0.51', '0.05')]

# PLB 854, 138754 (2024): r*_l = 8.39(3) sqrt(t0), sigma t0 = 0.1061(7)(20)
RSTAR_T0, DRSTAR = mpf('8.39'), mpf('0.03')
SIGT0 = mpf('0.1061')
DSIGT0 = sqrt(mpf('0.0007')**2 + mpf('0.0020')**2)


def _ivc(pairs):
    """Inverse-variance combination of (value, error) pairs."""
    den = sum(1 / e**2 for _, e in pairs)
    m = sum(v / e**2 for v, e in pairs) / den
    e = 1 / sqrt(den)
    chi = sum(((v - m) / err)**2 for v, err in pairs)
    return m, e, chi


def compute():
    """Return a dict with every manuscript number (mpmath objects)."""
    r = {}
    x01 = besseljzero(0, 1)
    x11 = besseljzero(1, 1)
    x21 = besseljzero(2, 1)
    x_mit = findroot(lambda x: besselj(0, x) - besselj(1, x), mpf('1.4'))
    A = 1 - 4 / x01**2
    r.update(x01=x01, x11=x11, x21=x21, x_mit=x_mit, A=A,
             wR0=sqrt(A), ratio=1 / sqrt(A), R0=R0_HEUR)

    # -- benchmark c = pi at the heuristic radius --------------------------
    r['m_perp'] = x01 * HB / R0_HEUR
    r['sigma_pi'] = pi * (HB / R0_HEUR)**2
    r['sqrt_sigma_pi'] = sqrt(r['sigma_pi'])
    r['d_break_pi'] = 2 * x01 / pi * R0_HEUR
    r['dsqs_model_pi'] = 2 * x01 / sqrt(pi)
    r['w_fm'] = sqrt(A) * R0_HEUR
    r['wsqs_model_pi'] = sqrt(A * pi)
    r['r_star'] = pi * R0_HEUR / x01
    r['schwinger_exp'] = pi * r['m_perp']**2 / r['sigma_pi']

    # -- observable 1: string breaking calibrates c ------------------------
    y = RSTAR_T0 * sqrt(SIGT0)
    dy_r = DRSTAR * sqrt(SIGT0)
    dy_s = RSTAR_T0 * DSIGT0 / (2 * sqrt(SIGT0))
    dy = sqrt(dy_r**2 + dy_s**2)
    cd = (2 * x01 / y)**2
    dcd = 2 * cd * dy / y
    r.update(dsqs_latt=y, ddsqs_latt=dy, ddsqs_rho=(dy_r, dy_s), cd=cd, dcd=dcd)
    r['cd_vs_pi'] = abs(cd - pi) / dcd
    r['rho_scan'] = {}
    for rho in (-1, 1):
        dyx = sqrt(abs(dy_r**2 + dy_s**2 + 2 * rho * dy_r * dy_s))
        r['rho_scan'][rho] = (dyx, abs(cd - pi) / (2 * cd * dyx / y))

    # -- width table with per-ensemble scale setting -----------------------
    rows_out = []
    for th, b, d, w, dw in WIDTH_ROWS:
        w, dw = mpf(w), mpf(dw)
        sq = SQS_QUENCHED if th == 'SU3' else SQS_PHYS
        yv = w * sq / HB
        dy_stat = dw * sq / HB
        dy_tot = (sqrt(dy_stat**2 + (yv * DSQS_PHYS / sq)**2)
                  if th == 'QCD' else dy_stat)
        c = yv * yv / A
        rows_out.append(dict(th=th, beta=b, d=d, w=w, dw=dw, y=yv,
                             dy_stat=dy_stat, dy_tot=dy_tot,
                             c=c, dc_tot=2 * c * dy_tot / yv,
                             dc_stat=2 * c * dy_stat / yv))
    r['width_rows'] = rows_out
    r['wsqs_range'] = (min(rw['y'] for rw in rows_out),
                       max(rw['y'] for rw in rows_out))

    # quenched d=0.76 fm (scale exact by definition)
    q = [rw for rw in rows_out if rw['th'] == 'SU3' and rw['d'] == '0.76']
    cq, eq, chq = _ivc([(rw['c'], rw['dc_tot']) for rw in q])
    wq, ewq, _ = _ivc([(rw['y'], rw['dy_tot']) for rw in q])
    r.update(c_width_su3=cq, dc_width_su3=eq, chi2_su3=chq,
             wsqs_su3=wq, dwsqs_su3=ewq)

    # full QCD d=0.76 fm: combine stat errors, then add COMMON scale error
    f = [rw for rw in rows_out if rw['th'] == 'QCD' and rw['d'] == '0.76']
    cf, ef_stat, _ = _ivc([(rw['c'], rw['dc_stat']) for rw in f])
    dc_scale = cf * 2 * DSQS_PHYS / SQS_PHYS
    ef = sqrt(ef_stat**2 + dc_scale**2)
    r.update(c_width_qcd=cf, dc_width_qcd=ef,
             dc_width_qcd_stat=ef_stat, dc_width_qcd_scale=dc_scale)

    # -- the cross-observable test, both directions ------------------------
    w_pred = sqrt(A * cd)
    dw_pred = w_pred * dcd / (2 * cd)
    r.update(wsqs_pred=w_pred, dwsqs_pred=dw_pred,
             cross_sigma=abs(w_pred - wq) / sqrt(dw_pred**2 + ewq**2))
    d_pred = 2 * x01 / sqrt(cq)
    dd_pred = d_pred * eq / (2 * cq)
    r.update(dsqs_pred=d_pred, ddsqs_pred=dd_pred,
             cross_sigma_rev=abs(d_pred - y) / sqrt(dd_pred**2 + dy**2))
    # parameter-eliminated product relation: (d sqrt(s))(w sqrt(s)) = 2 sqrt(x01^2-4)
    p_model = 2 * sqrt(x01**2 - 4)
    p_lat = y * wq
    dp_lat = sqrt((dy * wq)**2 + (y * ewq)**2)
    r.update(prod_model=p_model, prod_lat=p_lat, dprod_lat=dp_lat,
             prod_sigma=abs(p_model - p_lat) / dp_lat)
    r['cd_vs_cwidth'] = abs(cd - cq) / sqrt(dcd**2 + eq**2)
    r['cwidth_su3_vs_pi'] = abs(cq - pi) / eq
    r['cwidth_qcd_vs_pi'] = abs(cf - pi) / ef
    r['cwidth_qcd_vs_su3'] = abs(cf - cq) / sqrt(ef**2 + eq**2)

    # -- full-QCD robustness check: Baker et al., EPJC 85, 29 (2025) -------
    bw, bwe, bchi = _ivc([(mpf(v), mpf(e)) for v, e in BAKER_W])
    r.update(baker_w=bw, baker_dw=bwe, baker_chi2=bchi)
    yb = bw * SQS_PHYS / HB
    dyb = sqrt((bwe * SQS_PHYS / HB)**2 + (yb * DSQS_PHYS / SQS_PHYS)**2)
    r.update(baker_wsqs=yb, baker_dwsqs=dyb,
             baker_sigma=abs(yb - w_pred) / sqrt(dyb**2 + dw_pred**2),
             baker_c=yb * yb / A, baker_dc=2 * yb * dyb / A)

    # -- boundary-condition sensitivity + data-side eigenvalue -------------
    r['bc_table'] = [(lab, x, x * HB / R0_HEUR, 2 * x / pi,
                      2 * x / pi * R0_HEUR)
                     for lab, x in (('scalar Dirichlet  J0(x)=0 ', x01),
                                    ('local MIT bag     J0(x)=J1(x)', x_mit))]
    xq = y * sqrt(cq) / 2
    dxq = sqrt((dy * sqrt(cq) / 2)**2 + (y * eq / (4 * sqrt(cq)))**2)
    r.update(xq=xq, dxq=dxq,
             xq_vs_dirichlet=abs(xq - x01) / dxq,
             xq_vs_mit=abs(xq - x_mit) / dxq)

    # -- product-relation instantiations (cross-study and mixed) -----------
    P0 = p_model
    w21 = sqrt(A * cf)
    dw21 = w21 * ef / (2 * cf)
    P21 = y * w21
    dP21 = sqrt((dy * w21)**2 + (y * dw21)**2)
    # conservative variant: the shared sigma*t0 input enters both factors
    # (d*sqrt(s) ~ sqrt(sig t0); sqrt(sigma_phys) ~ sqrt(sig t0)) and is
    # propagated COHERENTLY; remaining parts (r*, w-stat, t0-phys) in quadrature
    rel_sigt0 = DSIGT0 / SIGT0                     # coherent, enters P linearly
    rel_rstar = DRSTAR / RSTAR_T0
    rel_wstat = ef_stat / cf / 2                   # from c-stat -> w-stat
    # sqrt(t0) = 0.1443(7)(13) fm (Bulava et al., their Ref. [22]);
    # enters sqrt(sigma_phys) as 1/sqrt(t0), independent of sigma*t0
    rel_t0phys = sqrt(mpf(7)**2 + mpf(13)**2) / mpf(1443)
    dP21c = P21 * sqrt(rel_sigt0**2 + rel_rstar**2 +
                       rel_wstat**2 + rel_t0phys**2)
    r['dP_qcd_cons'] = dP21c
    r['P_qcd_sigma_cons'] = None  # set after P0 known
    PB = y * yb
    dPB = sqrt((dy * yb)**2 + (y * dyb)**2)
    r['P_qcd_sigma_cons'] = abs(P0 - P21) / dP21c
    r.update(P0=P0, wsqs_qcd=w21, dwsqs_qcd=dw21,
             P_qcd=P21, dP_qcd=dP21, P_qcd_sigma=abs(P0 - P21) / dP21,
             P_baker=PB, dP_baker=dPB, P_baker_sigma=abs(P0 - PB) / dPB,
             P_mix_ratio=p_lat / P0, dP_mix_ratio=dp_lat / P0,
             xq_gE=p_lat / 2, dxq_gE=dp_lat / 2)

    # -- boundary-stiffness (Robin) reading of x_q -------------------------
    # x J1(x)/J0(x) = beta ; Dirichlet = beta -> inf ; massive exterior:
    # beta = y K1(y)/K0(y), y = R0/lambda (Clem penetration length)
    def beta_of(x):
        return x * besselj(1, x) / besselj(0, x)
    r['beta_central'] = beta_of(xq)
    r['beta_low1s'] = beta_of(xq - dxq)
    r['hardwall_sigma'] = r['xq_vs_dirichlet']        # same determination
    from mpmath import besselk
    def x_of_lambda(lam_fm):
        yv = R0_HEUR / mpf(lam_fm)
        rhs = yv * besselk(1, yv) / besselk(0, yv)
        return yv, findroot(lambda x: beta_of(x) - rhs, mpf('2.05'))
    r['robin_lam'] = {lam: x_of_lambda(lam) for lam in ('0.141', '0.174')}

    # -- dimensional outputs at the physical-point tension -----------------
    R0_out = sqrt(pi) * HB / SQS_PHYS
    r['R0_out'] = R0_out
    r['dR0_out'] = R0_out * DSQS_PHYS / SQS_PHYS
    db = 2 * x01 * HB / (sqrt(pi) * SQS_PHYS)
    r['d_break_out'] = db
    r['dd_break_out'] = db * DSQS_PHYS / SQS_PHYS
    return r


def report(r, out=print):
    """Human-readable dump of every manuscript number."""
    n = nstr
    out('=' * 74)
    out('Closed-form constants of the J0 Dirichlet mode (Table 1)')
    out('=' * 74)
    out('x_{0,1} = %s   x_{1,1} = %s   x_{2,1} = %s'
        % (n(r['x01'], 8), n(r['x11'], 8), n(r['x21'], 8)))
    out('A = 1-4/x01^2 = %s ; w/R0 = sqrt(A) = %s ; R0/w = %s'
        % (n(r['A'], 8), n(r['wR0'], 8), n(r['ratio'], 7)))
    out('R0 = %s fm (heuristic) -> m_perp = %s GeV ; sigma(c=pi) = %s GeV^2'
        % (n(r['R0'], 3), n(r['m_perp'], 6), n(r['sigma_pi'], 6)))
    out('sqrt(sigma) = %s GeV ; d_break(chiral) = %s fm ; r* = %s fm'
        % (n(r['sqrt_sigma_pi'], 6), n(r['d_break_pi'], 6), n(r['r_star'], 6)))
    out('d_break*sqrt(sigma)|_{c=pi} = %s ; w*sqrt(sigma)|_{c=pi} = %s'
        % (n(r['dsqs_model_pi'], 6), n(r['wsqs_model_pi'], 6)))
    out('w = %s fm ; Schwinger suppression exp(-%s) = %s'
        % (n(r['w_fm'], 5), n(r['schwinger_exp'], 4),
           n(exp(-r['schwinger_exp']), 3)))
    out('')
    out('=' * 74)
    out('Observable 1 -- string breaking calibrates c  (PLB 854, 138754)')
    out('=' * 74)
    out('d_break*sqrt(sigma) = %s +- %s'
        % (n(r['dsqs_latt'], 6), n(r['ddsqs_latt'], 4)))
    out('c_d = %s +- %s ; |c_d - pi| = %s sigma'
        % (n(r['cd'], 6), n(r['dcd'], 4), n(r['cd_vs_pi'], 3)))
    for rho, (dyx, sig) in sorted(r['rho_scan'].items()):
        out('  rho = %+d : dy = %s -> |c_d - pi| = %s sigma'
            % (rho, n(dyx, 4), n(sig, 3)))
    out('')
    out('=' * 74)
    out('Width table (PRD 95, 114511; per-ensemble scale setting)')
    out('=' * 74)
    out('%-4s %-6s %-5s %-12s %-14s %-12s'
        % ('th', 'beta', 'd', 'w [fm]', 'w*sqrt(s)', 'c implied'))
    for rw in r['width_rows']:
        out('%-4s %-6s %-5s %s(%s)   %s(%s)    %s(%s)'
            % (rw['th'], rw['beta'], rw['d'], n(rw['w'], 4), n(rw['dw'], 3),
               n(rw['y'], 4), n(rw['dy_tot'], 3),
               n(rw['c'], 4), n(rw['dc_tot'], 3)))
    out('range of w*sqrt(sigma) across all nine rows: %s -- %s'
        % (n(r['wsqs_range'][0], 4), n(r['wsqs_range'][1], 4)))
    out('quenched d=0.76: w*sqrt(sigma) = %s +- %s ; '
        'c_width^SU3 = %s +- %s (chi2/dof = %s/1)'
        % (n(r['wsqs_su3'], 4), n(r['dwsqs_su3'], 3),
           n(r['c_width_su3'], 4), n(r['dc_width_su3'], 3),
           n(r['chi2_su3'], 3)))
    out('full QCD d=0.76 (common 445(7) MeV scale error): '
        'c_width^{2+1} = %s +- %s  (stat %s, scale %s)'
        % (n(r['c_width_qcd'], 4), n(r['dc_width_qcd'], 3),
           n(r['dc_width_qcd_stat'], 3), n(r['dc_width_qcd_scale'], 3)))
    out('')
    out('=' * 74)
    out('The parameter-eliminated cross-observable test')
    out('=' * 74)
    out('PRODUCT RELATION (c-free): 2*sqrt(x01^2-4) = %s' % n(r['prod_model'], 8))
    out('  lattice product = %s +- %s  -> %s sigma'
        % (n(r['prod_lat'], 5), n(r['dprod_lat'], 3), n(r['prod_sigma'], 3)))
    out('PREDICT width from c_d : w*sqrt(sigma) = %s +- %s'
        % (n(r['wsqs_pred'], 4), n(r['dwsqs_pred'], 3)))
    out('MEASURED (quenched, d=0.76 fm)        = %s +- %s   -> %s sigma'
        % (n(r['wsqs_su3'], 4), n(r['dwsqs_su3'], 3), n(r['cross_sigma'], 3)))
    out('REVERSE: predict breaking from c_width = %s +- %s vs %s +- %s '
        '-> %s sigma'
        % (n(r['dsqs_pred'], 4), n(r['ddsqs_pred'], 3),
           n(r['dsqs_latt'], 4), n(r['ddsqs_latt'], 3),
           n(r['cross_sigma_rev'], 3)))
    out('c_d vs c_width^SU3 : %s sigma ; c_width^SU3 vs pi : %s sigma ; '
        'c_width^{2+1} vs pi : %s sigma ; QCD vs SU3 : %s sigma'
        % (n(r['cd_vs_cwidth'], 3), n(r['cwidth_su3_vs_pi'], 3),
           n(r['cwidth_qcd_vs_pi'], 3), n(r['cwidth_qcd_vs_su3'], 3)))
    out('')
    out('=' * 74)
    out('Full-QCD robustness check (Baker et al., EPJC 85, 29 (2025))')
    out('=' * 74)
    out('six widths at d = 0.723-0.760 fm, inverse-variance combined:')
    out('  w = %s +- %s fm   (chi2/dof = %s/5)'
        % (n(r['baker_w'], 4), n(r['baker_dw'], 3), n(r['baker_chi2'], 3)))
    out('  w*sqrt(sigma) [445(7) MeV, scale error correlated] = %s +- %s'
        % (n(r['baker_wsqs'], 4), n(r['baker_dwsqs'], 3)))
    out('  vs cross-prediction %s +- %s  ->  %s sigma'
        % (n(r['wsqs_pred'], 4), n(r['dwsqs_pred'], 3),
           n(r['baker_sigma'], 3)))
    out('  implied c = %s +- %s' % (n(r['baker_c'], 4), n(r['baker_dc'], 3)))
    out('')
    out('=' * 74)
    out('Boundary-condition sensitivity (Table 3) and x_q determination')
    out('=' * 74)
    for lab, x, mperp, dR0, dfm in r['bc_table']:
        out('%s  x = %-9s m_perp = %s GeV  d/R0 = %s  d_break = %s fm'
            % (lab, n(x, 6), n(mperp, 4), n(dR0, 4), n(dfm, 4)))
    out('data-side quark eigenvalue x_q = (1/2) d*sqrt(sigma)*sqrt(c_width) '
        '= %s +- %s' % (n(r['xq'], 4), n(r['dxq'], 3)))
    out('  vs scalar Dirichlet: %s sigma ; vs local MIT: %s sigma'
        % (n(r['xq_vs_dirichlet'], 3), n(r['xq_vs_mit'], 3)))
    out('product-relation instantiations (P0 = %s):' % n(r['P0'], 8))
    out('  cross-study Nf=2+1 (Cea rows): w*sqrt(s) = %s +- %s -> '
        'P = %s +- %s (uncorr) / +- %s (shared sig*t0 coherent) '
        '-> %s / %s sigma'
        % (n(r['wsqs_qcd'], 4), n(r['dwsqs_qcd'], 3),
           n(r['P_qcd'], 4), n(r['dP_qcd'], 3), n(r['dP_qcd_cons'], 3),
           n(r['P_qcd_sigma'], 3), n(r['P_qcd_sigma_cons'], 3)))
    out('  cross-study physical-mass (Baker, IVW; intra-study covariance '
        'unavailable): P = %s +- %s  (approx %s sigma under this treatment)'
        % (n(r['P_baker'], 4), n(r['dP_baker'], 3), n(r['P_baker_sigma'], 3)))
    out('  mixed (sharpest, B2-conditional): P/P0 = %s +- %s ; '
        'x_q*g_E = P/2 = %s +- %s'
        % (n(r['P_mix_ratio'], 4), n(r['dP_mix_ratio'], 3),
           n(r['xq_gE'], 4), n(r['dxq_gE'], 3)))
    out('boundary-stiffness (Robin) reading, x J1(x)/J0(x) = beta:')
    out('  beta(x_q = %s) = %s ; beta(x_q - 1sigma = %s) = %s ; '
        'hard-wall limit within %s sigma (same single determination)'
        % (n(r['xq'], 4), n(r['beta_central'], 3),
           n(r['xq'] - r['dxq'], 4), n(r['beta_low1s'], 3),
           n(r['hardwall_sigma'], 2)))
    for lam, (yv, xv) in sorted(r['robin_lam'].items()):
        out('  massive exterior, lambda = %s fm: y = R0/lambda = %s -> '
            'x = %s' % (lam, n(yv, 4), n(xv, 4)))
    out('')
    out('=' * 74)
    out('Dimensional outputs at the physical-point tension 0.445(7) GeV')
    out('=' * 74)
    out('R0 (c=pi)      = %s +- %s fm'
        % (n(r['R0_out'], 5), n(r['dR0_out'], 3)))
    out('d_break (c=pi) = %s +- %s fm  vs reported 1.211 +- 0.013 fm'
        % (n(r['d_break_out'], 5), n(r['dd_break_out'], 3)))


if __name__ == '__main__':
    report(compute())

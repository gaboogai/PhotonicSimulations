# Fit analytic functions r_magnitude(kx, ky) and r_phase(kx, ky) to the field
# reflectivity of a photonic-crystal slab, from the long-format k-space CSV
# written by kSpaceConvert.py (columns:
#   kx,ky,k_parallel,polar_deg,azi_deg,magnitude,phase_deg).
#
# WHAT FUNCTION IS EXPECTED
# The magnitude has an interior peak pinned at unity (|r| -> 1) on a high
# background, i.e. a guided-mode / Fano resonance of the slab; the phase rolls
# smoothly through it. The physically-expected model is therefore a single
# coupled-mode-theory (Fano) resonance whose detuning is quadratic in the
# in-plane wavevector, kappa^2 = kx^2 + ky^2 (this ONE complex function gives
# both magnitude and phase). The hexagonal lattice adds a weak 6-fold term, the
# lowest C6v invariant Re[(kx + i ky)^6]. Because a clean single answer is only
# approximate, two "reasonable estimate" families are ALSO fit and compared:
#   A) coupled-mode / Fano resonance in kappa^2  (physically expected, headline)
#      - complex CMT:  r = r_bg + c / (1 + i (kappa^2 - kr^2)/gamma)
#      - lossless transmission-Fano magnitude:  |t|^2 Fano, |r| = sqrt(1-|t|^2)
#   B) symmetry-adapted polynomial surrogates in {1, rho^2, rho^4, rho^6, h6}
#      fit separately to magnitude and phase (fast, robust, purely descriptive).
#
# Fits are done in a normalized wavevector k~ = k / k_scale (k_scale = max
# k_parallel in the data) for conditioning; the exported callables accept the
# SAME kx, ky units as the input CSV and normalize internally.
#
# Usage:
#   python fitPhotonicTest.py  path/to/...KSpace.csv
#
# Reusable:
#   from fitPhotonicTest import load_callables
#   fns = load_callables('..._fitParams.json')
#   fns['r_magnitude'](kx, ky), fns['r_phase'](kx, ky), fns['r_complex'](kx, ky)

import json
import os
import sys

import numpy as np
from scipy.optimize import curve_fit, least_squares

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')


def unique_path(path):
    """Return path, or path with an incrementing _N suffix if it already exists,
    so existing output files are never overwritten."""
    base, ext = os.path.splitext(path)
    candidate = path
    n = 1
    while os.path.exists(candidate):
        candidate = f'{base}_{n}{ext}'
        n += 1
    return candidate


# ============================================================
#   MODEL EVALUATION (all in NORMALIZED wavevector k~ = k/k_scale)
# ============================================================

def hex6(kxn, kyn):
    """Lowest-order C6v-invariant anisotropy term Re[(kx + i ky)^6]."""
    return np.real((kxn + 1j * kyn) ** 6)


def poly_design(kxn, kyn, n_rho, include_h6):
    """Design matrix over the symmetry-adapted basis {1, rho^2, rho^4, ...,
    rho^(2 n_rho)} (+ the 6-fold term). rho^2 = kx^2 + ky^2 (normalized)."""
    rho2 = kxn ** 2 + kyn ** 2
    cols = [np.ones_like(rho2)]
    for p in range(1, n_rho + 1):
        cols.append(rho2 ** p)
    if include_h6:
        cols.append(hex6(kxn, kyn))
    return np.column_stack(cols)


def poly_eval(coeffs, kxn, kyn, n_rho, include_h6):
    return poly_design(kxn, kyn, n_rho, include_h6) @ coeffs


def fano_t2(rho2, T0, kz2, kr2, gamma):
    """Transmission Fano lineshape (single guided resonance, transmission zero
    at kappa^2 = kz2, resonance at kappa^2 = kr2)."""
    return T0 * (rho2 - kz2) ** 2 / ((rho2 - kr2) ** 2 + gamma ** 2)


def fano_mag(rho2, T0, kz2, kr2, gamma):
    """Lossless-slab magnitude |r| = sqrt(1 - |t|^2), pinned to <= 1."""
    return np.sqrt(np.clip(1.0 - fano_t2(rho2, T0, kz2, kr2, gamma), 0.0, 1.0))


def cmt_r(kxn, kyn, rbg, c, kr2, gamma, a6=0.0):
    """Complex coupled-mode-theory reflection: a single Lorentzian resonance in
    kappa^2 on a complex background, with an optional 6-fold shift of the
    resonance position (a6)."""
    rho2 = kxn ** 2 + kyn ** 2
    kr2_eff = kr2 - a6 * hex6(kxn, kyn)
    return rbg + c / (1.0 + 1j * (rho2 - kr2_eff) / gamma)


# ============================================================
#   GOODNESS OF FIT
# ============================================================

def gof(data, model):
    """Return (R^2, RMS residual, max abs residual)."""
    resid = np.asarray(model) - np.asarray(data)
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((data - np.mean(data)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return r2, float(np.sqrt(np.mean(resid ** 2))), float(np.max(np.abs(resid)))


# ============================================================
#   READ THE k-SPACE CSV
# ============================================================

def read_kspace(path):
    """Read the long-format k-space CSV (named header). Returns a dict of
    columns as float arrays."""
    with open(path) as f:
        header = f.readline().strip().split(',')
    cols = {name: [] for name in header}
    with open(path) as f:
        next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            vals = line.split(',')
            for name, v in zip(header, vals):
                cols[name].append(float(v))
    return {name: np.array(v) for name, v in cols.items()}


# ============================================================
#   THE FIT
# ============================================================

def fit_all(csv_path):
    d = read_kspace(csv_path)
    kx, ky = d['kx'], d['ky']
    kpar = d['k_parallel']
    mag, phase = d['magnitude'], d['phase_deg']
    polar = d['polar_deg']

    k_scale = float(np.max(kpar))
    kxn, kyn = kx / k_scale, ky / k_scale
    rho2 = kxn ** 2 + kyn ** 2
    rho2_max = float(np.max(rho2))
    z = mag * np.exp(1j * np.deg2rad(phase))          # complex data
    polar_max = float(np.max(polar))

    results = {'k_scale': k_scale, 'polar_max_deg': polar_max,
               'n_points': int(mag.size)}

    # ----- Model B: polynomial surrogates (linear least-squares) -----
    A_mag = poly_design(kxn, kyn, n_rho=3, include_h6=True)   # 1,r2,r4,r6,h6
    mag_coeffs = np.linalg.lstsq(A_mag, mag, rcond=None)[0]
    mag_poly = A_mag @ mag_coeffs

    A_ph = poly_design(kxn, kyn, n_rho=2, include_h6=True)    # 1,r2,r4,h6
    ph_coeffs = np.linalg.lstsq(A_ph, phase, rcond=None)[0]
    ph_poly = A_ph @ ph_coeffs

    results['poly'] = {
        'mag_coeffs': mag_coeffs.tolist(), 'mag_basis': '1,rho2,rho4,rho6,h6',
        'phase_coeffs': ph_coeffs.tolist(), 'phase_basis': '1,rho2,rho4,h6',
        'mag_gof': gof(mag, mag_poly), 'phase_gof': gof(phase, ph_poly)}

    # ----- Model A (magnitude): lossless transmission-Fano -----
    kr2_0 = float(rho2[np.argmax(mag)])               # peak location
    t2_far = float(np.max(1.0 - mag ** 2))            # largest |t|^2 in data
    fano = None
    try:
        p0 = [t2_far, kr2_0, kr2_0, 0.3]
        bounds = ([0.0, 0.0, 0.0, 1e-8], [100.0, 4.0, 4.0, 100.0])
        popt, _ = curve_fit(fano_mag, rho2, mag, p0=p0, bounds=bounds,
                            maxfev=200000)
        fano = dict(T0=popt[0], kz2=popt[1], kr2=popt[2], gamma=popt[3])
        # |r| peaks (-> 1) at the transmission zero kappa^2 = kz2, so that is the
        # reflection-resonance angle (not kr2, the |t|^2 denominator centre).
        res_polar = float(np.degrees(np.arcsin(
            np.sqrt(max(popt[1], 0.0)) * np.sin(np.radians(polar_max)))))
        identifiable = (popt[3] <= 1.5 * rho2_max
                        and 0.0 <= popt[1] <= 1.2 * rho2_max)
        results['fano_mag'] = {
            'params': fano, 'gof': gof(mag, fano_mag(rho2, *popt)),
            'resonance_polar_deg': res_polar,
            'identifiable': bool(identifiable),
            'caveat': ('' if identifiable else
                       f"gamma={popt[3]:.2f} >> data span {rho2_max:.2f} in "
                       f"kappa^2; resonance width not resolved by this window")}
    except Exception as e:                            # noqa: BLE001
        results['fano_mag'] = {'error': str(e)}

    # ----- Model A (complex): coupled-mode resonance, joint mag+phase -----
    # multi-start over the resonance width for robustness
    i_peak = int(np.argmax(mag))
    i_far = int(np.argmax(rho2))
    rbg0 = z[i_far]
    best = None
    for g0 in (0.1, 0.3, 0.6, 1.0):
        c0 = z[i_peak] - rbg0
        p0 = [rbg0.real, rbg0.imag, c0.real, c0.imag, kr2_0, g0]
        lb = [-2, -2, -2, -2, 0.0, 1e-6]
        ub = [2, 2, 2, 2, 4.0, 100.0]

        def resid(p):
            m = cmt_r(kxn, kyn, p[0] + 1j * p[1], p[2] + 1j * p[3], p[4], p[5])
            return np.concatenate([(m - z).real, (m - z).imag])

        try:
            sol = least_squares(resid, p0, bounds=(lb, ub), max_nfev=200000)
            if best is None or sol.cost < best.cost:
                best = sol
        except Exception:                             # noqa: BLE001
            continue

    # ----- Model A (complex) + 6-fold anisotropy -----
    if best is not None:
        p = best.x
        cmt = dict(rbg=complex(p[0], p[1]), c=complex(p[2], p[3]),
                   kr2=p[4], gamma=p[5], a6=0.0)
        m = cmt_r(kxn, kyn, cmt['rbg'], cmt['c'], cmt['kr2'], cmt['gamma'])
        results['cmt'] = _cmt_report(cmt, mag, phase, m, polar_max, rho2_max)

        def resid6(q):
            mm = cmt_r(kxn, kyn, q[0] + 1j * q[1], q[2] + 1j * q[3],
                       q[4], q[5], q[6])
            return np.concatenate([(mm - z).real, (mm - z).imag])

        try:
            q0 = list(p) + [0.0]
            lb6, ub6 = [-2, -2, -2, -2, 0.0, 1e-6, -4], [2, 2, 2, 2, 4.0, 100.0, 4]
            sol6 = least_squares(resid6, q0, bounds=(lb6, ub6), max_nfev=200000)
            q = sol6.x
            cmt6 = dict(rbg=complex(q[0], q[1]), c=complex(q[2], q[3]),
                        kr2=q[4], gamma=q[5], a6=q[6])
            m6 = cmt_r(kxn, kyn, cmt6['rbg'], cmt6['c'], cmt6['kr2'],
                       cmt6['gamma'], cmt6['a6'])
            results['cmt_hex'] = _cmt_report(cmt6, mag, phase, m6, polar_max,
                                             rho2_max)
        except Exception as e:                        # noqa: BLE001
            results['cmt_hex'] = {'error': str(e)}

    results['_data'] = dict(kxn=kxn, kyn=kyn, rho2=rho2, mag=mag, phase=phase,
                            polar=polar, azi=d['azi_deg'], z=z)
    return results


def _cmt_report(cmt, mag, phase, model, polar_max, rho2_max):
    model_mag = np.abs(model)
    model_phase = np.degrees(np.angle(model))
    res_polar = float(np.degrees(np.arcsin(
        np.sqrt(max(cmt['kr2'], 0.0)) * np.sin(np.radians(polar_max)))))
    bg = float(np.abs(cmt['rbg']))
    # The resonance parameters are only identifiable if the resonance is at
    # least partly resolved by the swept window: its width gamma must not vastly
    # exceed the data span in kappa^2, its centre kr2 must sit in/near the data,
    # and the background must be physical (|r_bg| ~ 1 for a high-reflect slab).
    identifiable = (cmt['gamma'] <= 1.5 * rho2_max
                    and 0.0 <= cmt['kr2'] <= 1.2 * rho2_max
                    and 0.5 <= bg <= 1.05)
    reason = ''
    if not identifiable:
        if cmt['gamma'] > 1.5 * rho2_max:
            reason = (f"resonance width gamma={cmt['gamma']:.2f} >> data span "
                      f"{rho2_max:.2f} in kappa^2 (window sees a near-linear "
                      f"slice; kr2/gamma/Q not recoverable)")
        elif not 0.5 <= bg <= 1.05:
            reason = f"non-physical background |r_bg|={bg:.3f}"
        else:
            reason = f"resonance centre kr2={cmt['kr2']:.3f} outside data"
    return {
        'params': {'rbg_re': cmt['rbg'].real, 'rbg_im': cmt['rbg'].imag,
                   'c_re': cmt['c'].real, 'c_im': cmt['c'].imag,
                   'kr2': cmt['kr2'], 'gamma': cmt['gamma'], 'a6': cmt['a6']},
        'mag_gof': gof(mag, model_mag), 'phase_gof': gof(phase, model_phase),
        'background_mag': bg,
        'resonance_polar_deg': res_polar,
        'angular_sharpness': float(cmt['kr2'] / (2.0 * cmt['gamma'])),
        'identifiable': bool(identifiable), 'caveat': reason}


# ============================================================
#   CALLABLE BUILDERS (accept kx, ky in the INPUT CSV's units)
# ============================================================

def _build_callables(params):
    """Build kx,ky -> value callables. The default r_magnitude / r_phase use
    whichever model was selected (params['mag_model'], params['phase_model']);
    named per-model variants are always available. kx, ky are in the input
    CSV's units (rad/um for K_UNITS='physical')."""
    k_scale = params['k_scale']
    poly = params['poly']
    cmt = params.get('cmt_hex_params') or params.get('cmt_params')
    fano = params.get('fano_params')

    def _n(kx, ky):
        return np.asarray(kx, float) / k_scale, np.asarray(ky, float) / k_scale

    # --- per-model magnitude evaluators ---
    def r_magnitude_poly(kx, ky):
        kxn, kyn = _n(kx, ky)
        return poly_eval(np.array(poly['mag_coeffs']), kxn, kyn, 3, True)

    def r_magnitude_cmt(kx, ky):
        kxn, kyn = _n(kx, ky)
        return np.abs(cmt_r(kxn, kyn, complex(cmt['rbg_re'], cmt['rbg_im']),
                            complex(cmt['c_re'], cmt['c_im']),
                            cmt['kr2'], cmt['gamma'], cmt.get('a6', 0.0)))

    def r_magnitude_fano(kx, ky):
        kxn, kyn = _n(kx, ky)
        return fano_mag(kxn ** 2 + kyn ** 2, fano['T0'], fano['kz2'],
                        fano['kr2'], fano['gamma'])

    # --- per-model phase evaluators ---
    def r_phase_poly(kx, ky):
        kxn, kyn = _n(kx, ky)
        return poly_eval(np.array(poly['phase_coeffs']), kxn, kyn, 2, True)

    def r_phase_cmt(kx, ky):
        kxn, kyn = _n(kx, ky)
        return np.degrees(np.angle(cmt_r(
            kxn, kyn, complex(cmt['rbg_re'], cmt['rbg_im']),
            complex(cmt['c_re'], cmt['c_im']),
            cmt['kr2'], cmt['gamma'], cmt.get('a6', 0.0))))

    mag_fns = {'poly': r_magnitude_poly, 'cmt': r_magnitude_cmt,
               'cmt_hex': r_magnitude_cmt, 'fano': r_magnitude_fano}
    phase_fns = {'poly': r_phase_poly, 'cmt': r_phase_cmt,
                 'cmt_hex': r_phase_cmt}

    r_magnitude = mag_fns[params.get('mag_model', 'poly')]
    r_phase = phase_fns[params.get('phase_model', 'poly')]

    def r_complex(kx, ky):
        return r_magnitude(kx, ky) * np.exp(1j * np.deg2rad(r_phase(kx, ky)))

    out = {'r_complex': r_complex, 'r_magnitude': r_magnitude,
           'r_phase': r_phase, 'r_magnitude_poly': r_magnitude_poly,
           'r_phase_poly': r_phase_poly}
    if cmt is not None:
        out['r_magnitude_cmt'] = r_magnitude_cmt
        out['r_phase_cmt'] = r_phase_cmt
    if fano is not None:
        out['r_magnitude_fano'] = r_magnitude_fano
    return out


def load_callables(json_path):
    """Rebuild the fitted callables from a saved ..._fitParams.json."""
    with open(json_path) as f:
        params = json.load(f)
    return _build_callables(params)


# ============================================================
#   PLOTTING
# ============================================================

def _to_grid(polar, azi, q):
    pu, au = np.unique(polar), np.unique(azi)
    g = np.full((len(pu), len(au)), np.nan)
    ip = np.searchsorted(pu, polar)
    ja = np.searchsorted(au, azi)
    g[ip, ja] = q
    return pu, au, g


def make_plots(results, outfile, fns):
    D = results['_data']
    polar, azi = D['polar'], D['azi']
    k_scale = results['k_scale']
    polar_max = results['polar_max_deg']
    au = np.unique(azi)
    azi_lo = au[0]
    azi_mid = au[np.argmin(np.abs(au - 30.0))]

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # fine polar line, evaluate models (convert polar,azi -> physical k like the
    # CSV did: k_parallel = k_scale * sin(polar)/sin(polar_max))
    pf = np.linspace(0.0, polar_max, 300)
    s_ratio = np.sin(np.deg2rad(pf)) / np.sin(np.deg2rad(polar_max))

    for col, azi_val in enumerate((azi_lo, azi_mid)):
        m = np.isclose(azi, azi_val)
        order = np.argsort(polar[m])
        pol_s = polar[m][order]

        kpar_f = k_scale * s_ratio
        kx_f = kpar_f * np.cos(np.deg2rad(azi_val))
        ky_f = kpar_f * np.sin(np.deg2rad(azi_val))

        # magnitude
        ax[col, 0].scatter(pol_s, D['mag'][m][order], s=14, color='k',
                           zorder=3, label='data')
        ax[col, 0].plot(pf, fns['r_magnitude'](kx_f, ky_f), color='C0',
                       label='CMT/Fano')
        ax[col, 0].plot(pf, fns['r_magnitude_poly'](kx_f, ky_f), '--',
                       color='C1', label='polynomial')
        ax[col, 0].set_xlabel('Polar angle (deg)')
        ax[col, 0].set_ylabel('|r|')
        ax[col, 0].set_title(f'|r| vs polar at azi = {azi_val:g} deg')
        ax[col, 0].legend(fontsize=8)

        # phase
        ax[col, 1].scatter(pol_s, D['phase'][m][order], s=14, color='k',
                           zorder=3, label='data')
        ax[col, 1].plot(pf, fns['r_phase'](kx_f, ky_f), color='C0',
                       label='CMT')
        ax[col, 1].plot(pf, fns['r_phase_poly'](kx_f, ky_f), '--', color='C1',
                       label='polynomial')
        ax[col, 1].set_xlabel('Polar angle (deg)')
        ax[col, 1].set_ylabel('arg(r) (deg)')
        ax[col, 1].set_title(f'arg(r) vs polar at azi = {azi_val:g} deg')
        ax[col, 1].legend(fontsize=8)

    # residual heatmaps (headline model) over the (azi, polar) grid
    mag_res = fns['r_magnitude'](D['kxn'] * k_scale, D['kyn'] * k_scale) - D['mag']
    ph_res = fns['r_phase'](D['kxn'] * k_scale, D['kyn'] * k_scale) - D['phase']
    for row, (res, label) in enumerate(((mag_res, '|r| residual'),
                                        (ph_res, 'phase residual (deg)'))):
        pu, auu, g = _to_grid(polar, azi, res)
        vmax = np.nanmax(np.abs(g)) or 1.0
        mesh = ax[row, 2].pcolormesh(auu, pu, g, shading='auto', cmap='RdBu_r',
                                     vmin=-vmax, vmax=vmax)
        fig.colorbar(mesh, ax=ax[row, 2])
        ax[row, 2].set_xlabel('Azimuthal angle (deg)')
        ax[row, 2].set_ylabel('Polar angle (deg)')
        ax[row, 2].set_title(f'selected-model {label}')

    fig.suptitle(os.path.basename(outfile))
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)


# ============================================================
#   MAIN
# ============================================================

def _fmt_gof(g):
    return f'R2={g[0]:.6f}  RMS={g[1]:.3e}  max={g[2]:.3e}'


if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.exit('usage: python fitPhotonicTest.py path/to/...KSpace.csv')

    csv_path = sys.argv[1]
    results = fit_all(csv_path)

    # ----- collect per-quantity candidates: (model_name, R^2, gof) -----
    mag_cands = [('poly', results['poly']['mag_gof'])]
    phase_cands = [('poly', results['poly']['phase_gof'])]
    if 'cmt' in results:
        mag_cands.append(('cmt', results['cmt']['mag_gof']))
        phase_cands.append(('cmt', results['cmt']['phase_gof']))
    if results.get('cmt_hex', {}).get('params'):
        mag_cands.append(('cmt_hex', results['cmt_hex']['mag_gof']))
        phase_cands.append(('cmt_hex', results['cmt_hex']['phase_gof']))
    if 'params' in results.get('fano_mag', {}):
        mag_cands.append(('fano', results['fano_mag']['gof']))

    # pick the best-fitting model per quantity (highest R^2)
    mag_model = max(mag_cands, key=lambda t: t[1][0])[0]
    phase_model = max(phase_cands, key=lambda t: t[1][0])[0]

    # assemble the JSON-serializable parameter record + selected callables
    params = {'k_scale': results['k_scale'],
              'polar_max_deg': results['polar_max_deg'],
              'source_csv': os.path.basename(csv_path),
              'mag_model': mag_model, 'phase_model': phase_model,
              'poly': results['poly']}
    if 'cmt' in results:
        params['cmt_params'] = results['cmt']['params']
    if results.get('cmt_hex', {}).get('params'):
        params['cmt_hex_params'] = results['cmt_hex']['params']
    if 'params' in results.get('fano_mag', {}):
        params['fano_params'] = results['fano_mag']['params']
    fns = _build_callables(params)

    # ----- report -----
    print(f"\nsource            : {os.path.basename(csv_path)}")
    print(f"points            : {results['n_points']}  "
          f"(polar <= {results['polar_max_deg']:g} deg)")
    print(f"k_scale (max kpar): {results['k_scale']:.6f}")

    print("\n--- Model A: coupled-mode / Fano resonance (physically expected) ---")
    for key, label in (('cmt', 'complex CMT (isotropic)'),
                       ('cmt_hex', 'complex CMT + 6-fold Re[(kx+i ky)^6]')):
        c = results.get(key, {})
        if 'params' not in c:
            continue
        print(f"  {label}:")
        print(f"    magnitude   {_fmt_gof(c['mag_gof'])}")
        print(f"    phase       {_fmt_gof(c['phase_gof'])}")
        if c['identifiable']:
            print(f"    background |r_bg| = {c['background_mag']:.5f}, "
                  f"resonance polar = {c['resonance_polar_deg']:.4f} deg, "
                  f"angular sharpness = {c['angular_sharpness']:.3f}")
        else:
            print(f"    [resonance params UNDER-CONSTRAINED: {c['caveat']}]")
    fm = results.get('fano_mag', {})
    if 'params' in fm:
        print(f"  lossless transmission-Fano magnitude:")
        tag = (f"resonance polar = {fm['resonance_polar_deg']:.4f} deg"
               if fm['identifiable'] else f"UNDER-CONSTRAINED: {fm['caveat']}")
        print(f"    magnitude   {_fmt_gof(fm['gof'])}  ({tag})")

    print("\n--- Model B: symmetry-adapted polynomial surrogates ---")
    print(f"  magnitude  (1,rho2,rho4,rho6,h6)  {_fmt_gof(results['poly']['mag_gof'])}")
    print(f"  phase      (1,rho2,rho4,h6)       {_fmt_gof(results['poly']['phase_gof'])}")

    # ----- save params -----
    stem = csv_path.rsplit('.', 1)[0]
    if stem.endswith('KSpace'):
        stem = stem[:-len('KSpace')]
    json_path = unique_path(stem + '_fitParams.json')
    with open(json_path, 'w') as f:
        json.dump(params, f, indent=2)
    print(f"\nsaved parameters  : {json_path}")

    # ----- plots -----
    png_path = unique_path(stem + '_fit.png')
    make_plots(results, png_path, fns)
    print(f"saved diagnostics : {png_path}")

    # ----- what-function-is-expected + data-driven recommendation -----
    cmt_ok = results.get('cmt_hex', {}).get('identifiable') or \
        results.get('cmt', {}).get('identifiable')
    print("\nExpected function : guided-mode / Fano resonance in kappa^2=kx^2+ky^2")
    print("  r(kx,ky) = r_bg + c / (1 + i(kappa^2 - kr^2)/gamma), plus a weak "
          "6-fold\n  Re[(kx+i ky)^6] hexagonal term (jointly gives |r| and arg(r)).")
    if not cmt_ok:
        print("  -> This sweep only spans a near-linear SLICE of that resonance "
              "(gamma\n     >> data span), so kr2/gamma/Q are not recoverable "
              "here. Widen the\n     photonicTest.py polar range to pin them.")
    print(f"\nSelected for r_magnitude(kx,ky): '{mag_model}'   "
          f"r_phase(kx,ky): '{phase_model}'   (best R^2 on this data)")
    print(f"NOTE: fit domain is polar <= {results['polar_max_deg']:g} deg, "
          "azimuthal within the swept range only.")

# Extend a ...KSpace.csv (from kSpaceConvert.py) around a full 360 degrees of
# azimuth, by tiling the swept wedge with the lattice's rotational symmetry.
#
# photonicTest.py only sweeps one symmetry period in azimuth - 60 deg for a
# hexagonal lattice (C6v), 90 deg for a square one (C4v) - because the
# reflectivity repeats beyond that. This script fills in the rest of the circle
# so plotKSpace.py shows a full disc instead of a pie slice.
#
# Why plain rotation is enough: S4 rotates the E,H,k frame with the azimuth, so
# the structure AND the incident polarisation rotate together. The co-polarised
# r is therefore invariant under a rotation by the lattice period - magnitude
# and phase both carry over unchanged, with no mirror operations needed.
#
# The period is detected from the azimuthal span of the file (60 -> hex,
# 90 -> square) and then checked: the slice at the first azimuth must match the
# slice at the last one, since those are the same physical direction. That check
# is what catches a mis-typed file (e.g. a square sweep run over 60 deg, which
# would otherwise look hexagonal), because for a square lattice r(0) != r(60).
#
# Input/output CSV layout (identical, so the result is a drop-in for
# plotKSpace.py and fitPhotonicTest.py):
#   kx,ky,k_parallel,polar_deg,azi_deg,magnitude,phase_deg
#
# Usage:
#   python extendKspace.py path/to/...KSpace.csv [output.csv] [plot=true] [period=60]

import os
import subprocess
import sys

import numpy as np


# --------------------------------------------------------
#   PARAMETERS
# --------------------------------------------------------

# azimuthal span (deg) -> lattice whose rotational symmetry it corresponds to
SUPPORTED_PERIODS = {60.0: 'hex', 90.0: 'square'}

SPAN_TOL = 1e-6      # deg, tolerance when matching the span to a known period

# the seam (first vs last azimuth) must agree to within this fraction of the
# data's own range before the span is accepted as a true symmetry period
SEAM_RTOL = 1e-3

# drop rows that land on an identical (kx, ky); mainly the origin, where every
# azimuth at polar = 0 collapses to the same point and would otherwise be
# stacked once per tile (and over-weighted in any least-squares fit)
DEDUPLICATE = True


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


def read_kspace(path):
    """Read the long-format k-space CSV (named header) into a dict of arrays."""
    with open(path) as f:
        header = f.readline().strip().split(',')
    cols = {name: [] for name in header}
    with open(path) as f:
        next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            for name, v in zip(header, line.split(',')):
                cols[name].append(float(v))
    return {name: np.array(v) for name, v in cols.items()}


def seam_error(polar, values, m_lo, m_hi):
    """Max |difference| between the first- and last-azimuth slices, paired by
    polar angle. Returns None if the two slices do not share a polar grid."""
    o_lo = np.argsort(polar[m_lo])
    o_hi = np.argsort(polar[m_hi])
    p_lo, p_hi = polar[m_lo][o_lo], polar[m_hi][o_hi]
    if len(p_lo) != len(p_hi) or not np.allclose(p_lo, p_hi):
        return None
    return float(np.max(np.abs(values[m_lo][o_lo] - values[m_hi][o_hi])))


# --------------------------------------------------------
#   1/ ARGUMENTS
# --------------------------------------------------------

# positional: input.csv [output.csv]; optional plot=true / period=<deg> flags
positional = []
do_plot = False
period_override = None
for a in sys.argv[1:]:
    low = a.lower()
    if low.startswith('plot='):
        do_plot = low.split('=', 1)[1].strip() in ('true', '1', 'yes', 'on')
    elif low.startswith('period='):
        try:
            period_override = float(low.split('=', 1)[1])
        except ValueError:
            sys.exit(f'period= must be a number in degrees, got {a!r}')
    else:
        positional.append(a)

if len(positional) < 1:
    sys.exit('usage: python extendKspace.py path/to/...KSpace.csv '
             '[output.csv] [plot=true] [period=60]')

infile = positional[0]

# --------------------------------------------------------
#   2/ READ AND DETECT THE SYMMETRY PERIOD
# --------------------------------------------------------

d = read_kspace(infile)
for col in ('kx', 'ky', 'k_parallel', 'polar_deg', 'azi_deg', 'magnitude',
            'phase_deg'):
    if col not in d:
        sys.exit(f'{infile} is missing the {col!r} column - is it a '
                 'kSpaceConvert.py ...KSpace.csv file?')

azi, polar, kpar = d['azi_deg'], d['polar_deg'], d['k_parallel']
mag, phase = d['magnitude'], d['phase_deg']

azi_u = np.unique(azi)
azi_min, azi_max = float(azi_u[0]), float(azi_u[-1])
span = azi_max - azi_min

if period_override is not None:
    period = period_override
    lattice = SUPPORTED_PERIODS.get(period, 'custom')
    source = 'from period= override'
else:
    match = [p for p in SUPPORTED_PERIODS if abs(span - p) < SPAN_TOL]
    if not match:
        sys.exit(
            f'unsupported azimuthal span: {span:g} deg '
            f'(azi {azi_min:g} to {azi_max:g}).\n'
            f'  extendKspace.py expects one full symmetry period: '
            f'{", ".join(f"{p:g} deg ({n})" for p, n in SUPPORTED_PERIODS.items())}.\n'
            f'  Re-run photonicTest.py over a full period, or force one with '
            f'period=<deg>.')
    period = match[0]
    lattice = SUPPORTED_PERIODS[period]
    source = 'auto-detected from the azimuthal span'

n_tiles = 360.0 / period
if abs(n_tiles - round(n_tiles)) > 1e-9:
    sys.exit(f'period {period:g} deg does not divide 360 deg evenly')
n_tiles = int(round(n_tiles))

# --------------------------------------------------------
#   3/ VALIDATE THE SEAM (first azimuth must equal last azimuth)
# --------------------------------------------------------

m_lo = np.isclose(azi, azi_min)
m_hi = np.isclose(azi, azi_max)

mag_seam = seam_error(polar, mag, m_lo, m_hi)
phase_seam = seam_error(polar, phase, m_lo, m_hi)
if mag_seam is None or phase_seam is None:
    sys.exit('the first and last azimuth slices do not share the same polar '
             'grid - the file does not look like a regular photonicTest sweep')

mag_range = float(np.max(mag) - np.min(mag))
phase_range = float(np.max(phase) - np.min(phase))
mag_tol = max(1e-9, SEAM_RTOL * mag_range)
phase_tol = max(1e-9, SEAM_RTOL * phase_range)

if mag_seam > mag_tol or phase_seam > phase_tol:
    sys.exit(
        f'seam check FAILED: azi {azi_min:g} and azi {azi_max:g} deg should be '
        f'the same physical direction\n'
        f'  but differ by {mag_seam:.3e} in |r| (tolerance {mag_tol:.3e}) and '
        f'{phase_seam:.3e} deg in phase (tolerance {phase_tol:.3e}).\n'
        f'  A span of {span:g} deg is therefore NOT a symmetry period of this '
        f'structure - it may be the wrong lattice type\n'
        f'  (a square lattice needs 90 deg, a hexagonal one 60 deg). Override '
        f'with period=<deg> if you are sure.')

# --------------------------------------------------------
#   4/ TILE AROUND THE FULL 360 DEGREES
# --------------------------------------------------------

# keep [azi_min, azi_max): the last slice repeats the first, so including it
# would duplicate a seam once per tile
inner = azi < (azi_max - SPAN_TOL)

azi_out, kpar_out, polar_out, mag_out, phase_out = [], [], [], [], []
for k in range(n_tiles):
    azi_out.append((azi[inner] + k * period) % 360.0)
    kpar_out.append(kpar[inner])
    polar_out.append(polar[inner])
    mag_out.append(mag[inner])
    phase_out.append(phase[inner])

azi_out = np.concatenate(azi_out)
kpar_out = np.concatenate(kpar_out)
polar_out = np.concatenate(polar_out)
mag_out = np.concatenate(mag_out)
phase_out = np.concatenate(phase_out)

# rebuild kx, ky from k_parallel and the new azimuth (rather than rotating the
# stored kx, ky) so precision and units match kSpaceConvert.py exactly
kx_out = kpar_out * np.cos(np.deg2rad(azi_out))
ky_out = kpar_out * np.sin(np.deg2rad(azi_out))

n_tiled = len(kx_out)

if DEDUPLICATE:
    _, keep = np.unique(np.round(np.column_stack([kx_out, ky_out]), 12),
                        axis=0, return_index=True)
    keep = np.sort(keep)
    kx_out, ky_out, kpar_out = kx_out[keep], ky_out[keep], kpar_out[keep]
    polar_out, azi_out = polar_out[keep], azi_out[keep]
    mag_out, phase_out = mag_out[keep], phase_out[keep]

# --------------------------------------------------------
#   5/ WRITE THE EXTENDED CSV
# --------------------------------------------------------

if len(positional) > 1:
    outfile = unique_path(positional[1])
else:
    stem = infile.rsplit('.', 1)[0]
    outfile = unique_path(stem + ('360.csv' if stem.endswith('KSpace')
                                  else 'KSpace360.csv'))

os.makedirs(os.path.dirname(outfile) or '.', exist_ok=True)

with open(outfile, 'w') as f:
    f.write('kx,ky,k_parallel,polar_deg,azi_deg,magnitude,phase_deg\n')
    for i in range(len(kx_out)):
        f.write(f'{kx_out[i]},{ky_out[i]},{kpar_out[i]},{polar_out[i]},'
                f'{azi_out[i]},{mag_out[i]},{phase_out[i]}\n')

print(f'read {len(azi)} rows from {infile}')
print(f'azimuthal span   : {azi_min:g} to {azi_max:g} deg ({span:g} deg)')
print(f'symmetry period  : {period:g} deg -> {lattice} lattice, '
      f'{n_tiles} tiles ({source})')
print(f'seam check       : |r| {mag_seam:.3e} (tol {mag_tol:.3e}), '
      f'phase {phase_seam:.3e} deg (tol {phase_tol:.3e})')
print(f'azimuthal variation in the data: |r| {mag_range:.3e}, '
      f'phase {phase_range:.3e} deg')
print(f'tiled            : {n_tiled} rows'
      + (f' -> {len(kx_out)} after removing duplicate (kx, ky)'
         if DEDUPLICATE else ''))
print(f'azimuth coverage : {np.min(azi_out):g} to {np.max(azi_out):g} deg, '
      f'{len(np.unique(np.round(azi_out, 9)))} unique')
print(f'kx range         : {np.min(kx_out):.5f} - {np.max(kx_out):.5f}')
print(f'ky range         : {np.min(ky_out):.5f} - {np.max(ky_out):.5f}')
print(f'wrote {len(kx_out)} rows to {outfile}')

# --------------------------------------------------------
#   6/ OPTIONALLY PLOT THE EXTENDED k-SPACE MAPS
# --------------------------------------------------------

if do_plot:
    plot_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'plotKSpace.py')
    print(f'plotting {outfile}')
    subprocess.run([sys.executable, plot_script, outfile], check=False)

# Convert a photonicTest.py field-reflectivity pair (Magnitude + Phase CSVs, on
# a polar x azimuthal incidence-angle grid) into a long-format k-space table
# whose independent variables are the in-plane wavevector components kx, ky.
#
# This is the coordinate-adjustment step that feeds fitPhotonicTest.py, which
# fits analytic functions r_magnitude(kx, ky) and r_phase(kx, ky).
#
# Angle -> wavevector convention (S4, incidence medium = air, n = 1):
#   IncidenceAngles=(polar, azi) in degrees, so with k0 = 2*pi/lam,
#     k_parallel = k0 * sin(polar)
#     kx = k0 * sin(polar) * cos(azi)
#     ky = k0 * sin(polar) * sin(azi)
#     kz = k0 * cos(polar)   (propagation direction, not written out)
# (The polar angle rotates the k-vector about y, the azimuthal about z; see
# photonicTest.py and the S4 SetExcitationPlanewave reference.)
#
# Input CSV layout (from photonicTest.py):
#   - header row : azimuthal angles (deg), the COLUMNS
#   - first column of each row : polar angle (deg), the ROWS
#   - cells : |r| (Magnitude.csv) or arg(r) in degrees (Phase.csv)
#
# Output CSV: one row per grid cell, with a named header:
#   kx,ky,k_parallel,polar_deg,azi_deg,magnitude,phase_deg
#
# Usage:
#   python kSpaceConvert.py Magnitude.csv Phase.csv [output.csv]

import os
import subprocess
import sys

import numpy as np


# --------------------------------------------------------
#   PARAMETERS
# --------------------------------------------------------

lam = 1.064          # wavelength (um); k0 = 2*pi/lam, incidence medium is air

# "physical"   -> kx, ky in rad/um (k0 = 2*pi/lam)
# "normalized" -> kx, ky divided by k0, i.e. the direction cosines sin(t)cos(p)
K_UNITS = "physical"

UNWRAP_PHASE = True  # unwrap arg(r) along the polar axis before writing, so the
                     # phase is a smooth surface for fitting (no-op if it never
                     # wraps, as in the near-normal-incidence sweeps)


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


def read_grid(path):
    """Read one photonicTest CSV: header row = azimuthal angles, first column =
    polar angles, cells = values (NaN where blank)."""
    with open(path) as f:
        header = f.readline().strip().split(',')
    azis = np.array([float(v) for v in header[1:] if v.strip() != ''])

    polars = []
    values = []
    with open(path) as f:
        next(f)  # skip header row
        for line in f:
            cells = line.strip().split(',')
            if cells[0].strip() == '':
                continue
            polars.append(float(cells[0]))
            row = [float(v) if v.strip() != '' else np.nan
                   for v in cells[1:1 + len(azis)]]
            values.append(row)

    return azis, np.array(polars), np.array(values)


# --------------------------------------------------------
#   1/ READ THE MAGNITUDE / PHASE PAIR
# --------------------------------------------------------

# positional args: Magnitude.csv Phase.csv [output.csv]; optional plot=true
# flag anywhere on the command line auto-plots the k-space maps at the end
positional = []
do_plot = False
for a in sys.argv[1:]:
    if a.lower().startswith('plot='):
        do_plot = a.split('=', 1)[1].strip().lower() in ('true', '1', 'yes', 'on')
    else:
        positional.append(a)

if len(positional) < 2:
    sys.exit('usage: python kSpaceConvert.py Magnitude.csv Phase.csv '
             '[output.csv] [plot=true]')

mag_file, phase_file = positional[0], positional[1]

azis, polars, mag = read_grid(mag_file)
azis_p, polars_p, phase = read_grid(phase_file)

if (len(azis) != len(azis_p) or len(polars) != len(polars_p)
        or not np.allclose(azis, azis_p) or not np.allclose(polars, polars_p)):
    sys.exit('grid mismatch: the two CSVs cover different angle grids')

# unwrap the phase along polar (axis 0) so it is smooth for downstream fitting
if UNWRAP_PHASE:
    phase = np.rad2deg(np.unwrap(np.deg2rad(phase), axis=0))

# --------------------------------------------------------
#   2/ CONVERT ANGLES -> WAVEVECTOR
# --------------------------------------------------------

k0 = 2.0 * np.pi / lam                       # rad/um
k_prefactor = k0 if K_UNITS == "physical" else 1.0
if K_UNITS not in ("physical", "normalized"):
    sys.exit(f"K_UNITS must be 'physical' or 'normalized', got {K_UNITS!r}")

theta = np.deg2rad(polars)                    # polar, one per row
phi = np.deg2rad(azis)                        # azimuthal, one per column

# broadcast to the full (polar, azi) grid
kpar = k_prefactor * np.sin(theta)[:, None] * np.ones_like(phi)[None, :]
kx = kpar * np.cos(phi)[None, :]
ky = kpar * np.sin(phi)[None, :]

# --------------------------------------------------------
#   3/ WRITE THE LONG-FORMAT k-SPACE CSV
# --------------------------------------------------------

stem = mag_file.rsplit('.', 1)[0]
if stem.endswith('Magnitude'):
    stem = stem[:-len('Magnitude')]

if len(positional) > 2:
    outfile = unique_path(positional[2])
else:
    outfile = unique_path(stem + 'KSpace.csv')

os.makedirs(os.path.dirname(outfile) or '.', exist_ok=True)

with open(outfile, 'w') as f:
    f.write('kx,ky,k_parallel,polar_deg,azi_deg,magnitude,phase_deg\n')
    for i, p in enumerate(polars):
        for j, a in enumerate(azis):
            if np.isnan(mag[i, j]) or np.isnan(phase[i, j]):
                continue
            f.write(f'{kx[i, j]},{ky[i, j]},{kpar[i, j]},{p},{a},'
                    f'{mag[i, j]},{phase[i, j]}\n')

n_rows = int(np.sum(~(np.isnan(mag) | np.isnan(phase))))
print(f'read {mag.shape[0]} polar x {mag.shape[1]} azimuthal angles')
print(f'wavevector units : {K_UNITS} '
      f'({"rad/um, k0 = %.5f" % k0 if K_UNITS == "physical" else "kx/k0"})')
print(f'k_parallel range : {np.nanmin(kpar):.5f} - {np.nanmax(kpar):.5f}')
print(f'kx range         : {np.nanmin(kx):.5f} - {np.nanmax(kx):.5f}')
print(f'ky range         : {np.nanmin(ky):.5f} - {np.nanmax(ky):.5f}')
print(f'wrote {n_rows} rows to {outfile}')

# --------------------------------------------------------
#   4/ OPTIONALLY PLOT THE k-SPACE MAPS (|r| and phase vs kx, ky)
# --------------------------------------------------------

if do_plot:
    plot_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'plotKSpace.py')
    print(f'plotting {outfile}')
    subprocess.run([sys.executable, plot_script, outfile], check=False)

# --------------------------------------------------------
#   5/ AUTO-RUN THE FITTER ON THE k-SPACE CSV
# --------------------------------------------------------

fit_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'fitPhotonicTest.py')
print(f'fitting {outfile}')
subprocess.run([sys.executable, fit_script, outfile], check=False)

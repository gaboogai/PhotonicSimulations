# Plot a ...KSpace.csv (from kSpaceConvert.py) as two 2D maps over the in-plane
# wavevector plane:
#   left  : field reflectivity magnitude |r|   vs (kx, ky)
#   right : field reflectivity phase   arg(r)  vs (kx, ky)
#
# The CSV is long-format with a named header:
#   kx,ky,k_parallel,polar_deg,azi_deg,magnitude,phase_deg
# The (kx, ky) samples form a wedge (the polar x azimuthal sweep mapped into
# k-space, not a rectangular grid), so the maps are drawn with a triangulation
# (tripcolor). Both colour scales autoscale to their own data (dynamic range).
#
# Axis unit note: kx, ky are labelled rad/um, matching kSpaceConvert.py's
# default K_UNITS='physical'. If you converted with K_UNITS='normalized', the
# values are kx/k0 (dimensionless); edit K_AXIS_LABEL below.
#
# Usage:
#   python plotKSpace.py path/to/...KSpace.csv [output.png]

import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

K_AXIS_LABEL = ('kx (rad/um)', 'ky (rad/um)')


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


# --------------------------------------------------------
#   1/ CHOOSE THE OUTPUT TO PLOT
# --------------------------------------------------------

if len(sys.argv) < 2:
    sys.exit('usage: python plotKSpace.py path/to/...KSpace.csv [output.png]')

infile = sys.argv[1]
if len(sys.argv) > 2:
    outfile = unique_path(sys.argv[2])
else:
    outfile = unique_path(infile.rsplit('.', 1)[0] + '.png')

# --------------------------------------------------------
#   2/ READ AND DE-DUPLICATE THE POINTS
# --------------------------------------------------------

d = read_kspace(infile)
kx, ky = d['kx'], d['ky']
mag, phase = d['magnitude'], d['phase_deg']

# every azimuth collapses to the same point at normal incidence (kx=ky=0), so
# drop duplicate (kx, ky) points before triangulating
_, keep = np.unique(np.round(np.column_stack([kx, ky]), 12), axis=0,
                    return_index=True)
keep = np.sort(keep)
kx, ky, mag, phase = kx[keep], ky[keep], mag[keep], phase[keep]

# --------------------------------------------------------
#   3/ PLOT AND SAVE
# --------------------------------------------------------

fig, (ax_mag, ax_ph) = plt.subplots(1, 2, figsize=(12, 5))

mesh = ax_mag.tripcolor(kx, ky, mag, shading='gouraud', cmap='viridis')
cbar = fig.colorbar(mesh, ax=ax_mag)
cbar.set_label('|r|')
ax_mag.set_xlabel(K_AXIS_LABEL[0])
ax_mag.set_ylabel(K_AXIS_LABEL[1])
ax_mag.set_aspect('equal')
ax_mag.set_title('Field reflectivity magnitude')

mesh = ax_ph.tripcolor(kx, ky, phase, shading='gouraud', cmap='viridis')
cbar = fig.colorbar(mesh, ax=ax_ph)
cbar.set_label('arg(r) (deg)')
ax_ph.set_xlabel(K_AXIS_LABEL[0])
ax_ph.set_ylabel(K_AXIS_LABEL[1])
ax_ph.set_aspect('equal')
ax_ph.set_title('Field reflectivity phase')

fig.suptitle(os.path.basename(infile))
fig.tight_layout()
fig.savefig(outfile, dpi=150, bbox_inches='tight')
plt.close(fig)

print(f'read {len(d["kx"])} points ({len(kx)} unique) from {infile}')
print(f'|r| range   : {np.nanmin(mag):.6f} - {np.nanmax(mag):.6f}')
print(f'phase range : {np.nanmin(phase):.4f} - {np.nanmax(phase):.4f} deg')
print(f'saved {outfile}')

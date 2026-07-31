# Plot an angle-sweep field-reflectivity pair (the Magnitude + Phase CSVs
# written by photonicTest.py) as a 2x2 figure over the 2D
# (polar, azimuthal) incidence-angle grid:
#   top row    : |r| and arg(r) heatmaps, x = polar angle, y = azimuthal angle
#   bottom row : two 1D sections of |r| vs polar angle, taken at the first and
#                last azimuthal angles of the sweep
#
# The CSV layout (written by photonicTest.py) is:
#   - header row : azimuthal angles (deg), the COLUMNS
#   - first column of each row : polar angle (deg), the ROWS
#   - cells : |r| (Magnitude.csv) or arg(r) in degrees (Phase.csv)
#
# Both heatmaps autoscale their colour range to the data.
#
# Usage:
#   python plotPhotonicTest.py Magnitude.csv Phase.csv [output.png]

import os
import sys

import numpy as np

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
#   1/ CHOOSE THE OUTPUTS TO PLOT
# --------------------------------------------------------

if len(sys.argv) < 3:
    sys.exit('usage: python plotPhotonicTest.py Magnitude.csv Phase.csv '
             '[output.png]')

mag_file, phase_file = sys.argv[1], sys.argv[2]

if len(sys.argv) > 3:
    outfile = unique_path(sys.argv[3])
else:
    stem = mag_file.rsplit('.', 1)[0]
    if stem.endswith('Magnitude'):
        stem = stem[:-len('Magnitude')]
    outfile = unique_path(stem + 'MagnitudePhase.png')

# --------------------------------------------------------
#   2/ READ THE TWO MATRICES
# --------------------------------------------------------

azis, polars, mag = read_grid(mag_file)
azis_p, polars_p, phase = read_grid(phase_file)

if (len(azis) != len(azis_p) or len(polars) != len(polars_p)
        or not np.allclose(azis, azis_p) or not np.allclose(polars, polars_p)):
    sys.exit('grid mismatch: the two CSVs cover different angle grids')

# --------------------------------------------------------
#   3/ PLOT AND SAVE
# --------------------------------------------------------

fig, ((ax_mag, ax_phase), (ax_cut0, ax_cut1)) = plt.subplots(2, 2,
                                                             figsize=(11, 9))

# heatmaps: x = polar, y = azimuthal (matrices are (polar, azi) -> transpose)
mesh = ax_mag.pcolormesh(polars, azis, mag.T, shading='auto', cmap='viridis')
cbar = fig.colorbar(mesh, ax=ax_mag)
cbar.set_label('|r|')
ax_mag.set_xlabel('Polar angle (deg)')
ax_mag.set_ylabel('Azimuthal angle (deg)')
ax_mag.set_title('Field reflectivity magnitude')

mesh = ax_phase.pcolormesh(polars, azis, phase.T, shading='auto',
                           cmap='viridis')
cbar = fig.colorbar(mesh, ax=ax_phase)
cbar.set_label('arg(r) (deg)')
ax_phase.set_xlabel('Polar angle (deg)')
ax_phase.set_ylabel('Azimuthal angle (deg)')
ax_phase.set_title('Field reflectivity phase')

# 1D sections: |r| vs polar angle at the first and last azimuthal angles
ax_cut0.scatter(polars, mag[:, 0], s=12)
ax_cut0.set_xlabel('Polar angle (deg)')
ax_cut0.set_ylabel('|r|')
ax_cut0.set_title(f'|r| vs polar at azimuthal = {azis[0]:g} deg')

ax_cut1.scatter(polars, mag[:, 30], s=12)
ax_cut1.set_xlabel('Polar angle (deg)')
ax_cut1.set_ylabel('|r|')
ax_cut1.set_title(f'|r| vs polar at azimuthal = {azis[30]:g} deg')

fig.suptitle(os.path.basename(mag_file))
fig.tight_layout()
fig.savefig(outfile, dpi=150, bbox_inches='tight')
plt.close(fig)

print(f'read {mag.shape[0]} polar x {mag.shape[1]} azimuthal angles '
      f'from {mag_file}')
print(f'|r| range: {np.nanmin(mag):.4f} - {np.nanmax(mag):.4f}')
print(f'phase range: {np.nanmin(phase):.1f} - {np.nanmax(phase):.1f} deg')
print(f'saved {outfile}')

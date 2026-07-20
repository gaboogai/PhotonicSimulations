# Plot the reflectance from a single simbroad output (one cleanbroad<angle>.csv)
# as a 2D heatmap over the lattice / radius sweep grid, and save it as a PNG.
#
# The CSV layout (written by simbroad.py) is:
#   - header row : lattice parameters (um)
#   - first column of each row : hole radius (um)
#   - cells : reflectance R (blank where the holes overlap, lattice <= 2*radius)

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


# --------------------------------------------------------
#   1/ CHOOSE THE OUTPUT TO PLOT
# --------------------------------------------------------

infile = sys.argv[1] if len(sys.argv) > 1 else 'cleanbroad0.csv'
outfile = unique_path(infile.rsplit('.', 1)[0] + '.png')

# --------------------------------------------------------
#   2/ READ THE MATRIX
# --------------------------------------------------------

# header row gives the lattice values (skip the leading blank corner cell)
with open(infile) as f:
    header = f.readline().strip().split(',')
lattices = np.array([float(v) for v in header[1:] if v.strip() != ''])

# remaining rows: first column is the radius, the rest are reflectance values
radii = []
reflectance = []
with open(infile) as f:
    next(f)  # skip header row
    for line in f:
        cells = line.strip().split(',')
        if cells[0].strip() == '':
            continue
        radii.append(float(cells[0]))
        row = [float(v) if v.strip() != '' else np.nan
               for v in cells[1:1 + len(lattices)]]
        reflectance.append(row)

radii = np.array(radii)
reflectance = np.array(reflectance)  # shape (n_radii, n_lattices)

# --------------------------------------------------------
#   3/ PLOT AND SAVE
# --------------------------------------------------------

plt.figure()
mesh = plt.pcolormesh(lattices, radii, reflectance,
                      shading='auto', cmap='viridis')
cbar = plt.colorbar(mesh)
cbar.set_label('Reflectance')
plt.xlabel('Lattice parameter a (um)')
plt.ylabel('Hole radius (um)')
plt.title(f'Reflectance - {infile}')
plt.autoscale(tight=True)
plt.savefig(outfile, transparent=False, dpi=150, bbox_inches='tight')
plt.close()

print(f'read {reflectance.shape[0]} radii x {reflectance.shape[1]} lattices '
      f'from {infile}')
print(f'reflectance range: {np.nanmin(reflectance):.4f} - '
      f'{np.nanmax(reflectance):.4f}')
print(f'saved {outfile}')

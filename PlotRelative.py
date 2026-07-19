# Plot the reflectance from a single HexRelative.py output as a 2D heatmap.
#
# The HexRelative CSV layout is:
#   - header row : radius ratios (r / a), the COLUMNS
#   - first column of each row : lattice constant a (um), the ROWS
#   - cells : reflectance R (blank where the holes overlap)
#
# Usage:
#   python PlotRelative.py [path/to/HexRelative...csv]
# With no argument it plots the most recent HexRelative*.csv in RawData/.

import glob
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

infile = sys.argv[1]

outfile = unique_path(infile.rsplit('.', 1)[0] + '.png')

# --------------------------------------------------------
#   2/ READ THE MATRIX
# --------------------------------------------------------

# header row gives the radius ratios (skip the leading blank corner cell)
with open(infile) as f:
    header = f.readline().strip().split(',')
ratios = np.array([float(v) for v in header[1:] if v.strip() != ''])

# remaining rows: first column is the lattice constant, the rest reflectance
lattices = []
reflectance = []
with open(infile) as f:
    next(f)  # skip header row
    for line in f:
        cells = line.strip().split(',')
        if cells[0].strip() == '':
            continue
        lattices.append(float(cells[0]))
        row = [float(v) if v.strip() != '' else np.nan
               for v in cells[1:1 + len(ratios)]]
        reflectance.append(row)

lattices = np.array(lattices)
reflectance = np.array(reflectance)  # shape (n_lattices, n_ratios)

# --------------------------------------------------------
#   3/ PLOT AND SAVE
# --------------------------------------------------------

plt.figure()
mesh = plt.pcolormesh(ratios, lattices, reflectance,
                      shading='auto', cmap='viridis')
cbar = plt.colorbar(mesh)
cbar.set_label('Reflectance')
plt.xlabel('Hole radius (fraction of lattice, r/a)')
plt.ylabel('Lattice parameter a (um)')
plt.title(f'Reflectance - {os.path.basename(infile)}')
plt.autoscale(tight=True)
plt.savefig(outfile, transparent=False, dpi=150, bbox_inches='tight')
plt.close()

print(f'read {reflectance.shape[0]} lattices x {reflectance.shape[1]} ratios '
      f'from {infile}')
print(f'reflectance range: {np.nanmin(reflectance):.4f} - '
      f'{np.nanmax(reflectance):.4f}')
print(f'saved {outfile}')

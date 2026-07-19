# Feasibility plot: compare two HexRelative.py outputs (typically the same
# lattice/ratio grid at two different incidence angles) and show how much the
# reflectivity changes with angle - but only in the regions that are actually
# useful, i.e. where the FIRST file already reflects strongly (R > 0.15).
#
# For every (radius ratio, lattice constant) cell it plots
#     R(file_a) - R(file_b)
# masked to cells where R(file_a) > THRESHOLD; all other cells are left blank.
# A small difference in a high-reflectivity region means an angle-robust design.
#
# Each input CSV has the HexRelative layout:
#   - header row : radius ratios (r / a), the COLUMNS
#   - first column of each row : lattice constant a (um), the ROWS
#   - cells : reflectance R (blank where the holes overlap)
#
# Usage:
#   python FeasibilityPlot.py file_a.csv file_b.csv [output.png]
# file_a is the one the R > 0.15 threshold is applied to.

import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

THRESHOLD = 0.15   # only show cells where R(file_a) exceeds this


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


def read_reflectance(infile):
    """Read a HexRelative CSV into (ratios, lattices, reflectance matrix)."""
    with open(infile) as f:
        header = f.readline().strip().split(',')
        ratios = np.array([float(v) for v in header[1:] if v.strip() != ''])

        lattices = []
        reflectance = []
        for line in f:
            cells = line.strip().split(',')
            if cells[0].strip() == '':
                continue
            lattices.append(float(cells[0]))
            row = [float(v) if v.strip() != '' else np.nan
                   for v in cells[1:1 + len(ratios)]]
            reflectance.append(row)

    return ratios, np.array(lattices), np.array(reflectance)


# --------------------------------------------------------
#   1/ CHOOSE THE TWO OUTPUTS TO COMPARE
# --------------------------------------------------------

if len(sys.argv) < 3:
    sys.exit('usage: python FeasibilityPlot.py file_a.csv file_b.csv '
             '[output.png]')

file_a = sys.argv[1]
file_b = sys.argv[2]
outfile = sys.argv[3] if len(sys.argv) > 3 else 'feasibility_diff.png'
outfile = unique_path(outfile)

# --------------------------------------------------------
#   2/ READ BOTH MATRICES
# --------------------------------------------------------

ratios_a, lattices_a, R_a = read_reflectance(file_a)
ratios_b, lattices_b, R_b = read_reflectance(file_b)

if R_a.shape != R_b.shape:
    sys.exit(f'grid mismatch: {file_a} is {R_a.shape}, {file_b} is {R_b.shape}')
if not (np.allclose(ratios_a, ratios_b) and np.allclose(lattices_a, lattices_b)):
    sys.exit('the two outputs are on different ratio/lattice grids')

# difference, masked to cells where file_a reflects strongly enough
diff = R_a - R_b
mask = ~(R_a > THRESHOLD)           # True where we DON'T want to show a value
diff_masked = np.where(mask, np.nan, diff)

n_feasible = int(np.sum(R_a > THRESHOLD))

# --------------------------------------------------------
#   3/ PLOT AND SAVE
# --------------------------------------------------------

# symmetric colour scale about zero (diverging map), robust to an all-masked grid
finite = diff_masked[np.isfinite(diff_masked)]
vmax = np.max(np.abs(finite)) if finite.size else 1.0
if vmax == 0:
    vmax = 1.0

plt.figure()
mesh = plt.pcolormesh(ratios_a, lattices_a, diff_masked,
                      shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
cbar = plt.colorbar(mesh)
cbar.set_label(f'Reflectance difference (R_a - R_b), R_a > {THRESHOLD:g} only')
plt.xlabel('Hole radius (fraction of lattice, r/a)')
plt.ylabel('Lattice parameter a (um)')
plt.title(f'{os.path.basename(file_a)} - {os.path.basename(file_b)}')
plt.autoscale(tight=True)
plt.savefig(outfile, transparent=False, dpi=150, bbox_inches='tight')
plt.close()

print(f'compared {R_a.shape[0]} lattices x {R_a.shape[1]} ratios')
print(f'cells with R(file_a) > {THRESHOLD:g}: {n_feasible} of {R_a.size}')
if finite.size:
    print(f'difference range (over shown cells): '
          f'{np.nanmin(diff_masked):.4f} - {np.nanmax(diff_masked):.4f}')
else:
    print('no cells exceed the threshold - nothing plotted')
print(f'saved {outfile}')

# Plot the difference in reflectance between two simbroad outputs
# (two cleanbroad<angle>.csv matrices) as a 2D heatmap, and save it as a PNG.
#
# The CSV layout (written by simbroad.py) is:
#   - header row : lattice parameters (um)
#   - first column of each row : hole radius (um)
#   - cells : reflectance R (blank where the holes overlap, lattice <= 2*radius)
#
# Usage:
#   python plotReflectanceDiff.py file_a.csv file_b.csv [output.png]
# The plotted quantity is R(file_a) - R(file_b).

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


def read_reflectance(infile):
    """Read a cleanbroad CSV into (lattices, radii, reflectance matrix)."""
    with open(infile) as f:
        header = f.readline().strip().split(',')
        lattices = np.array([float(v) for v in header[1:] if v.strip() != ''])

        radii = []
        reflectance = []
        for line in f:
            cells = line.strip().split(',')
            if cells[0].strip() == '':
                continue
            radii.append(float(cells[0]))
            row = [float(v) if v.strip() != '' else np.nan
                   for v in cells[1:1 + len(lattices)]]
            reflectance.append(row)

    return lattices, np.array(radii), np.array(reflectance)


# --------------------------------------------------------
#   1/ CHOOSE THE TWO OUTPUTS TO COMPARE
# --------------------------------------------------------

if len(sys.argv) < 3:
    sys.exit('usage: python plotReflectanceDiff.py file_a.csv file_b.csv '
             '[output.png]')

file_a = sys.argv[1]
file_b = sys.argv[2]
outfile = sys.argv[3] if len(sys.argv) > 3 else 'reflectance_diff.png'

# --------------------------------------------------------
#   2/ READ BOTH MATRICES
# --------------------------------------------------------

lattices_a, radii_a, R_a = read_reflectance(file_a)
lattices_b, radii_b, R_b = read_reflectance(file_b)

if R_a.shape != R_b.shape:
    sys.exit(f'grid mismatch: {file_a} is {R_a.shape}, '
             f'{file_b} is {R_b.shape}')
if not (np.allclose(lattices_a, lattices_b) and np.allclose(radii_a, radii_b)):
    sys.exit('the two outputs are on different lattice/radius grids')

diff = R_a - R_b

# --------------------------------------------------------
#   3/ PLOT AND SAVE
# --------------------------------------------------------

# symmetric colour scale about zero so the diverging map is centred
vmax = np.nanmax(np.abs(diff))

plt.figure()
mesh = plt.pcolormesh(lattices_a, radii_a, diff,
                      shading='auto', cmap='RdBu_r',
                      vmin=-vmax, vmax=vmax)
cbar = plt.colorbar(mesh)
cbar.set_label('Reflectance difference')
plt.xlabel('Lattice parameter (um)')
plt.ylabel('Hole radius (um)')
plt.title(f'Reflectance: {file_a} - {file_b}')
plt.autoscale(tight=True)
outfile = unique_path(outfile)
plt.savefig(outfile, transparent=False, dpi=150, bbox_inches='tight')
plt.close()

print(f'compared {R_a.shape[0]} radii x {R_a.shape[1]} lattices')
print(f'difference range: {np.nanmin(diff):.4f} - {np.nanmax(diff):.4f}')
print(f'saved {outfile}')

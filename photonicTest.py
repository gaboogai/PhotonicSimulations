# Angle sweep of the complex FIELD reflectivity for a FIXED lattice/radius point.
#
# FieldReflectanceSweep.py sweeps the (lattice, radius-ratio) grid at fixed
# incidence angles; this script is the complementary cut: the lattice constant
# and hole radius are held fixed and the two incidence angles are swept over a
# 2D (polar x azimuthal) grid. Per S4's Python API, IncidenceAngles=(polar, azi)
# is in degrees: the polar angle rotates the incident frame about y, and the
# azimuthal angle then rotates it about z. (Older scripts' comments say the
# azimuthal rotation is about x; the S4 API reference says z.)
#
# The per-cell field reflectivity is computed by FieldReflectance.field_reflectance
# (imported, so the physics stays identical to the single-point script). Each cell
# yields two numbers - magnitude |r| and phase arg(r) - so TWO CSVs are written
# in the usual matrix layout (header row = azimuthal angles, first column =
# polar angles):
#   ...Magnitude...csv : |r|
#   ...Phase...csv     : arg(r) in degrees
# Both are plotted together as a two-panel heatmap by plotPhotonicTest.py at
# the end.
#
# Usage:
#   python photonicTest.py

import os
import subprocess
import sys
import time

import numpy as np

import FieldReflectance as FR   # reuse the verified single-point field solver


def unique_path(path):
    """Return path, or path with an incrementing _N suffix if it already exists,
    so existing output files are never overwritten."""
    base, ext = os.path.splitext(path)
    candidate = path
    n = 1
    while os.path.exists(candidate):
        candidate = f"{base}_{n}{ext}"
        n += 1
    return candidate


# --------------------------------------------------------
#   1/ SET THE PARAMETERS
# --------------------------------------------------------

lam = 1.064          # wavelength (um)

shape = "hex"        # "hex" or "square" - selectable geometry

# the fixed lattice/radius point to evaluate (um)
lattice = 1.1125     # lattice constant a
radius = 0.37        # hole radius r

thickness = 0.1      # slab thickness (um)
n_orders = 20        # number of Fourier orders

# polar angle range (degrees)
polar_start = 0.0
polar_stop = 1.0
n_polar = 50
polars = np.linspace(polar_start, polar_stop, n_polar + 1)

# azimuthal angle range (degrees)
azi_start = 0.0
azi_stop = 60.0
n_azi = 60
azis = np.linspace(azi_start, azi_stop, n_azi + 1)

# push the chosen physics settings into the imported single-point solver so both
# scripts compute identically
FR.shape = shape
FR.lam = lam
FR.thickness = thickness
FR.n_orders = n_orders


# --------------------------------------------------------
#   2/ RUN AND WRITE TWO CSVs (magnitude, phase)
#      (polar-outer rows / azimuthal-inner columns)
# --------------------------------------------------------

if __name__ == "__main__":

    if lattice <= 2 * radius or lattice < 0:
        sys.exit(f"holes overlap: lattice a = {lattice:g} um must exceed "
                 f"2 * radius = {2 * radius:g} um")

    tic = time.time()

    stem = (f"AngleSweep/{shape}FieldAngleL:{lattice:g}R:{radius:g}"
            f"Pstart:{polar_start:g}Pstop:{polar_stop:g}"
            f"Astart:{azi_start:g}Astop:{azi_stop:g}")
    mag_file = unique_path(f"{stem}Magnitude.csv")
    phase_file = unique_path(f"{stem}Phase.csv")
    os.makedirs(os.path.dirname(mag_file), exist_ok=True)

    with open(mag_file, "w") as fm, open(phase_file, "w") as fp:

        # header row: azimuthal angles (deg)
        fm.write(" ,")
        fp.write(" ,")
        for azi in azis:
            fm.write(f"{azi},")
            fp.write(f"{azi},")
        fm.write("\n")
        fp.write("\n")

        # one row per polar angle
        for polar in polars:
            fm.write(f"{polar},")
            fp.write(f"{polar},")

            for azi in azis:
                d = FR.field_reflectance(polar, azi, lattice, radius)
                fm.write(f"{d['mag']},")
                fp.write(f"{d['phase_deg']},")

            fm.write("\n")
            fp.write("\n")

            done = 100 * (polar - polar_start + (polars[1] - polars[0])) \
                / (polar_stop - polar_start)
            print(f"polar:{polar:g} - completion: {done:.1f}%")

    toc = time.time()
    print("Total time: " + str((toc - tic) / 60))

    # plot the magnitude/phase pair as a two-panel heatmap
    plot_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "plotPhotonicTest.py")
    print(f"plotting {mag_file} + {phase_file}")
    subprocess.run([sys.executable, plot_script, mag_file, phase_file],
                   check=False)

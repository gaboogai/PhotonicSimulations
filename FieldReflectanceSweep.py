# Sweep version of FieldReflectance.py: complex FIELD reflectivity over a
# (lattice, radius-ratio) grid, for each (polar, azimuthal) incidence angle.
#
# Like the *Relative* scripts, the hole radius is defined RELATIVE to the lattice
# (r = ratio * a) and the loops are lattice-outer / ratio-inner. Geometry is
# selectable (hex or square) via `shape`, and both incidence angles are swept.
#
# The per-cell field reflectivity is computed by FieldReflectance.field_reflectance
# (imported, so the physics stays identical to the single-point script). Each cell
# yields two numbers - magnitude |r| and phase arg(r) - so for every (polar, azi)
# angle TWO CSVs are written in the HexRelative layout that PlotRelative.py reads:
#   ...Magnitude...csv : |r|
#   ...Phase...csv     : arg(r) in degrees
# Both are auto-plotted with PlotRelative.py at the end.
#
# Usage:
#   python FieldReflectanceSweep.py

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

# incidence angles (degrees); polar = rot about y, azimuthal = rot about x
angles_polar = [0.313]
angles_azi = [0,30.0]

thickness = 0.03     # slab thickness (um)
n_orders = 20        # number of Fourier orders

# lattice parameter range (um)
lattice_start = 1.05
lattice_stop = 1.3
n_lattices = 100
lattices = np.linspace(lattice_start, lattice_stop, n_lattices + 1)

# radius range, as a fraction of the current lattice constant (r = ratio * a)
ratio_start = 0.075
ratio_stop = 0.5
n_radius = 100
ratios = np.linspace(ratio_start, ratio_stop, n_radius + 1)

# push the chosen physics settings into the imported single-point solver so both
# scripts compute identically
FR.shape = shape
FR.lam = lam
FR.thickness = thickness
FR.n_orders = n_orders


# --------------------------------------------------------
#   2/ RUN AND WRITE TWO CSVs (magnitude, phase) PER ANGLE
#      (lattice-outer / ratio-inner: radius = ratio * lattice)
# --------------------------------------------------------

if __name__ == "__main__":

    tic = time.time()

    outputs = []

    for angle_azi in angles_azi:
        for angle_polar in angles_polar:

            stem = (f"ThicknessTestHex/{shape}FieldPolar:{angle_polar:g}Azi:{angle_azi:g}"
                    f"Lstart:{lattice_start:g}Lstop:{lattice_stop:g}"
                    f"Rstart{ratio_start:g}Rstop{ratio_stop:g}Thickness{thickness:g}")
            mag_file = unique_path(f"{stem}Magnitude.csv")
            phase_file = unique_path(f"{stem}Phase.csv")
            os.makedirs(os.path.dirname(mag_file), exist_ok=True)
            outputs.extend((mag_file, phase_file))

            with open(mag_file, "w") as fm, open(phase_file, "w") as fp:

                # header row: radius ratios (r / a)
                fm.write(" ,")
                fp.write(" ,")
                for ratio in ratios:
                    fm.write(f"{ratio},")
                    fp.write(f"{ratio},")
                fm.write("\n")
                fp.write("\n")

                # one row per lattice constant
                for lattice in lattices:
                    fm.write(f"{lattice},")
                    fp.write(f"{lattice},")

                    for ratio in ratios:
                        radius = ratio * lattice
                        if lattice <= 2 * radius or lattice < 0:
                            fm.write(" ,")       # holes overlap -> blank, as in Lua
                            fp.write(" ,")
                        else:
                            d = FR.field_reflectance(angle_polar, angle_azi,
                                                     lattice, radius)
                            fm.write(f"{d['mag']},")
                            fp.write(f"{d['phase_deg']},")

                    fm.write("\n")
                    fp.write("\n")

                    done = 100 * (lattice - lattice_start + (lattices[1] - lattices[0])) \
                        / (lattice_stop - lattice_start)
                    print(f"polar:{angle_polar:g} azi:{angle_azi:g} - completion: {done:.1f}%")

    toc = time.time()
    print("Total time: " + str((toc - tic) / 60))

    # automatically plot each CSV that was written (magnitude and phase heatmaps)
    plot_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "PlotRelative.py")
    for output in outputs:
        print(f"plotting {output}")
        subprocess.run([sys.executable, plot_script, output], check=False)

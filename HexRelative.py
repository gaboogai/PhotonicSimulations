# Variant of simbroad.py with the hole radius defined RELATIVE to the lattice.
#
# Sweeps a thin SiN slab patterned with a hexagonal lattice of circular vacuum
# holes (suspended in air, p-polarised plane wave). Instead of sweeping an
# absolute radius range, the radius is taken as a fixed set of fractions of the
# current lattice constant (0.35 - 0.49 of a). Because the radius now depends on
# the lattice, the loops in section 3 are ordered lattice-outer / ratio-inner.
#
# For every angle it writes one CSV whose COLUMNS are the radius ratios (r/a)
# and whose ROWS are the lattice constants, each cell the reflectance.

import os
import subprocess
import sys
import time

import numpy as np

import S4  # native S4 Python extension, used by reflectance()


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

# angle of incidence range (degrees)
angles_polar = [0] #pretty sure this is incident, rot about y 0.313
angles_azi = [0.0,15.0,30.0,45.0,60.0] #,75.0,90.0, 105.0,120.0,135.0,150.0,165.0,180.0] #azimuthal, rot about x

thickness = 0.03     # slab thickness (um)
shape = "hex"        # pattern type (hexagonal lattice)
material_name = "SiN"

# lattice parameter range (um)
lattice_start = 1.188
lattice_stop = 1.228
n_lattices = 100
lattices = np.linspace(lattice_start, lattice_stop, n_lattices + 1)

# radius range, as a fraction of the current lattice constant (r = ratio * a)
ratio_start = 0.075
ratio_stop = 0.5
n_radius = 100
ratios = np.linspace(ratio_start, ratio_stop, n_radius + 1)


n_orders = 20        # number of Fourier orders (SetNumG in the Lua script)


class SellmeierSiN:
    """SiN with the same Sellmeier dispersion used in simbroad copy.lua
    (solcore's database Si3N4 gives n = 2.0013 at 1064 nm vs 2.0114 here)."""

    material_string = "SiN"

    def __str__(self):
        return self.material_string

    def n(self, wl):
        lam2 = (np.asarray(wl) * 1e6) ** 2  # wl is in m, formula is in um
        eps = (1 + 3.0249 * lam2 / (lam2 - 0.1353406**2)
               + 40314 * lam2 / (lam2 - 1239.842**2))
        return np.sqrt(eps)

    def k(self, wl):
        return np.zeros_like(np.asarray(wl, dtype=float))


# --------------------------------------------------------
#   2/ SWEEP RANGES
# --------------------------------------------------------


SiN = SellmeierSiN()
wavelengths = np.array([lam]) * 1e-6


def reflectance(angle_polar, angle_azi, lattice, radius):
    """Reflectance for one (angle, lattice, radius), driven purely through the
    native S4 Python API (no solcore/rayflare). Lengths in nm, frequency =
    1/wavelength, reflectance = Re(-back/forward) of the incidence-medium flux.
    """
    wl_nm = lam * 1e3                       # wavelength in nm
    a_nm = lattice * 1e3                    # lattice constant in nm

    # hexagonal lattice
    S = S4.New(Lattice=((a_nm, 0.0),
                        (np.cos(np.pi / 3) * a_nm, np.sin(np.pi / 3) * a_nm)),
               NumBasis=n_orders)

    # same solver options RayFlare applies to its S4 simulations
    S.SetOptions(LatticeTruncation="Circular", DiscretizedEpsilon=False,
                 DiscretizationResolution=8, PolarizationDecomposition=False,
                 PolarizationBasis="Default", LanczosSmoothing=False,
                 SubpixelSmoothing=False, ConserveMemory=False,
                 WeismannFormulation=False, Verbosity=0)

    # permittivities: air = 1, SiN from the same Sellmeier dispersion (eps = n^2)
    S.SetMaterial(Name="Air", Epsilon=1.0 + 0j)
    S.SetMaterial(Name="SiN", Epsilon=complex(SiN.n(lam * 1e-6)) ** 2)

    # incidence medium / patterned slab / transmission medium; semi-infinite
    # media have zero thickness in S4, slab thickness in nm
    S.AddLayer(Name="incidence", Thickness=0.0, Material="Air")
    S.AddLayer(Name="slab", Thickness=thickness * 1e3, Material="SiN")
    S.AddLayer(Name="transmission", Thickness=0.0, Material="Air")

    # circular vacuum hole in the slab
    S.SetRegionCircle(Layer="slab", Material="Air", Center=(0.0, 0.0),
                      Radius=radius * 1e3)

    # p-polarised plane wave; IncidenceAngles is the same (angle, 0) tuple
    # RayFlare passes through from (theta_in, phi_in)
    S.SetExcitationPlanewave(IncidenceAngles=(angle_polar, angle_azi),
                             sAmplitude=0.0, pAmplitude=1.0, Order=0)
    S.SetFrequency(1.0 / wl_nm)

    forw, back = S.GetPowerFlux(Layer="incidence")
    return float((-back / forw).real)


# --------------------------------------------------------
#   3/ RUN AND WRITE ONE CSV PER ANGLE
#   (lattice-outer / ratio-inner: radius = ratio * lattice)
# --------------------------------------------------------

if __name__ == "__main__":

    tic = time.time()

    outputs = []

    for angle_azi in angles_azi:
        for angle_polar in angles_polar:

            filename = unique_path(
                f"Two-Angle-Hex30nm/HexRelativePolar:{angle_polar:g}Azi:{angle_azi:g}"
                f"Lstart:{lattice_start:g}Lstop:{lattice_stop:g}"
                f"Rstart{ratio_start:g}Rstop{ratio_stop:g}.csv")
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            outputs.append(filename)
            with open(filename, "w") as f:

                # header row: radius ratios (r / a)
                f.write(" ,")
                for ratio in ratios:
                    f.write(f"{ratio},")
                f.write("\n")

                # one row per lattice constant
                for lattice in lattices:
                    f.write(f"{lattice},")

                    for ratio in ratios:
                        radius = ratio * lattice
                        if lattice <= 2 * radius or lattice < 0:
                            f.write(" ,")        # holes overlap -> blank, as in Lua
                        else:
                            f.write(f"{reflectance(angle_polar,angle_azi , lattice, radius)},")

                    f.write("\n")

                    done = 100 * (lattice - lattice_start + (lattices[1] - lattices[0])) \
                        / (lattice_stop - lattice_start)
                    print(f"polar:{angle_polar:g} azi:{angle_azi:g} - completion: {done:.1f}%")

    toc = time.time()
    print("Total time: " + str((toc - tic) / 60))

    # automatically plot each CSV that was written
    plot_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "PlotRelative.py")
    for output in outputs:
        print(f"plotting {output}")
        subprocess.run([sys.executable, plot_script, output], check=False)

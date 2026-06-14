# Full recreation of "simbroad copy.lua" using the python S4 extension via
# RayFlare (see Example.md).
#
# Sweeps a thin SiN slab patterned with a hexagonal lattice of circular vacuum
# holes (suspended in air, p-polarised plane wave) over the full radius and
# lattice ranges, for each angle of incidence. For every angle it writes
# cleanbroad<angle>.csv: a matrix whose columns are lattice parameters and whose
# rows are hole radii, with each cell the reflectance. Cells where the holes
# would overlap (lattice <= 2*radius) are left blank, exactly as in the Lua.

import time

import numpy as np
from solcore import material
from solcore.structure import Layer
from solcore.solar_cell import SolarCell

from rayflare.rigorous_coupled_wave_analysis import rcwa_structure
from rayflare.options import default_options

import S4  # native S4 Python extension, used by S4reflectance()

# --------------------------------------------------------
#   1/ SET THE PARAMETERS
# --------------------------------------------------------

lam = 1.064          # wavelength (um)

# angle of incidence range (degrees)
angle_start = 4.0
angle_stop = 5.0
n_angles = 1
angles = np.linspace(angle_start, angle_stop, n_angles + 1)

thickness = 0.02     # slab thickness (um)
shape = "hex"        # pattern type (hexagonal lattice)
material_name = "SiN"

# lattice parameter range (um)
lattice_start = 1.21
lattice_stop = 1.23
n_lattices = 100
lattices = np.linspace(lattice_start, lattice_stop, n_lattices + 1)

# radius range (um)
radius_start = 0.2
radius_stop = 0.6
n_radius = 100
radii = np.linspace(radius_start, radius_stop, n_radius + 1)


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
Air = material("Air")()
wavelengths = np.array([lam]) * 1e-6


def reflectance(angle, lattice, radius):
    """Reflectance for one (angle, lattice, radius), matching the Lua per-point
    S4 simulation. RayFlare uses m for wavelengths and nm for in-plane sizes."""
    a_nm = lattice * 1e3
    size = ((a_nm, 0),
            (np.cos(np.pi / 3) * a_nm, np.sin(np.pi / 3) * a_nm))

    options = default_options()
    options.wavelength = wavelengths
    options.orders = n_orders
    options.RCWA_method = "S4"
    options.pol = "p"                     # {0,0} s / {1,0} p excitation
    options.theta_in = np.deg2rad(angle)
    options.phi_in = 0.0

    geometry = [{"type": "circle", "mat": Air, "center": (0, 0),
                 "radius": radius * 1e3}]
    structure = SolarCell([Layer(material=SiN, width=thickness * 1e-6,
                                 geometry=geometry)])
    S4_setup = rcwa_structure(structure, size=size, options=options,
                              incidence=Air, transmission=Air)
    return S4_setup.calculate(options)["R"][0]


def S4reflectance(angle, lattice, radius):
    """Identical calculation to reflectance(), but driven purely through the
    native S4 Python API (no solcore/rayflare). It mirrors exactly how RayFlare
    sets up S4: lengths in nm, frequency = 1/wavelength, and the reflectance
    taken as Re(-back/forward) of the incidence-medium power flux.

    Not used by the sweep below (which calls reflectance()); provided as a
    pure-S4 equivalent.
    """
    wl_nm = lam * 1e3                       # wavelength in nm
    a_nm = lattice * 1e3                    # lattice constant in nm

    # hexagonal lattice, same vectors as reflectance()
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
    S.SetExcitationPlanewave(IncidenceAngles=(angle, 0.0),
                             sAmplitude=0.0, pAmplitude=1.0, Order=0)
    S.SetFrequency(1.0 / wl_nm)

    forw, back = S.GetPowerFlux(Layer="incidence")
    return float((-back / forw).real)


# --------------------------------------------------------
#   3/ RUN AND WRITE ONE CSV PER ANGLE
# --------------------------------------------------------

if __name__ == "__main__":

    tic = time.time()

    for angle in angles:

        filename = f"cleanbroad{angle:g}.csv"
        with open(filename, "w") as f:

            # header row: lattice parameters
            f.write(" ,")
            for lattice in lattices:
                f.write(f"{lattice},")
            f.write("\n")

            # one row per radius
            for radius in radii:
                f.write(f"{radius},")

                for lattice in lattices:
                    if lattice <= 2 * radius or lattice < 0:
                        f.write(" ,")        # holes overlap -> blank, as in Lua
                    else:
                        f.write(f"{reflectance(angle, lattice, radius)},")

                f.write("\n")

                done = 100 * (radius - radius_start + (radii[1] - radii[0])) \
                    / (radius_stop - radius_start)
                print(f"angle {angle:g} deg - completion: {done:.1f}%")

    toc = time.time()
    print("Total time: " + str((toc - tic) / 60))

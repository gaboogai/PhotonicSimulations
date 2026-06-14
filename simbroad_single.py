# Recreation of "simbroad copy.lua" for a single radius, lattice and angle,
# using the python S4 extension through RayFlare (see Example.md).
#
# Structure: a thin SiN slab patterned with a hexagonal lattice of circular
# vacuum holes, suspended in air, illuminated with a p-polarised plane wave.
# Output: the reflectance (the Lua script's back_r from GetPowerFlux('AirAbove'),
# which is negative there because backward flux points away from the source),
# written to cleanbroad<angle>.csv in the same header-row/header-column layout.

import time

import numpy as np
from solcore import material
from solcore.structure import Layer
from solcore.solar_cell import SolarCell

from rayflare.rigorous_coupled_wave_analysis import rcwa_structure
from rayflare.options import default_options

# --------------------------------------------------------
#   1/ SET THE PARAMETERS
# --------------------------------------------------------

lam = 1.064      # wavelength (um)
angle = 0.0      # angle of incidence (degrees)
thickness = 0.02 # slab thickness (um)
lattice = 1.21   # hexagonal lattice parameter (um)
radius = 0.2     # hole radius (um)
n_orders = 20    # number of Fourier orders (SetNumG in the Lua script)


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
#   2/ REFLECTIVITY CALCULATION
# --------------------------------------------------------
    
tic = time.time()

SiN = SellmeierSiN()
Air = material("Air")()

# RayFlare works in m for wavelengths and nm for lattice/shape dimensions
wavelengths = np.array([lam]) * 1e-6

# hexagonal lattice, same vectors as S:SetLattice in the Lua script
a_nm = lattice * 1e3
size = ((a_nm, 0), (np.cos(np.pi / 3) * a_nm, np.sin(np.pi / 3) * a_nm))

options = default_options()
options.wavelength = wavelengths
options.orders = n_orders
options.RCWA_method = "S4"
options.pol = "p"                      # {0,0} s / {1,0} p excitation in the Lua
options.theta_in = np.deg2rad(angle)
options.phi_in = 0.0

if lattice <= 2 * radius:
    raise ValueError("holes overlap: lattice <= 2*radius")

# SiN slab with a circular vacuum hole at the centre of the unit cell
geometry = [{"type": "circle", "mat": Air, "center": (0, 0),
             "radius": radius * 1e3}]

structure = SolarCell([Layer(material=SiN, width=thickness * 1e-6,
                             geometry=geometry)])

S4_setup = rcwa_structure(structure, size=size, options=options,
                          incidence=Air, transmission=Air)
RAT = S4_setup.calculate(options)

R = RAT["R"][0]
T = RAT["T"][0]
A = RAT["A_per_layer"][0, 0]

# --------------------------------------------------------
#   3/ OUTPUT
# --------------------------------------------------------

filename = f"cleanbroad{angle:g}.csv"
with open(filename, "w") as f:
    f.write(f" ,{lattice},\n")
    f.write(f"{radius},{R},\n")

print(f"lambda = {lam} um, angle = {angle} deg, thickness = {thickness} um, "
      f"lattice = {lattice} um, radius = {radius} um")
print(f"R = {R:.6f}, T = {T:.6f}, A = {A:.3e}")
print(f"written to {filename}")

toc = time.time()
print("Total time: " + str((toc - tic) / 60))

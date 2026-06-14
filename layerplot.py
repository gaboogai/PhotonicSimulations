# Plot the dielectric constant of the patterned SiN slab for a specified
# lattice parameter and hole radius, using only the native S4 Python API, and
# save it as a PNG.
#
# The slab is the same structure used in simbroad.py: a SiN layer with a
# central circular vacuum hole on a hexagonal lattice. The dielectric is sampled
# with S4's GetEpsilon over the unit cell at mid-slab depth; as with RayFlare's
# get_fourier_epsilon this is the epsilon reconstructed from the truncated
# G-vector set, so it looks 'fuzzy' (more Fourier orders sharpen it).
#
# The plot can be tiled over several unit cells (cells x cells) to show the
# periodic hexagonal arrangement of holes; GetEpsilon is periodic over the
# lattice, so a wider sampling window tiles automatically.
#
# Usage:
#   python layerplot.py [lattice_um] [radius_um] [cells] [output.png]

import sys

import numpy as np
import S4

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------
#   1/ SET THE PARAMETERS
# --------------------------------------------------------

lam = 1.064          # wavelength (um)
thickness = 0.02     # slab thickness (um)
n_orders = 20        # number of Fourier orders (controls reconstruction detail)
grid_per_cell = 140  # sampling resolution per axis, per unit cell

lattice = float(sys.argv[1]) if len(sys.argv) > 1 else 1.21   # lattice (um)
radius = float(sys.argv[2]) if len(sys.argv) > 2 else 0.2     # hole radius (um)
cells = int(sys.argv[3]) if len(sys.argv) > 3 else 3          # cells per axis
outfile = sys.argv[4] if len(sys.argv) > 4 else \
    f"epsilon_l{lattice:g}_r{radius:g}_{cells}x{cells}.png"


def sin_epsilon():
    """SiN permittivity at lam, from the Sellmeier dispersion in simbroad."""
    lam2 = lam ** 2
    eps = (1 + 3.0249 * lam2 / (lam2 - 0.1353406**2)
           + 40314 * lam2 / (lam2 - 1239.842**2))
    return complex(eps)


# --------------------------------------------------------
#   2/ BUILD THE STRUCTURE AND SAMPLE THE DIELECTRIC
# --------------------------------------------------------

if lattice <= 2 * radius:
    sys.exit(f"holes overlap: lattice ({lattice}) <= 2*radius ({2 * radius})")

a_nm = lattice * 1e3

# hexagonal lattice, same vectors as simbroad.py
S = S4.New(Lattice=((a_nm, 0.0),
                    (np.cos(np.pi / 3) * a_nm, np.sin(np.pi / 3) * a_nm)),
           NumBasis=n_orders)

S.SetOptions(LatticeTruncation="Circular", DiscretizedEpsilon=False,
             DiscretizationResolution=8, PolarizationDecomposition=False,
             PolarizationBasis="Default", LanczosSmoothing=False,
             SubpixelSmoothing=False, ConserveMemory=False,
             WeismannFormulation=False, Verbosity=0)

S.SetMaterial(Name="Air", Epsilon=1.0 + 0j)
S.SetMaterial(Name="SiN", Epsilon=sin_epsilon())

S.AddLayer(Name="incidence", Thickness=0.0, Material="Air")
S.AddLayer(Name="slab", Thickness=thickness * 1e3, Material="SiN")
S.AddLayer(Name="transmission", Thickness=0.0, Material="Air")

# circular vacuum hole in the slab
S.SetRegionCircle(Layer="slab", Material="Air", Center=(0.0, 0.0),
                  Radius=radius * 1e3)

S.SetFrequency(1.0 / (lam * 1e3))

# sample the reconstructed dielectric at mid-slab depth, over a window spanning
# `cells` unit cells along each axis (the Fourier reconstruction is periodic)
grid = grid_per_cell * cells
half = cells * a_nm / 2
xs = np.linspace(-half, half, grid)
ys = np.linspace(-half, half, grid)
z_mid = thickness * 1e3 / 2
eps = np.array([[S.GetEpsilon(x, y, z_mid) for x in xs] for y in ys])

# --------------------------------------------------------
#   3/ PLOT AND SAVE
# --------------------------------------------------------

fig, (ax_re, ax_im) = plt.subplots(1, 2, figsize=(11, 4.5))

m_re = ax_re.pcolormesh(xs, ys, eps.real, shading="auto", cmap="viridis")
ax_re.set_title("Re(epsilon)")
ax_re.set_xlabel("x (nm)")
ax_re.set_ylabel("y (nm)")
ax_re.set_aspect("equal")
fig.colorbar(m_re, ax=ax_re)

m_im = ax_im.pcolormesh(xs, ys, eps.imag, shading="auto", cmap="viridis")
ax_im.set_title("Im(epsilon)")
ax_im.set_xlabel("x (nm)")
ax_im.set_ylabel("y (nm)")
ax_im.set_aspect("equal")
fig.colorbar(m_im, ax=ax_im)

fig.suptitle(f"Dielectric constant - lattice {lattice:g} um, "
             f"radius {radius:g} um, {cells}x{cells} cells")

fig.savefig(outfile, dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"lattice = {lattice} um, radius = {radius} um, "
      f"eps(SiN) = {sin_epsilon().real:.4f}")
print(f"Re(eps) range: {eps.real.min():.4f} - {eps.real.max():.4f}")
print(f"saved {outfile}")

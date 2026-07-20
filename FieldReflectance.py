# Field (amplitude) reflectivity for a single lattice/radius point.
#
# The existing *Reflectance.py / *Relative*.py scripts only report the POWER
# reflectance R = Re(-back/forw) from S4's GetPowerFlux. This script instead
# reports the complex specular FIELD reflection coefficient r, i.e. its
# magnitude |r| and phase arg(r), obtained from S4's GetAmplitudes (which returns
# the complex mode amplitudes; the first entry of each polarisation half is the
# specular order). r = back_specular / forw_specular.
#
# Geometry is selectable (hex or square) and both incidence angles (polar and
# azimuthal) are used, matching the conventions in the existing scripts:
#   - hex   : lattice ((a,0),(cos60 a, sin60 a)); SiN via Sellmeier dispersion
#   - square: lattice ((a,0),(0,a));              SiN constant index n = 1.9861
#
# As a built-in correctness check, |r_co|^2 + |r_cross|^2 (both specular
# polarisation channels share the same order/medium/kz) is compared against the
# power reflectance from GetPowerFlux; they should agree to numerical precision.
#
# Usage:
#   python FieldReflectance.py

import numpy as np

import S4  # native S4 Python extension


# --------------------------------------------------------
#   1/ SET THE PARAMETERS
# --------------------------------------------------------

lam = 1.064          # wavelength (um)

shape = "hex"        # "hex" or "square" - selectable geometry

# incidence angles (degrees); polar = rot about y, azimuthal = rot about x
angle_polar = 0.0
angle_azi = 0.0

# the single point to evaluate (um)
lattice = 1.1125#1.1125        # lattice constant a
radius = 0.37#0.37        # hole radius r

thickness = 0.1     # slab thickness (um)
n_orders = 20        # number of Fourier orders (SetNumG in the Lua script)


class SellmeierSiN:
    """SiN with the same Sellmeier dispersion used in HexReflectance.py."""

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


SiN = SellmeierSiN()


# --------------------------------------------------------
#   2/ COMPLEX FIELD REFLECTION COEFFICIENT
# --------------------------------------------------------

def field_reflectance(angle_polar, angle_azi, lattice, radius):
    """Complex specular field reflection coefficient r for one
    (angle_polar, angle_azi, lattice, radius). Same S4 setup as the existing
    scripts; only the output quantity differs (GetAmplitudes instead of just
    GetPowerFlux). Returns a dict with the co-pol r, |r|, phase (deg), the
    field-derived reflectance R_field, and the power reflectance R_power.
    """
    wl_nm = lam * 1e3                       # wavelength in nm
    a_nm = lattice * 1e3                    # lattice constant in nm

    # lattice vectors and SiN permittivity depend on the selected geometry
    if shape == "hex":
        lattice_vecs = ((a_nm, 0.0),
                        (np.cos(np.pi / 3) * a_nm, np.sin(np.pi / 3) * a_nm))
        eps_SiN = complex(SiN.n(lam * 1e-6)) ** 2
    elif shape == "square":
        lattice_vecs = ((a_nm, 0.0), (0.0, a_nm))
        eps_SiN = 1.9861**2 + 0j
    else:
        raise ValueError(f"shape must be 'hex' or 'square', got {shape!r}")

    S = S4.New(Lattice=lattice_vecs, NumBasis=n_orders)

    # same solver options RayFlare applies to its S4 simulations
    S.SetOptions(LatticeTruncation="Circular", DiscretizedEpsilon=False,
                 DiscretizationResolution=8, PolarizationDecomposition=False,
                 PolarizationBasis="Default", LanczosSmoothing=False,
                 SubpixelSmoothing=False, ConserveMemory=False,
                 WeismannFormulation=False, Verbosity=0)

    S.SetMaterial(Name="Air", Epsilon=1.0 + 0j)
    S.SetMaterial(Name="SiN", Epsilon=eps_SiN)

    # incidence / patterned slab / transmission; semi-infinite media have zero
    # thickness in S4, slab thickness in nm
    S.AddLayer(Name="incidence", Thickness=0.0, Material="Air")
    S.AddLayer(Name="slab", Thickness=thickness * 1e3, Material="SiN")
    S.AddLayer(Name="transmission", Thickness=0.0, Material="Air")

    # circular vacuum hole in the slab
    S.SetRegionCircle(Layer="slab", Material="Air", Center=(0.0, 0.0),
                      Radius=radius * 1e3)

    # p-polarised plane wave
    S.SetExcitationPlanewave(IncidenceAngles=(angle_polar, angle_azi),
                             sAmplitude=0.0, pAmplitude=1.0, Order=0)
    S.SetFrequency(1.0 / wl_nm)

    # complex mode amplitudes in the incidence layer: forward = incident,
    # backward = reflected. Each diffraction order has two polarisation/field
    # components: order g sits at basis index p and at p + nb (nb = number of
    # orders). The specular (0,0) order is located from the basis set; under
    # azimuthal rotation the incident p-wave splits across BOTH of its components,
    # so the incident polarisation is the VECTOR (forw[s0], forw[s1]).
    forw, back = S.GetAmplitudes(Layer="incidence")
    forw = np.asarray(forw, dtype=complex)
    back = np.asarray(back, dtype=complex)

    basis = S.GetBasisSet()                  # integer reciprocal coords per order
    nb = len(basis)
    p0 = next(i for i, g in enumerate(basis) if tuple(g) == (0, 0))
    s0, s1 = p0, p0 + nb                      # the two specular components

    # incident specular polarisation vector and its power
    p_in = abs(forw[s0]) ** 2 + abs(forw[s1]) ** 2

    # co-polarised field reflection = reflected specular vector projected onto the
    # incident polarisation direction (reduces to back/forw when azimuth = 0)
    r_co = (np.conj(forw[s0]) * back[s0]
            + np.conj(forw[s1]) * back[s1]) / p_in

    # total specular reflectance from the amplitudes (co- + cross-pol of the
    # specular order); cross-pol power is the remainder
    R_field = (abs(back[s0]) ** 2 + abs(back[s1]) ** 2) / p_in
    R_cross = max(R_field - abs(r_co) ** 2, 0.0)

    # independent power checks: specular order alone (valid at all angles) and the
    # total over every propagating order (exceeds specular when diffraction occurs)
    P = S.GetPowerFluxByOrder(Layer="incidence")
    k = int(np.argmax([abs(p[0].real) for p in P]))   # incident (specular) order
    R_spec_power = float((-P[k][1].real) / P[k][0].real)
    n_prop = sum(1 for p in P if abs(p[1].real) > 1e-12)

    f, b = S.GetPowerFlux(Layer="incidence")
    R_total = float((-b / f).real)

    return {
        "r": r_co,
        "mag": abs(r_co),
        "phase_deg": float(np.degrees(np.angle(r_co))),
        "R_cross": R_cross,
        "R_field": R_field,
        "R_spec_power": R_spec_power,
        "R_total": R_total,
        "n_prop": n_prop,
    }


# --------------------------------------------------------
#   3/ RUN AND REPORT
# --------------------------------------------------------

if __name__ == "__main__":

    res = field_reflectance(angle_polar, angle_azi, lattice, radius)

    print(f"geometry            : {shape}")
    print(f"incidence angles    : polar {angle_polar:g} deg, azi {angle_azi:g} deg")
    print(f"lattice a           : {lattice:g} um")
    print(f"hole radius r       : {radius:g} um  (r/a = {radius / lattice:.4f})")
    print()
    print(f"field reflection r  : {res['r'].real:+.6f} {res['r'].imag:+.6f}j")
    print(f"  magnitude |r|     : {res['mag']:.6f}")
    print(f"  phase arg(r)      : {res['phase_deg']:+.4f} deg")
    print(f"  cross-pol R       : {res['R_cross']:.6e}")
    print()
    # cross-check: the specular reflectance from the amplitudes should match the
    # specular order's power at all angles; the total reflectance additionally
    # counts any propagating diffraction orders.
    diff = res["R_field"] - res["R_spec_power"]
    print(f"check |r|^2 (spec)   : {res['R_field']:.8f}")
    print(f"      power (spec)   : {res['R_spec_power']:.8f}  (diff {diff:+.1e})")
    print(f"      power (total)  : {res['R_total']:.8f}  ({res['n_prop']} propagating order(s))")

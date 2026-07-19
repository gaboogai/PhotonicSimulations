# How S4 works: the physics of RCWA / the Fourier Modal Method

This note explains the electromagnetic method behind **S4** (Stanford Stratified
Structure Solver) — the **Rigorous Coupled-Wave Analysis (RCWA)**, also called the
**Fourier Modal Method (FMM)**. It is written to accompany the scripts in this
repository (which drive S4's native Python API), and it ties the theory back to the
specific functions those scripts call: `New`, `SetExcitationPlanewave`,
`GetAmplitudes`, `GetBasisSet`, `GetPowerFlux`, and `GetFields`.

---

## 1. The problem S4 solves

S4 computes how light scatters from a **stack of layers that are periodic in the
two in-plane directions (x, y) and uniform along the propagation direction (z)
within each layer**. Each layer can be patterned (e.g. a lattice of holes), but its
cross-section does not change with z; a 3-D structure is built as a *stack* of such
z-invariant layers.

For our scripts the stack is:

```
  incidence (semi-infinite air)
  slab       (SiN, thickness t, patterned with a circular air hole per unit cell)
  transmission (semi-infinite air)
```

Given an incident plane wave, S4 returns the reflected and transmitted fields,
their power, and — crucially for `FieldReflectance.py` — the **complex** amplitude
(magnitude **and** phase) of each diffraction order.

This is a *frequency-domain*, *rigorous* solver: it solves Maxwell's equations
exactly (up to a controllable Fourier truncation), with no ray or scalar
approximation. It is extremely efficient for **layered periodic** media, and is not
a general-geometry solver like FDTD or FEM.

---

## 2. Physical picture: periodicity forces discrete diffraction orders

Because the permittivity is periodic in-plane with lattice vectors
$\mathbf{a}_1, \mathbf{a}_2$, the **Bloch–Floquet theorem** applies: the field is a
periodic function times a plane-wave phase. A periodic medium can only exchange
in-plane momentum with the wave in discrete units of the **reciprocal lattice
vectors** $\mathbf{G}_1, \mathbf{G}_2$ (defined by
$\mathbf{a}_i\!\cdot\!\mathbf{G}_j = 2\pi\,\delta_{ij}$).

So an incident wave with in-plane wavevector $\mathbf{k}_\parallel$ can only scatter
into a **discrete set of diffraction orders**:

$$
\mathbf{k}_\parallel^{(m,n)} = \mathbf{k}_\parallel + m\,\mathbf{G}_1 + n\,\mathbf{G}_2,
\qquad (m,n)\in\mathbb{Z}^2 .
$$

The **specular order** is $(m,n)=(0,0)$ — the ordinary reflection/transmission that
keeps the incident in-plane direction. All other orders are diffracted beams.

Each order has an out-of-plane wavevector fixed by the dispersion relation in a
given medium of index $n_\text{med}$:

$$
k_z^{(m,n)} = \sqrt{\,(n_\text{med}\,k_0)^2 - \big|\mathbf{k}_\parallel^{(m,n)}\big|^2\,},
\qquad k_0 = \frac{2\pi}{\lambda}.
$$

- If $k_z$ is **real**, the order **propagates** and carries power to infinity.
- If $k_z$ is **imaginary**, the order is **evanescent** — it decays away from the
  structure and carries no net power (but still stores energy and matters near
  interfaces).

**This is the single most important idea for interpreting our results.** At normal
incidence on a sub-wavelength lattice, only the specular order propagates, so *all*
reflected power is specular. As you tilt the beam or make the period approach the
wavelength, higher orders cross from evanescent to propagating (the **Rayleigh–Wood
anomaly**) and start carrying power. That is exactly why, in `FieldReflectance.py`,
the specular $|r|^2$ equals the total power reflectance only when `n_prop == 1`, and
is smaller than it once diffraction orders open up.

---

## 3. Maxwell's equations in a Fourier basis

Within one z-invariant layer, S4 expands the periodic permittivity and all field
components as **Fourier series over the reciprocal lattice**:

$$
\varepsilon(x,y) = \sum_{\mathbf{G}} \hat\varepsilon_{\mathbf{G}}\,
e^{i\mathbf{G}\cdot \mathbf{r}_\parallel},
\qquad
\mathbf{E}(\mathbf{r}) = \sum_{\mathbf{G}}
\mathbf{E}_{\mathbf{G}}(z)\, e^{i(\mathbf{k}_\parallel+\mathbf{G})\cdot \mathbf{r}_\parallel}.
$$

The infinite sum is **truncated** to a finite set of $N$ retained $\mathbf{G}$
vectors — this is the `NumBasis` (a.k.a. `NumG`) parameter in `S4.New(...)`. The set
of retained $\mathbf{G}$'s is what `GetBasisSet()` reports.

Substituting these expansions into Maxwell's two curl equations turns them into a
set of **coupled ordinary differential equations in z** for the Fourier
coefficients. The coupling comes from the products like $\varepsilon\,\mathbf{E}$,
which become **convolutions** (Toeplitz matrices $\llbracket\varepsilon\rrbracket$)
in Fourier space. Schematically:

$$
\frac{d}{dz}
\begin{pmatrix}\mathbf{E}_\parallel\\ \mathbf{H}_\parallel\end{pmatrix}
= i\,\mathbf{M}\,
\begin{pmatrix}\mathbf{E}_\parallel\\ \mathbf{H}_\parallel\end{pmatrix},
$$

where $\mathbf{M}$ is a dense matrix built from the truncated $\mathbf{G}$'s, the
convolution matrices of $\varepsilon$ (and $1/\varepsilon$), and $k_0$.

---

## 4. Eigenmodes: the "coupled waves"

Because a layer is uniform in z, the z-dependence of each eigen-solution of the
system above is a pure exponential $e^{i q z}$. Solving the eigenvalue problem for
$\mathbf{M}$ gives, per layer:

- a set of **eigenvalues** $q$ → the effective out-of-plane wavevectors of the
  layer's Bloch modes (propagating or evanescent), and
- **eigenvectors** → the transverse field profiles of those modes.

These eigenmodes are the "coupled waves" of RCWA. Any field in the layer is a linear
combination of **forward** ($+z$) and **backward** ($-z$) modes with complex
**amplitudes**. For a *uniform* layer (like the incidence/transmission air regions),
the eigenmodes are simply the individual diffraction orders — which is why
`GetAmplitudes` in those layers is directly interpretable as per-order amplitudes.

---

## 5. Fourier factorization and convergence (why the options exist)

Naively Fourier-expanding products of discontinuous functions (fields and
$\varepsilon$ both jump at material boundaries) converges *slowly and non-uniformly*,
especially for the field component **normal** to an interface and for TM/p
polarization. The fix is the **correct Fourier factorization rules** (Li's rules):
use the direct rule for products continuous across the interface and the *inverse
rule* (with $\llbracket 1/\varepsilon\rrbracket^{-1}$) for the discontinuous normal
component.

For 2-D crossed gratings this requires resolving field components **relative to the
local boundary normal**, i.e. building a smooth vector field tangent/normal to the
pattern edges. In S4 this is the **`PolarizationDecomposition`** /
`PolarizationBasis` option (the "normal vector method"), and
`SubpixelSmoothing` / `DiscretizedEpsilon` / `LanczosSmoothing` are related knobs
that trade accuracy against speed/robustness. Enabling proper factorization
dramatically improves convergence with respect to `NumBasis`. (Our scripts follow
the RayFlare defaults and rely on a moderately large `NumBasis` instead.)

See S4's *FMM formulations* documentation page for the exact options and the theory
they implement.

---

## 6. Connecting layers: the scattering matrix

Each layer gives eigenmodes and amplitudes; the layers are joined by matching the
tangential $\mathbf{E}$ and $\mathbf{H}$ at every interface. Doing this by
multiplying transfer (T) matrices is numerically **unstable** because evanescent
modes produce exponentially growing terms. S4 instead uses the **scattering-matrix
(S-matrix) method**, which relates *outgoing* amplitudes to *incoming* amplitudes
and keeps only decaying exponentials — stable for arbitrarily thick or lossy stacks.
The excitation (`SetExcitationPlanewave`) sets the incoming amplitude of the
specular order in the incidence medium; solving the S-matrix yields all reflected
and transmitted mode amplitudes.

---

## 7. What the output functions actually give

| S4 function | Returns | In our scripts |
|---|---|---|
| `GetBasisSet()` | The retained reciprocal-lattice orders as integer pairs $(m,n)$. Length = number of orders actually kept, `nb`. | Locate the specular order `(0,0)` and get the block stride `nb`. |
| `GetAmplitudes(Layer)` | `(forw, back)`: complex amplitudes of the forward/backward eigenmodes in that layer. In a uniform layer these are the per-order, per-polarization amplitudes. | Extract the complex specular reflection coefficient. |
| `GetPowerFlux(Layer)` | Integrated z-component of the (complex) Poynting vector — **total** power summed over all orders. | Cross-check; this is what the other `*Reflectance` scripts report as `R`. |
| `GetPowerFluxByOrder(Layer)` | Poynting flux split **per diffraction order**. | Verify the specular-order power matches `\|r\|^2`. |
| `GetFields(x,y,z)` | Complex $\mathbf{E}, \mathbf{H}$ at a point (real-space Fourier reconstruction). | (Not used, but an alternative route to phase.) |

### The amplitude array layout (the subtle part)

`GetAmplitudes` returns a flat array of length `2 × NumBasis`. It is **two stacked
blocks** — one per transverse field/polarization component (S4 works in an
$E_x$-/$E_y$-like Fourier basis, *not* physical s/p). Order sitting at basis
position `p` therefore appears **twice**: at index `p` and at index `p + nb`, where
`nb = len(GetBasisSet())`.

Three consequences, all verified empirically in `FieldReflectance.py`:

1. **The specular order is not index 0.** Its position `p0` is wherever `(0,0)`
   falls in `GetBasisSet()` (e.g. index 19 for our hex lattice, 13 for square).
2. **The block stride is `nb`, not `len(forw)//2`.** With circular truncation the
   kept-order count (`nb = 19` for hex) is less than `NumBasis = 20`, and the array
   is padded to `2·NumBasis`, so the two specular components are at `p0` and
   `p0 + nb`.
3. **A tilted (azimuthal) p-wave splits across both components.** The incident
   polarization is a *vector* $(f_{s_0}, f_{s_1})$ over the two specular indices
   (we measured $\cos\phi, -\sin\phi$). So the co-polarized field reflection
   coefficient is the reflected specular vector **projected onto the incident
   direction**:

$$
r = \frac{f_{s_0}^{*}\,b_{s_0} + f_{s_1}^{*}\,b_{s_1}}
         {|f_{s_0}|^2 + |f_{s_1}|^2},
\qquad
R_\text{spec} = \frac{|b_{s_0}|^2 + |b_{s_1}|^2}{|f_{s_0}|^2 + |f_{s_1}|^2}.
$$

$|r|$ is the field-reflectivity magnitude and $\arg(r)$ its phase (referenced to the
plane $z=0$, the top of the incidence layer). The phase is well defined because the
real, positive per-mode power-normalization factor is identical for the forward and
backward of the same order/medium and cancels in the ratio. `R_spec` computed this
way matches `GetPowerFluxByOrder` for the specular order to ~1e-7 across all angles
and both lattices.

**Amplitude vs. power.** `GetPowerFlux` gives a real total reflectance and discards
phase and the per-order split; `GetAmplitudes` retains the full complex, per-order
information. "Field reflectivity, magnitude **and** angle" needs the latter.

---

## 8. Truncation, the lattice, and accuracy

- `NumBasis` sets how many $\mathbf{G}$'s are kept. Cost scales roughly as
  $\mathcal{O}(N^3)$ in time and $\mathcal{O}(N^2)$ in memory. Convergence should be
  checked by increasing `NumBasis` until outputs stop changing.
- `LatticeTruncation="Circular"` keeps the shortest $\mathbf{G}$'s (by length),
  which is why the kept count (`nb`) is not a perfect square and can be slightly
  below `NumBasis`.
- The lattice vectors are set in `S4.New(Lattice=...)`; our hex script uses
  $((a,0),(a\cos 60^\circ, a\sin 60^\circ))$ and the square script $((a,0),(0,a))$.
- Sharp features and metals converge more slowly and benefit from the factorization
  options in §5.

---

## 9. Assumptions and limitations

- **Linear, non-magnetic** ($\mu=1$), local, frequency-domain response; materials
  enter only through $\varepsilon(\omega)$.
- Structure must be a **stack of in-plane-periodic, z-uniform layers**. Aperiodic or
  fully arbitrary 3-D geometry is not the tool's domain (use FDTD/FEM there).
- Results are **exact only in the limit of infinite `NumBasis`**; always do a
  convergence study.
- Truncation-induced artifacts (Gibbs ringing) are worst at high-contrast, metallic,
  or sharp-cornered interfaces.

---

## 10. Relation to the scripts in this repo

- `HexReflectance.py`, `SquareReflectance.py`, `HexRelative.py`,
  `SquareRelativeReflective.py` → **power** reflectance via `GetPowerFlux`
  ($R=\mathrm{Re}(-\text{back}/\text{forw})$).
- `FieldReflectance.py` → **complex field** reflectivity (magnitude + phase) via
  `GetAmplitudes` + `GetBasisSet`, cross-checked against `GetPowerFluxByOrder`.
- `FieldReflectanceSweep.py` → the same field quantity swept over a
  (lattice, radius-ratio) grid, writing magnitude and phase CSVs and heat-mapping
  them with `PlotRelative.py`.

---

## 11. Sources for further study

**The S4 method and software**
- V. Liu and S. Fan, "S4: A free electromagnetic solver for layered periodic
  structures," *Computer Physics Communications* **183**, 2233–2244 (2012).
  (The S4 paper — describes exactly this implementation.)
- S4 documentation and *FMM formulations* page:
  <https://web.stanford.edu/group/fan/S4/> and the local `python_api.rst`.

**Foundational RCWA / FMM papers**
- M. G. Moharam and T. K. Gaylord, "Rigorous coupled-wave analysis of planar-grating
  diffraction," *JOSA* **71**, 811 (1981).
- M. G. Moharam, E. B. Grann, D. A. Pommet, T. K. Gaylord, "Formulation for stable
  and efficient implementation of the rigorous coupled-wave analysis of binary
  gratings," *JOSA A* **12**, 1068 (1995).
- M. G. Moharam, D. A. Pommet, E. B. Grann, T. K. Gaylord, "Stable implementation of
  the rigorous coupled-wave analysis … enhanced transmittance matrix approach,"
  *JOSA A* **12**, 1077 (1995).

**Fourier factorization and convergence (the "Li rules")**
- L. Li, "Use of Fourier series in the analysis of discontinuous periodic
  structures," *JOSA A* **13**, 1870 (1996).
- L. Li, "New formulation of the Fourier modal method for crossed surface-relief
  gratings," *JOSA A* **14**, 2758 (1997).
- P. Lalanne and G. M. Morris, "Highly improved convergence of the coupled-wave
  method for TM polarization," *JOSA A* **13**, 779 (1996).
- G. Granet and B. Guizal, "Efficient implementation of the coupled-wave method for
  metallic lamellar gratings in TM polarization," *JOSA A* **13**, 1019 (1996).

**Scattering-matrix formulation (numerical stability)**
- D. M. Whittaker and I. S. Culshaw, "Scattering-matrix treatment of patterned
  multilayer photonic structures," *Physical Review B* **60**, 2610 (1999).
- L. Li, "Formulation and comparison of two recursive matrix algorithms for modeling
  layered diffraction gratings," *JOSA A* **13**, 1024 (1996).

**Books and tutorials**
- E. Popov (ed.), *Gratings: Theory and Numeric Applications* (2012; freely
  available online) — comprehensive, method-by-method.
- R. C. Rumpf, computational electromagnetics lectures and RCWA tutorials,
  <https://empossible.net/academics/> — practical, step-by-step derivations.
- J. D. Joannopoulos et al., *Photonic Crystals: Molding the Flow of Light*, 2nd ed.
  (2008) — background on Bloch modes, diffraction, and photonic periodic media.

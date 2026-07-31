# PhotonicSimulations
Simulate photonic crystals using S4 python extension

## Scripts

### simbroad.py
Sweeps a hole-patterned SiN slab over the full radius and lattice ranges (for each incidence angle) and writes the reflectance grid to `cleanbroad<angle>.csv`.

```
python simbroad.py
```

### simbroad_single.py
Runs the same calculation as `simbroad.py` for a single radius, lattice and angle (set at the top of the file) and writes the one-point reflectance to `cleanbroad<angle>.csv`.

```
python simbroad_single.py
```

### plotReflectance.py
Plots the reflectance from a single `cleanbroad*.csv` output as a heatmap over the lattice/radius grid and saves it as a PNG.

```
python plotReflectance.py [cleanbroad0.csv]
```

### plotReflectanceDiff.py
Plots the difference in reflectance between two `cleanbroad*.csv` outputs as a diverging heatmap and saves it as a PNG.

```
python plotReflectanceDiff.py file_a.csv file_b.csv [output.png]
```

### layerplot.py
Plots the dielectric constant (real and imaginary components) of the patterned SiN slab for a specified lattice and hole radius, optionally tiled over several unit cells, and saves it as a PNG.

```
python layerplot.py [lattice_um] [radius_um] [cells] [output.png]
```

### HexRelative.py
Variant of `simbroad.py` where the hole radius is defined relative to the lattice (radius = ratio × a). Sweeps the lattice and radius-ratio ranges (for each angle) and writes a CSV to `RawData/` whose columns are radius ratios (r/a) and whose rows are lattice constants.

```
python HexRelative.py
```

### PlotRelative.py
Plots the reflectance from a single `HexRelative.py` output as a heatmap (radius ratio vs lattice constant) and saves it as a PNG. With no argument it plots the most recent `RawData/HexRelative*.csv`.

```
python PlotRelative.py [path/to/HexRelative...csv]
```

### FeasibilityPlot.py
Plots the reflectance difference between two `HexRelative.py` outputs (typically two different angles), masked to cells where the first file's reflectance exceeds 0.15, and saves it as a PNG. Highlights designs that are both highly reflective and angle-robust.

```
python FeasibilityPlot.py file_a.csv file_b.csv [output.png]
```

### FieldReflectance.py
Computes the complex specular field reflection coefficient r — its magnitude `|r|` and phase `arg(r)` — for a single lattice, radius and incidence angle (polar and azimuthal), for a selectable hex or square lattice. Uses S4's mode amplitudes (`GetAmplitudes`) and cross-checks `|r|²` against the power reflectance.

```
python FieldReflectance.py
```

### FieldReflectanceSweep.py
Sweep version of `FieldReflectance.py`: computes the complex field reflectivity over the lattice and radius-ratio grid (radius = ratio × a) for each incidence angle, writing two CSVs per angle (magnitude and phase) to `FieldSweep/`, and auto-plots each with `PlotRelative.py`.

```
python FieldReflectanceSweep.py
```

### photonicTest.py
Fixes the lattice constant and hole radius and sweeps the two incidence angles (polar × azimuthal) on a 2D grid, writing the field reflectivity magnitude and phase to two CSVs in `AngleSweep/`. Auto-plots the pair with `plotPhotonicTest.py`.

```
python photonicTest.py
```

### plotPhotonicTest.py
Plots a `photonicTest.py` magnitude/phase pair: heatmaps of `|r|` and `arg(r)` over the polar/azimuthal grid, plus 1D sections of `|r|` vs polar angle at the first and last azimuthal angles. Saves a PNG.

```
python plotPhotonicTest.py Magnitude.csv Phase.csv [output.png]
```

### kSpaceConvert.py
Converts a `photonicTest.py` magnitude/phase pair from incidence angles (polar, azimuthal) to in-plane wavevector components (kx, ky), writing a long-format `...KSpace.csv`. Passing `plot=true` also runs `plotKSpace.py`; either way it auto-runs `fitPhotonicTest.py`.

```
python kSpaceConvert.py Magnitude.csv Phase.csv [output.csv] [plot=true]
```

### plotKSpace.py
Plots a `...KSpace.csv` as two 2D maps over the wavevector plane: field reflectivity magnitude `|r|` vs (kx, ky) and phase `arg(r)` vs (kx, ky). Saves a PNG.

```
python plotKSpace.py path/to/...KSpace.csv [output.png]
```

### extendKspace.py
Extends a `...KSpace.csv` around a full 360° of azimuth by tiling the swept wedge with the lattice's rotational symmetry, writing a `...KSpace360.csv`. The period is detected from the azimuthal span (60° for hex, 90° for square) and checked by confirming the first and last azimuth slices match; it errors if the span is neither, or if that seam check fails. Passing `plot=true` also runs `plotKSpace.py`.

```
python extendKspace.py path/to/...KSpace.csv [output.csv] [plot=true] [period=60]
```

### fitPhotonicTest.py
Fits analytic functions of the field reflectivity, `r_magnitude(kx, ky)` and `r_phase(kx, ky)`, to a `...KSpace.csv` — a coupled-mode/Fano resonance model and symmetry-adapted polynomial surrogates. Reports goodness-of-fit, and saves the fit parameters (JSON) and a data-vs-fit diagnostic PNG.

```
python fitPhotonicTest.py path/to/...KSpace.csv
```
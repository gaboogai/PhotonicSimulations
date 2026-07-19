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
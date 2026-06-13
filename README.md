# PhotonicSimulations
Simulate photonic crystals using S4 python extension

## Scripts

### simbroad.py
Sweeps a hole-patterned SiN slab over the full radius and lattice ranges (for each incidence angle) and writes the reflectance grid to `cleanbroad<angle>.csv`.

```
python simbroad.py
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
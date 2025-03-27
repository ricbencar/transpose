# OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION
*(Enhanced Approach Incorporating Wave Obliquity and Refraction)*

## Overview

This program processes wave data from an input CSV file, computes nearshore wave parameters at a specified depth using linear wave theory enhanced with refraction and shoaling effects, and generates:

- **`output.csv`** – Contains the computed results for each time step.
- **`report.txt`** – Provides descriptive statistics of both input and computed variables, including annual maxima.

![transpose](https://github.com/user-attachments/assets/4e49cd6a-75bf-4e3d-9e5d-d905f913f08e)
---
## Computed Parameters

| Parameter         | Description |
|------------------|-------------|
| **L0** | Deep-water wavelength: `L0 = g * T² / (2π)` |
| **L** | Local wavelength, solved via Newton-Raphson from `L = L0 * tanh((2π * depth_d) / L)` |
| **kh** | Product of the wave number (`k = 2π / L`) and local depth (`depth_d`) |
| **alpha_offshore** | Signed offshore wave obliquity (crest-to-coast difference in degrees), considering approach relative to the coastline. |
| **alpha_local** | Local wave angle (degrees) after refraction, derived using Snell's law. |
| **mwd_local** | Local mean wave direction (degrees), adjusted from offshore `mwd` based on the change in wave angle (`alpha_offshore - alpha_local`). |
| **Ks** | Shoaling coefficient, `Ks = sqrt(Cg0 / Cg)`. |
| **Kr** | Refraction coefficient, `Kr = sqrt(cos(alpha_offshore) / cos(alpha_local))`. |
| **Hb** | Breaking wave height (Miche, 1944): `Hb = 0.142 * L * tanh(kh)`. |
| **swh_local** | Local significant wave height, calculated as the minimum of the transformed height (`swh_offshore * Ks * Kr`) and the breaking height (`Hb`). |

**Note:** Waves arriving from directions between `coast_dir` and `coast_dir + 180°` (clockwise, i.e., from the land side relative to the specified coastline orientation), or waves with non-positive offshore height (`swh <= 0`), result in `swh_local` and other derived local parameters (L, kh, alpha_local, mwd_local, Ks, Kr, Hb) being set to **zero** for that time step.

---
## Results (Report File)

The `report.txt` file includes:

- The exact command line used to invoke the program.
- **Descriptive Statistics:** For each variable (input and computed), the report provides:
    - Count
    - Mean
    - Standard Deviation
    - Minimum & Maximum
    - Median (50th percentile)
    - Percentiles: 1%, 10%, 25%, 75%, 90%, 99%
- **Annual Maxima:** A table showing the maximum `swh_offshore` and `swh_local` for each year present in the data, plus the overall maximum across all years.

**Important Note on Statistics:**
- For the variables **`alpha_local`**, **`swh_local`**, **`mwd_local`**, **`Ks`**, **`Kr`**, and **`Hb`**, the descriptive statistics (count, mean, stddev, min, max, percentiles) are calculated **excluding** any time steps where the computed `swh_local` is zero. This effectively removes waves originating from the land side or those with zero initial offshore height from these specific statistical summaries.
- Statistics for all other variables (`swh_offshore`, `mwd_offshore`, `pp1d`, `L0`, `L`, `kh`, `alpha_offshore`) include all valid input time steps.

**Directional Statistics (`mwd_offshore`, `mwd_local`):**
A hybrid approach is used:
- **Circular Mean** and **Circular Standard Deviation** are computed using the unit-vector method.
- Minimum, Maximum, Median, and Quantiles are calculated using ordinary linear statistics on the angles wrapped to the range [0, 360).
- For `mwd_local`, these statistics also exclude time steps where `swh_local` is zero, consistent with the note above.

---

## Usage

```sh
./transpose input_csv coast_dir depth_d
```

### Arguments:
- **`input_csv`** : CSV input file (with columns: `datetime`, `swh`, `mwd`, `pp1d`)
- **`coast_dir`** : Coastline orientation in degrees (clockwise from North)
- **`depth_d`** : Local depth (meters)

---

## CSV Input Format

The input CSV file should be comma-separated with at least the following columns:

```csv
datetime, swh, mwd, pp1d, [additional columns ignored]
```

---

## CSV Output Format

The generated `output.csv` will contain the following comma-separated columns:

```csv
datetime,swh_offshore,mwd_offshore,pp1d,L0,L,kh,alpha_offshore,alpha_local,swh_local,mwd_local,Ks,Kr,Hb
```

---

## Compilation

To compile the program, use the following command:

```sh
g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o transpose transpose.cpp
```

This command enables **optimizations** and includes several **compiler warnings** to ensure code quality.

---

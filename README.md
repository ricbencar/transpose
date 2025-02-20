# OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION
*(Simplified Approach Using Linear Wave Theory)*

## Overview

This program processes wave data from an input CSV file, computes nearshore wave parameters at a specified depth, and generates:

- **`output.csv`** – Contains the computed results.
- **`report.txt`** – Provides descriptive statistics of both input and computed variables.

For directional wave data (i.e. mwd_offshore and mwd_local), the program employs a hybrid approach: the circular mean and circular standard deviation are computed using the unit‑vector method, while the minimum, maximum, median, and quantiles are calculated using ordinary linear statistics on the wrapped angles (in the range [0,360)).

---

## USAGE

```sh
./transpose input_csv coast_dir depth_d
```

### Arguments:
- **`input_csv`** : CSV input file (with columns: `datetime`, `swh`, `mwd`, `pp1d`)
- **`coast_dir`** : Coastline orientation in degrees (clockwise from North)
- **`depth_d`** : Local depth (meters)

---

## CSV INPUT FORMAT

The input CSV file should be comma-separated with at least the following columns:

```csv
datetime, swh, mwd, pp1d, [additional columns ignored]
```

---

## OUTPUT CSV FORMAT

The generated `output.csv` will contain the following comma-separated columns:

```csv
datetime,swh_offshore,mwd_offshore,pp1d,L0,L,kh,alpha_offshore,alpha_local,swh_local,mwd_local,Ks,Kr,Hb
```

---

## Computed Parameters

| Parameter         | Description |
|------------------|-------------|
| **L0** | Deep-water wavelength: `L0 = g * T² / (2π)` |
| **L** | Local wavelength, solved from `L = L0 * tanh((2π * depth_d) / L)` |
| **kh** | Wave number (`k = 2π / L`) times local depth (`h`) |
| **alpha_offshore** | Offshore wave approach angle relative to coastline |
| **alpha_local** | Local wave angle after refraction |
| **mwd_local** | Local mean wave direction, adjusted from offshore `mwd` |
| **Ks** | Shoaling coefficient |
| **Kr** | Refraction coefficient |
| **Hb** | Breaking wave height (Miche, 1944): `Hb = 0.142 * L * tanh((2π * depth_d) / L)` |
| **swh_local** | Local significant wave height (minimum of `swh * Ks * Kr` and `Hb`) |

**Note:** Waves arriving from directions between `coast_dir` and `coast_dir + 180°` (i.e., from the land side) are set to **zero**.

---

## Report File Details

The report.txt file provides:

- **A descriptive statistics report for each output variable with additional percentiles at 1%, 10%, 25%, 50% (median), 75%, 90%, and 99%.
- **A table displaying the annual maxima for swh_offshore and swh_local, with the final row indicating the overall maximum for each variable.
The command line used to run the program at the top of the report.

---

## COMPILATION

To compile the program, use the following command:

```sh
g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o transpose transpose.cpp
```

Compile Options Explained:

    -O3: Enables high-level optimizations for maximum performance.
    -fopenmp: Enables OpenMP support for multi-threading.
    -march=native: Optimizes the code for the architecture of the compiling machine.
    -std=c++17: Uses the C++17 standard.
    -Wall -Wextra -pedantic: Activates a comprehensive set of compiler warnings to ensure code quality.
    -Wconversion: Warns about implicit type conversions.
    -Wsign-conversion: Warns about implicit sign conversions.
    -static -static-libgcc -static-libstdc++: Links libraries statically to enhance portability.

This command enables **optimizations** and includes several **compiler warnings** to ensure code quality.

---
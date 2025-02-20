--------------------------------------------------------------------\
  OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION\
  (Simplified Approach Using Linear Wave Theory and Hybrid Circular Statistics)

  Overview:\
    This program processes wave data from an input CSV file and computes\
    nearshore wave parameters at a specified depth. It generates two output files:\
      - "output.csv" -- Contains the computed nearshore wave parameters.\
      - "report.txt" -- Provides descriptive statistics for both the input and\
                       computed variables.\
    The report includes:\
      * The command line used to invoke the program.\
      * Descriptive statistics for each variable (count, mean, standard deviation,\
        minimum, maximum, median, and percentiles at 1%, 10%, 25%, 50%, 75%, 90%,\
        and 99%).\
      * A table of annual maxima for swh_offshore and swh_local, with a final row\
        indicating the overall maximum values.

    For directional wave data (mwd_offshore and mwd_local), a hybrid approach is used:\
      - The circular mean and circular standard deviation are computed using the\
        unit-vector method.\
      - The minimum, maximum, median, and quantiles are calculated using ordinary\
        linear statistics on the wrapped angles (in [0,360)).

  USAGE:\
    ./transpose input_csv coast_dir depth_d

    Where:\
      input_csv : CSV input file containing at least the following columns:\
                  datetime, swh, mwd, pp1d (additional columns are ignored)\
      coast_dir : Coastline orientation in degrees (clockwise from North)\
      depth_d   : Local depth in meters

  EXPECTED CSV INPUT FORMAT (comma-separated):\
      datetime, swh, mwd, pp1d, [additional columns ignored]

  OUTPUT CSV FORMAT (comma-separated):\
      datetime,swh_offshore,mwd_offshore,pp1d,L0,L,kh,alpha_offshore,\
      alpha_local,swh_local,mwd_local,Ks,Kr,Hb

  Explanation of computed parameters:\
      L0             : Deep-water wavelength, calculated as (g * T²) / (2π)\
      L              : Local wavelength, solved via Newton-Raphson from\
                       L = L0 * tanh((2π * depth_d) / L)\
      kh             : Product of the wave number (k = 2π / L) and local depth (h)\
      alpha_offshore : Offshore wave approach angle relative to the coastline\
      alpha_local    : Local wave angle after refraction\
      mwd_local      : Local mean wave direction, adjusted from the offshore mwd\
      Ks             : Shoaling coefficient\
      Kr             : Refraction coefficient\
      Hb             : Breaking wave height (per Miche, 1944), computed as\
                       Hb = 0.142 * L * tanh((2π * depth_d) / L)\
      swh_local      : Local significant wave height, computed as the minimum of\
                       (swh * Ks * Kr) and Hb

    Note: Waves arriving from directions between coast_dir and coast_dir+180°\
          (i.e., from the land side) are set to zero.

  Report File Details:\
    The report.txt file includes:\
      - The exact command line used to run the program.\
      - Detailed descriptive statistics for each variable, including count,\
        mean, standard deviation, minimum, maximum, median, and percentiles at\
        1%, 10%, 25%, 50%, 75%, 90%, and 99%.\
      - A table displaying the annual maxima for swh_offshore and swh_local,\
        with the final row indicating the overall maximum values.

  Compilation Details:\
    To compile the program, use the following command:

      g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic\
          -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++\
          -o transpose transpose.cpp

    Explanation of compile options:\
      - -O3                   : Enables high-level optimizations for maximum performance.\
      - -fopenmp              : Enables OpenMP support for multi-threading.\
      - -march=native         : Optimizes the code for the architecture of the compiling machine.\
      - -std=c++17            : Uses the C++17 standard.\
      - -Wall -Wextra -pedantic: Activates a broad set of compiler warnings to ensure code quality.\
      - -Wconversion          : Warns about implicit type conversions.\
      - -Wsign-conversion     : Warns about implicit sign conversions.\
      - -static, -static-libgcc, -static-libstdc++: Links libraries statically, enhancing portability.\
  --------------------------------------------------------------------

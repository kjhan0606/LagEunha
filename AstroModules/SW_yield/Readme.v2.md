# Stellar Wind Yield Library (SWL)

A lightweight, thread-safe C library for calculating metal yields from stellar winds in cosmological and galaxy formation simulations (e.g., **LagEunha**).

## Key Features

- **Massive & AGB Winds**: Covers the full mass range from low-mass AGB progenitors to Very Massive Stars (VMS).
- **Refined Physics (v2.1)**:
  - **WN/WC Phase Separation**: Distinguishes between Nitrogen-rich (early) and Carbon-rich (late) Wolf-Rayet winds.
  - **Pop III Support**: Handles extremely low metallicity environments by suppressing line-driven winds.
  - **Super-AGB Transition**: Smooth interpolation in the $8 - 10 M_{\odot}$ transition zone.
  - **Correct Low-Z Scaling**: Primary elements (C, O) are preserved at low Z, while Secondary elements (N) are scaled.
- **Multiple IMFs**: Supports Salpeter, Kroupa, Chabrier, and Baldry & Glazebrook.
- **OpenMP/MPI Ready**: Designed without global mutable state for safe parallel execution.

## Scientific Background

### 1. Primary vs. Secondary Yields
At low metallicities ($Z < 10^{-4}$), the library distinguishes between:
- **Primary Elements (C, O)**: Synthesized by the star's own fusion (Triple-alpha). Yields remain robust even at $Z \approx 0$.
- **Secondary Elements (N)**: Dependent on CNO seeds. Yields scale down linearly with initial $Z$.

### 2. Wolf-Rayet Phases (WN vs. WC)
For massive stars ($M > 25 M_{\odot}$), mass loss is time-separated:
- **WN Phase (Age ~90-95%)**: Envelope stripping releases He and N.
- **WC/WO Phase (Age >95%)**: Core exposure releases C and O.
This separation is critical for accurate cooling calculations in hydrodynamic simulations.

## Quick Start

### Build
```bash
gcc -c stellar_wind_lib.c -O3 -fPIC
ar rcs libstellar_wind.a stellar_wind_lib.o
```

## Usage example
```c
#include "stellar_wind_lib.h"

SWL_ClusterWindYield result;
swl_init();

// Calculate for a 10^6 Msun cluster, solar metallicity, at 4 Myr (WR phase)
swl_calculate_cluster_wind_yield(SWL_IMF_CHABRIER, 0.02, 1e6, 4.0e6, &result);

printf("Nitrogen Yield: %e Msun\n", result.element_yields.yields[SWL_EL_N]);
printf("Carbon Yield:   %e Msun\n", result.element_yields.yields[SWL_EL_C]);
```

# Exact Cooling Integrator

A high-performance C library for exact time integration of radiative cooling and heating in astrophysical gas simulations, based on [Zhu, Smith & Hernquist (2017), MNRAS 470, 1017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.1017Z).

## Overview

This package provides an **exact** (non-approximate) solution to the cooling/heating equation:

$$\frac{du}{dt} = -\Lambda(u)$$

where $u$ is the specific internal energy and $\Lambda$ is the net cooling rate. Unlike simple Euler integration, this method:

- **Preserves accuracy** even for timesteps larger than the cooling time
- **Handles stiff cooling** without numerical instability
- **Correctly captures equilibrium temperatures**

## Features

- #  Exact time integration using the Y-function transformation (Eq. 4-10 of ZSH17)
- #  Pre-computed lookup tables for O(1) fast interpolation
- #  OpenMP parallelization for table construction
- #  Thread-safe integration functions
- #  Coupled Saha equation solver for mean molecular weight
- #  Support for variable metallicity and redshift

## File Structure

```
+-- exact_cooling_integrator.c   # Main integration library
+-- exact_cooling_integrator.h   # Public API header
+-- cooling_rate_mpi.h           # GridParams structure definition
+-- getMU.c                      # Mean molecular weight calculator
+-- getMU.h                      # getMU header
+-- generate_sample_table.c      # Sample cooling table generator
+-- Makefile                     # Build system
+-- README.md                    # This file
```

## Requirements

- GCC or compatible C99 compiler
- OpenMP (optional, for parallel table building)
- Math library (`-lm`)

## Building

```bash
# Build test executable
make all

# Build static library only
make lib

# Clean build artifacts
make clean
```

## Quick Start

```c
#include "exact_cooling_integrator.h"
#include "getMU.h"

int main() {
    // 1. Initialize mean molecular weight grid
    init_mu_grid();
    
    // 2. Load cooling table (from CLOUDY or similar)
    load_cooling_table("coolheat.z=0.bin");
    
    // 3. Perform exact integration
    double nH = 0.1;           // Hydrogen number density [cm^-3]
    double Z = 1.0;            // Metallicity [Z_solar]
    double T = 1.0e6;          // Initial temperature [K]
    double dt = 0.1 * 3.156e13; // Timestep [seconds] (0.1 Myr)
    double T_final;
    
    exact_integrate(nH, Z, T, dt, &T_final);
    printf("T: %.2e -> %.2e K\n", T, T_final);
    
    // 4. Cleanup
    free_cooling_table();
    free_mu_grid();
    return 0;
}
```

## Fast Lookup Table Mode

For repeated calls (e.g., in hydrodynamic simulations), use the pre-computed lookup table:

```c
// One-time setup (can be parallelized with OpenMP)
load_cooling_table("coolheat.z=0.bin");
build_dT_lookup_table(dt);  // Build table for fixed timestep

// Fast repeated calls in simulation loop (~100x faster)
#pragma omp parallel for
for (int i = 0; i < N_particles; i++) {
    double T_final;
    interpolate_exact_integrate(nH[i], Z[i], T[i], &T_final);
    T[i] = T_final;
}

// Cleanup
free_dT_lookup_table();
free_cooling_table();
```

## API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `load_cooling_table(filename)` | Load binary cooling table |
| `exact_integrate(nH, Z, T, dt, &T_final)` | Exact integration (accurate) |
| `interpolate_exact_integrate(nH, Z, T, &T_final)` | Fast interpolation (approximate) |
| `euler_integrate(nH, Z, T, dt, &T_final)` | Simple Euler method (for comparison) |

### Lookup Table Functions

| Function | Description |
|----------|-------------|
| `build_dT_lookup_table(dt)` | Build lookup table for timestep `dt` |
| `free_dT_lookup_table()` | Free lookup table memory |
| `is_dT_table_loaded()` | Check if table is ready |
| `get_dT_table_timestep()` | Get timestep used for table |

### Utility Functions

| Function | Description |
|----------|-------------|
| `get_cooling_time(nH, Z, T)` | Get cooling time $t_{cool} = u/\Lambda$ |
| `get_net_cooling_rate(nH, Z, T)` | Get net cooling rate $\Lambda$ |
| `print_table_info()` | Print cooling table information |

### Mean Molecular Weight (getMU)

| Function | Description |
|----------|-------------|
| `init_mu_grid()` | Initialize # lookup table |
| `free_mu_grid()` | Free # table memory |
| `get_mean_molecular_weight(rho, T, Z)` | Get # from table |
| `calculate_mu(rho, T, Z)` | Direct Saha equation solve |

## Cooling Table Format

The binary cooling table has the following structure:

```c
typedef struct {
    int L, M, N;                    // Grid dimensions
    double redshift;                // Redshift
    double density_min, density_max;
    double metallicity_min, metallicity_max;
    double temp_min, temp_max;
    double temp[N_Temperature_Array]; // Temperature grid from CLOUDY
} GridParams;

// File layout:
// [GridParams]
// [cooling_rate: L*M*N doubles]  // erg/cm^3/s
// [heating_rate: L*M*N doubles]  // erg/cm^3/s
```

## Grid Resolution Recommendations

| Accuracy Target | L (Density) | M (Metallicity) | N (Temperature) | Memory |
|-----------------|-------------|-----------------|-----------------|--------|
| ~5% error | 40 | 30 | 128 | ~1.2 MB |
| **~1-2% error (recommended)** | **60** | **40** | **256** | **~4.7 MB** |
| <1% error | 80 | 50 | 512 | ~15.6 MB |

Current default: `N_Temperature_Array = 256` (defined in `cooling_rate_mpi.h`)

## Thread Safety

| Function | Thread-Safe | Notes |
|----------|-------------|-------|
| `load_cooling_table()` | #  | Call once at startup |
| `build_dT_lookup_table()` | #  | Call once after loading |
| `free_cooling_table()` | #  | Call once at shutdown |
| **`exact_integrate()`** | #  | Safe for parallel use |
| **`interpolate_exact_integrate()`** | #  | Safe for parallel use |
| `get_cooling_time()` | #  | Safe for parallel use |
| `get_net_cooling_rate()` | #  | Safe for parallel use |

## Test Results

Example output with 60×40×256 grid:

```
=== TEST 6: Lookup Table Accuracy ===
Grid point tests: 125, Max error: 1.95e-15 (machine precision)
Intermediate tests: 100, Max error: ~1-50% (depends on grid resolution)

=== TEST 7: Performance ===
exact_integrate():             4.00 #s/call
interpolate_exact_integrate(): ~0.01 #s/call
Speedup: ~400x
```

## How It Works

### Net Cooling Rate

The code computes the **net cooling rate** by combining cooling and heating:

```c
net_rate = cooling_rate - heating_rate
```

| Condition | Meaning | Result |
|-----------|---------|--------|
| `net_rate > 0` | Cooling dominant | Temperature decreases |
| `net_rate < 0` | Heating dominant | Temperature increases |
| `net_rate = 0` | **Equilibrium** | Cooling = Heating |

**Important**: Given infinite time, the gas will converge to the **equilibrium temperature** where `cooling = heating`.

### Integration Methods

#### `exact_integrate()` - Accurate Method

```
+-------------------------------------------------------------+
|  exact_integrate(nH, Z, T_init, dt, &T_final)               |
|                                                             |
|  1. Uses net_rate = cooling - heating (directly)            |
|  2. Performs Y-function integration (Eq. 4-10)              |
|  3. Computes T_final exactly                                |
|                                                             |
|  # Uses net_rate at every call                              |
|  # Slower (~4 #s/call)                                      |
+-------------------------------------------------------------+
```

The Y-function transformation (from ZSH17):

$$Y(u) = \frac{\Lambda_{ref}}{u_{ref}} \int \frac{du}{\Lambda(u)}$$

where $\Lambda(u)$ is the **net cooling rate** (cooling - heating).

#### `interpolate_exact_integrate()` - Fast Method

```
+-------------------------------------------------------------+
|  [Pre-computation] build_dT_lookup_table(dt)                |
|                                                             |
|  For all (nH, Z, T) grid points:                            |
|      Call exact_integrate()  ----# Store in lookup table    |
|      (uses net_rate internally)                             |
+-------------------------------------------------------------+
                          |
                          #
+-------------------------------------------------------------+
|  interpolate_exact_integrate(nH, Z, T_init, &T_final)       |
|                                                             |
|  1. Trilinear interpolation from lookup table               |
|  2. T_final = T_init × 10^(interpolated_log_ratio)          |
|                                                             |
|  # Does NOT compute net_rate directly                       |
|  # Uses pre-computed results (interpolation only)           |
|  # Very fast (~0.01 #s/call)                                |
+-------------------------------------------------------------+
```

### Comparison

| Function | Uses net_rate? | Method | Speed |
|----------|---------------|--------|-------|
| `exact_integrate()` | **Directly** | Y-function integration | ~4 #s |
| `interpolate_exact_integrate()` | **Indirectly** | Table interpolation | ~0.01 #s |

Both functions reflect the same physics (cooling - heating), but differ in computation:
- `exact_integrate()`: Computes exactly every time
- `interpolate_exact_integrate()`: Interpolates pre-computed results

### Equilibrium Convergence

Example test with `nH = 0.1 cm#³`, `Z = 0.1 Z#`, integrated for 1 Gyr:

| Initial T [K] | Final T [K] | Converged to T_eq? |
|---------------|-------------|-------------------|
| 1,000 | 8,263 | # |
| 10,000 | 8,263 | # |
| 100,000 | 8,263 | # |
| 1,000,000 | 8,263 | # |
| 10,000,000 | 8,263 | # |

All initial temperatures converge to the same equilibrium temperature!

## Physical Model

The code solves the coupled Saha ionization equations for H and He:

```
H   + e  #  H#   + 2e    (# = 13.6 eV)
He  + e  #  He#  + 2e    (# = 24.6 eV)
He# + e  #  He## + 2e    (# = 54.4 eV)
```

Mass fractions are computed from metallicity:
- $X = X_{\text{primordial}} \times (1 - Z)$
- $Y = Y_{\text{primordial}} + \Delta Y/\Delta Z \times Z$
- Primordial values: $X_p = 0.76$, $Y_p = 0.24$

## References

1. **Zhu, Smith & Hernquist (2017)**, "An Exact Integration Scheme for Radiative Cooling in Hydrodynamical Simulations", MNRAS 470, 1017 ([ADS](https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.1017Z))

2. **Townsend (2009)**, "An Exact Integration Scheme for Radiative Cooling in Hydrodynamical Simulations", ApJS 181, 391 ([ADS](https://ui.adsabs.harvard.edu/abs/2009ApJS..181..391T))

## License

This library is provided for academic and research use. Please cite the relevant papers when publishing results.

## Author

Prof. Juhan Kim at Korea Institute for Advanced Study (kjhan0606@gmail.com)

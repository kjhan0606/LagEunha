# Stellar Wind Yield Library

A C library for calculating metal yields from stellar winds ejected into the interstellar medium (ISM) for galaxy formation and evolution simulations.

## Features

- **Multiple IMF Models**: Salpeter (1955), Kroupa (2001), Chabrier (2003), Baldry & Glazebrook (2003)
- **Time-dependent Yields**: Calculate cumulative yields for star clusters at any age
- **Two Yield Types**:
  - **Net Yields**: For chemical evolution tracking (can be negative)
  - **Ejected Mass**: Actual mass injected into ISM (always positive)
- **Comprehensive Element Coverage**: H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, and total metals
- **Modern Yield Tables**: Based on latest nucleosynthesis calculations (2016-2023)

## Scientific References

### Yield Tables
- **Kobayashi et al. (2020)** ApJ 900, 179 - Comprehensive nucleosynthesis yields
- **Limongi & Chieffi (2018)** ApJS 237, 13 - Massive star yields with rotation
- **Higgins et al. (2023)** MNRAS 526, 534 - Very Massive Star (VMS) wind yields
- **Karakas & Lugaro (2016)** ApJ 825, 26 - AGB star yields

### Mass-Loss Prescriptions
- **Vink et al. (2001)** A&A 369, 574 - O/B star winds
- **Vink (2022)** ARA&A 60, 203 - Updated mass-loss theory review
- **Sabhahit et al. (2022)** MNRAS 514, 3736 - VMS enhanced winds

## Installation

### Prerequisites
- GCC or compatible C compiler
- GNU Make
- Math library (libm)

### Build

```bash
# Clone or download the source files
# Then build:
make

# This creates:
#   libstellar_wind.a  - Static library
#   test_stellar_wind  - Test/demo program
```

### Install (Optional)

```bash
# Install to /usr/local (requires sudo)
sudo make install

# Or install to custom location
make PREFIX=/your/path install
```

## Quick Start

### Basic Example

```c
#include "stellar_wind_lib.h"
#include <stdio.h>

int main() {
    SWL_ClusterWindYield result;
    
    // Initialize library
    swl_init();
    
    // Calculate wind yields for a star cluster
    swl_calculate_cluster_wind_yield(
        SWL_IMF_CHABRIER,    // IMF model
        0.02,                // Metallicity Z (solar # 0.014)
        1.0e6,               // Cluster mass: 10^6 M_sun
        1.0e7,               // Cluster age: 10 Myr
        &result
    );
    
    // Access results
    printf("Total wind mass ejected: %.4e M_sun\n", result.total_wind_mass);
    
    // For hydrodynamic simulations - actual mass going into ISM
    printf("\nEjected mass to ISM:\n");
    printf("  H:  %.4e M_sun\n", result.ejected_mass.yields[SWL_EL_H]);
    printf("  He: %.4e M_sun\n", result.ejected_mass.yields[SWL_EL_HE]);
    printf("  C:  %.4e M_sun\n", result.ejected_mass.yields[SWL_EL_C]);
    printf("  O:  %.4e M_sun\n", result.ejected_mass.yields[SWL_EL_O]);
    
    // For chemical evolution - net production
    printf("\nNet yields (for chemical evolution):\n");
    printf("  H:  %+.4e M_sun (consumed)\n", result.element_yields.yields[SWL_EL_H]);
    printf("  C:  %+.4e M_sun (produced)\n", result.element_yields.yields[SWL_EL_C]);
    
    // Current ejection rates
    printf("\nCurrent ejection rate: %.4e M_sun/yr\n", result.wind_mass_rate);
    
    swl_cleanup();
    return 0;
}
```

### Compile and Run

```bash
gcc -o my_program my_program.c -L. -lstellar_wind -lm
./my_program
```

## API Reference

### Main Functions

#### `swl_calculate_cluster_wind_yield()`

The primary function for galaxy simulations. Calculates cumulative wind yields for a star cluster.

```c
SWL_ErrorCode swl_calculate_cluster_wind_yield(
    SWL_IMFType imf_type,        // IMF model
    double metallicity,           // Initial Z (0 < Z < 0.1)
    double cluster_mass,          // Total mass in M_sun
    double cluster_age,           // Age in years
    SWL_ClusterWindYield* result  // Output structure
);
```

**Returns**: `SWL_SUCCESS` on success, error code otherwise.

#### `swl_calculate_star_wind_yield()`

Calculate wind yields for a single star.

```c
SWL_ErrorCode swl_calculate_star_wind_yield(
    double initial_mass,          // Stellar mass in M_sun
    double metallicity,           // Initial Z
    double stellar_age,           // Age in years (-1 for lifetime total)
    SWL_StarWindYield* result     // Output structure
);
```

### Result Structures

#### `SWL_ClusterWindYield`

```c
typedef struct {
    // Input parameters
    SWL_IMFType imf_type;
    double metallicity;
    double cluster_mass;
    double cluster_age;
    
    // Yield results - TWO DIFFERENT QUANTITIES
    SWL_ElementYields element_yields;  // NET yields (Ejected - Initial)
    SWL_ElementYields ejected_mass;    // ACTUAL mass ejected to ISM
    double total_wind_mass;            // Total wind mass (M_sun)
    
    // Mass budget
    double mass_in_living_stars;       // Mass still in stars
    double mass_in_remnants;           // Mass in WD/NS/BH
    double mass_returned_total;        // Total returned to ISM
    
    // Stellar population info
    double turnoff_mass;               // Current MS turnoff mass
    double n_stars_total;              // Total stars formed
    double n_stars_dead;               // Stars that have died
    
    // Instantaneous rates
    double wind_mass_rate;             // Current wind rate (M_sun/yr)
    SWL_ElementYields yield_rate;      // Net yield rate
    SWL_ElementYields ejection_rate;   // Ejection rate to ISM
} SWL_ClusterWindYield;
```

### Element Indices

```c
typedef enum {
    SWL_EL_H = 0,      // Hydrogen
    SWL_EL_HE,         // Helium
    SWL_EL_C,          // Carbon
    SWL_EL_N,          // Nitrogen
    SWL_EL_O,          // Oxygen
    SWL_EL_NE,         // Neon
    SWL_EL_MG,         // Magnesium
    SWL_EL_SI,         // Silicon
    SWL_EL_S,          // Sulfur
    SWL_EL_CA,         // Calcium
    SWL_EL_FE,         // Iron
    SWL_EL_Z_TOTAL     // Total metals
} SWL_Element;
```

### IMF Types

```c
typedef enum {
    SWL_IMF_SALPETER = 0,        // Salpeter (1955): #(m) # m^-2.35
    SWL_IMF_KROUPA,              // Kroupa (2001): 3-segment power law
    SWL_IMF_CHABRIER,            // Chabrier (2003): log-normal + power law
    SWL_IMF_BALDRY_GLAZEBROOK    // Baldry & Glazebrook (2003)
} SWL_IMFType;
```

### Utility Functions

```c
// Stellar lifetime (years)
double swl_stellar_lifetime(double mass, double metallicity);

// Main sequence turnoff mass for given age
double swl_turnoff_mass(double age, double metallicity);

// Remnant mass (WD/NS/BH)
double swl_remnant_mass(double initial_mass, double metallicity);

// Normalized IMF value
double swl_imf(double mass, SWL_IMFType imf_type);

// Solar abundance (mass fraction)
double swl_get_solar_abundance(SWL_Element element);

// Export results to CSV
SWL_ErrorCode swl_export_to_csv(const SWL_ClusterWindYield* result, 
                                 const char* filename);
```

## Understanding Net Yield vs Ejected Mass

### Why Two Different Quantities?

| Quantity | Definition | H value | Use Case |
|----------|------------|---------|----------|
| **Net Yield** | Ejected - Initial | Negative | Chemical evolution models |
| **Ejected Mass** | Actual mass to ISM | Positive | Hydrodynamic simulations |

### Physical Explanation

Stars convert hydrogen into heavier elements through nuclear fusion:
```
4H # He + energy
```

When a star ejects material via stellar winds:
- The ejected gas has **less H** and **more He/metals** than the initial composition
- **Net Yield of H is negative** because H was consumed as fuel
- **Ejected H mass is positive** because H is still present in the wind (just less than initially)

### Example (10# M_sun cluster at 10 Myr)

| Element | Net Yield | Ejected Mass | Initial in Wind |
|---------|-----------|--------------|-----------------|
| H | -6.28×10# | +3.73×10# | 1.00×10# M_sun |
| He | +3.84×10# | +7.74×10# | 3.91×10# M_sun |
| C | +5.32×10³ | +5.82×10³ | 4.95×10² M_sun |
| O | +7.69×10³ | +8.89×10³ | 1.21×10³ M_sun |

**Wind composition**: H 27% + He 55% + metals 20%  
(vs initial ISM: H 74% + He 24% + metals 2%)

## Example Results

### Time Evolution of Yields

For a 10# M_sun cluster with Chabrier IMF at solar metallicity:

| Age | Turnoff Mass | Wind Mass | C yield | N yield | O yield |
|-----|--------------|-----------|---------|---------|---------|
| 1 Myr | 13.1 M_sun | 1.10×10# | 3.84×10³ | 5.04×10³ | 6.46×10³ |
| 10 Myr | 6.4 M_sun | 1.39×10# | 5.32×10³ | 7.31×10³ | 7.69×10³ |
| 100 Myr | 3.3 M_sun | 2.04×10# | 9.75×10³ | 1.01×10# | 9.14×10³ |
| 1 Gyr | 1.8 M_sun | 2.72×10# | 1.17×10# | 1.10×10# | 9.54×10³ |
| 10 Gyr | 1.1 M_sun | 3.31×10# | 1.23×10# | 1.12×10# | 9.61×10³ |

### IMF Comparison

Net yields at 10 Myr for 10# M_sun cluster (Z = 0.02):

| Element | Salpeter | Kroupa | Chabrier | Baldry-Glazebrook |
|---------|----------|--------|----------|-------------------|
| C | 3.04×10³ | 5.02×10³ | 5.32×10³ | 4.14×10³ |
| N | 4.18×10³ | 6.89×10³ | 7.31×10³ | 5.70×10³ |
| O | 4.29×10³ | 7.25×10³ | 7.69×10³ | 5.86×10³ |
| Z_total | 1.44×10# | 2.40×10# | 2.54×10# | 1.96×10# |

## Running Tests

```bash
# Run full test suite
make test

# Run with custom parameters
./test_stellar_wind [imf] [Z] [mass] [age]

# Examples:
./test_stellar_wind 2 0.02 1e6 1e7    # Chabrier, solar Z, 10^6 Msun, 10 Myr
./test_stellar_wind 1 0.001 1e5 1e8   # Kroupa, low Z, 10^5 Msun, 100 Myr
```

## Integration with Simulation Codes

### For SPH/Grid-based Hydro Codes

Use `ejected_mass` for mass injection into ISM:

```c
// In your feedback routine
void inject_stellar_wind_feedback(Particle* star_particle, double dt) {
    SWL_ClusterWindYield result;
    swl_calculate_cluster_wind_yield(
        SWL_IMF_CHABRIER,
        star_particle->metallicity,
        star_particle->mass,
        star_particle->age,
        &result
    );
    
    // Mass to inject this timestep
    double dm_H  = result.ejection_rate.yields[SWL_EL_H] * dt;
    double dm_He = result.ejection_rate.yields[SWL_EL_HE] * dt;
    double dm_C  = result.ejection_rate.yields[SWL_EL_C] * dt;
    // ... etc
    
    // Inject into neighboring gas cells/particles
    inject_mass_to_neighbors(star_particle, dm_H, dm_He, dm_C, ...);
}
```

### For Chemical Evolution Models

Use `element_yields` (net yields) for tracking metallicity evolution:

```c
// Update ISM metallicity
void update_ism_composition(Galaxy* gal, double dt) {
    for (int i = 0; i < gal->n_star_particles; i++) {
        SWL_ClusterWindYield result;
        swl_calculate_cluster_wind_yield(...);
        
        // Net production adds to ISM metal content
        gal->ism_carbon += result.yield_rate.yields[SWL_EL_C] * dt;
        gal->ism_oxygen += result.yield_rate.yields[SWL_EL_O] * dt;
        // ...
    }
}
```

## File Structure

```
stellar_wind_lib/
+-- stellar_wind_lib.h    # Header file with API declarations
+-- stellar_wind_lib.c    # Implementation
+-- test_stellar_wind.c   # Test and demonstration program
+-- Makefile              # Build system
+-- README.md             # This file
+-- libstellar_wind.a     # Compiled static library
```

## License

This library is provided for academic and research use. Please cite the relevant yield table papers when publishing results.

## Version History

- **v2.0** (2025): Added ejected mass calculation, separate from net yields
- **v1.0** (2025): Initial release with net yield calculation

## Contact

For questions, bug reports, or feature requests, please open an issue on the repository.

---

*This library was developed for galaxy formation simulations requiring accurate treatment of stellar wind feedback and chemical enrichment.*

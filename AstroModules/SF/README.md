# Star Formation Module for Galaxy Simulations

A comprehensive C implementation of modern star formation models and IMF sampling for cosmological galaxy formation simulations.

## Overview

This module implements state-of-the-art star formation prescriptions based on recent observational and theoretical work. It is designed to be integrated into galaxy formation simulation codes and includes:

- Multiple star formation rate calculation methods
- Accurate IMF sampling (Chabrier, Kroupa, Salpeter)
- Physical dependencies on gas properties (density, temperature, metallicity, turbulence)
- Individual stellar mass sampling for population synthesis

## Key Features

### Star Formation Models

1. **Kennicutt-Schmidt Law** (Kennicutt 1998)
   - Empirical relation: Σ_SFR ∝ Σ_gas^1.4
   - Global scaling law observed in galaxies
   - Accounts for surface density effects

2. **Volumetric Model** (Krumholz & McKee 2005)
   - SFR = ε_ff × M_gas / t_ff
   - Free-fall time based approach
   - Turbulence-regulated efficiency (Padoan et al. 2012)

### Initial Mass Function (IMF)

Three widely-used IMF models are implemented:

1. **Chabrier (2003) IMF**
   - Log-normal for m < 1 M_sun
   - Power law (α = -2.3) for m ≥ 1 M_sun
   - Most commonly used in modern simulations

2. **Kroupa (2001) IMF**
   - Broken power law with three segments
   - α = -0.3 for m < 0.08 M_sun
   - α = -1.3 for 0.08 < m < 0.5 M_sun
   - α = -2.3 for m > 0.5 M_sun

3. **Salpeter (1955) IMF**
   - Simple power law: ξ(m) ∝ m^-2.35
   - Historical reference

## Physics Included

### Input Parameters
- **Gas density** (ρ): Volume density in g cm^-3
- **Temperature** (T): Gas temperature in K
- **Metallicity** (Z): Metal mass fraction
- **Velocity dispersion** (σ_v): 3D turbulent velocity in cm s^-1
- **Volume**: Gas element volume in cm^3

### Calculated Quantities
- Free-fall time: t_ff = sqrt(3π / (32Gρ))
- Jeans mass: M_J ∝ (T^3 / ρ)^1/2
- Mach number: M = σ_v / c_s
- Star formation efficiency: ε_ff(M, α_vir)

### Physical Thresholds
- Minimum density: ~100 H atoms/cm^3
- Maximum temperature: 10^4 K
- Only cold, dense gas forms stars

## Usage

### Basic Example

```c
#include "star_formation.h"

int main() {
    // Define gas element
    GasElement gas;
    gas.density = 1.0e-20;              // g cm^-3
    gas.temperature = 20.0;             // K
    gas.metallicity = METALLICITY_SOLAR;
    gas.velocity_dispersion = 1.0e5;    // cm s^-1
    gas.volume = pow(10.0 * PARSEC, 3); // cm^3
    
    // Form stars
    StarParticle sp;
    double timestep = 1.0e6; // years
    
    if (form_stars(&gas, timestep, IMF_CHABRIER_2003, &sp)) {
        printf("Formed %.2e M_sun in %d stars\n", 
               sp.mass, sp.num_stars);
        
        // Access individual stellar masses
        for (int i = 0; i < sp.num_stars; i++) {
            printf("Star %d: %.3f M_sun\n", 
                   i, sp.stellar_masses[i]);
        }
        
        // Clean up
        free_star_particle(&sp);
    }
    
    return 0;
}
```

### Compilation

```bash
# Compile with math library
gcc -o star_formation_example star_formation_example.c -lm -O2

# Or use the provided Makefile
make
make run
```

## Scientific References

### Star Formation Theory
1. **Kennicutt (1998)**, "The Global Schmidt Law in Star-forming Galaxies"
   - ApJ 498, 541
   - DOI: 10.1086/305588
   - Empirical Kennicutt-Schmidt relation

2. **Krumholz & McKee (2005)**, "A General Theory of Turbulence-regulated Star Formation"
   - ApJ 630, 250
   - DOI: 10.1086/431734
   - Volumetric SF efficiency model

3. **Padoan et al. (2012)**, "The Star Formation Rate of Supersonic MHD Turbulence"
   - ApJ 759, L27
   - DOI: 10.1088/2041-8205/759/2/L27
   - Turbulence effects on SF efficiency

4. **Ostriker & Kim (2022)**, "Pressure-Regulated, Feedback-Modulated Star Formation"
   - ApJ 936, 137
   - DOI: 10.3847/1538-4357/ac7de2
   - PRFM model from TIGRESS simulations

5. **Kim et al. (2024)**, "Star Formation Model from High-Resolution ISM Simulations"
   - arXiv:2502.13244
   - Latest TIGRESS-based implementations

### IMF Studies
6. **Salpeter (1955)**, "The Luminosity Function and Stellar Evolution"
   - ApJ 121, 161
   - DOI: 10.1086/145971
   - Original power-law IMF

7. **Kroupa (2001)**, "On the Variation of the Initial Mass Function"
   - MNRAS 322, 231
   - DOI: 10.1046/j.1365-8711.2001.04022.x
   - Multi-component broken power law

8. **Chabrier (2003)**, "Galactic Stellar and Substellar Initial Mass Function"
   - PASP 115, 763
   - DOI: 10.1086/376392
   - Log-normal + power law formulation

9. **Kroupa et al. (2024)**, "The Initial Mass Function of Stars"
   - arXiv:2410.07311
   - Recent comprehensive review

### Observational Constraints
10. **Kennicutt & Evans (2012)**, "Star Formation in the Milky Way and Nearby Galaxies"
    - ARA&A 50, 531
    - DOI: 10.1146/annurev-astro-081811-125610
    - Review of SF observations

## Implementation Notes

### Star Formation Rate Calculation
The module provides two methods:

1. **Kennicutt-Schmidt approach**: Better for resolved simulations where surface densities are well-defined
2. **Volumetric approach**: Better for SPH or grid-based codes with volume elements

Users can switch between methods in `form_stars()` function.

### IMF Sampling
- Uses inverse transform sampling for power laws
- Box-Muller transform for log-normal distributions
- Mass range: 0.01 - 100 M_sun
- Normalizes sampled masses to match total stellar mass exactly

### Numerical Considerations
- Maximum 10,000 stars sampled per gas element (configurable)
- For large stellar masses, consider using stellar population synthesis instead
- Velocity inheritance includes Gaussian random perturbations
- Maximum 50% gas-to-star conversion per timestep for stability

## Testing

The example program includes four test cases:

1. **Dense Molecular Cloud**: Realistic GMC conditions
2. **IMF Comparison**: Side-by-side comparison of all three IMFs
3. **Parameter Effects**: SFR dependence on metallicity and temperature
4. **Time Evolution**: Following gas depletion over time

Run tests with:
```bash
make run
```

## Future Enhancements

Potential additions for future versions:
- Magnetic field effects (Kretschmer & Teyssier 2020)
- Metallicity-dependent IMF variations
- Stellar feedback implementation
- Time-dependent stellar evolution
- Binary star populations
- Clustered star formation

## License

This code is provided for educational and research purposes.
Please cite the relevant papers when using this module in publications.

## Author

Created for galaxy formation simulation research.
Based on observational and theoretical work from 1955-2024.

## Version History

- v1.0 (2025): Initial implementation
  - Kennicutt-Schmidt and volumetric SF models
  - Three IMF options
  - Complete documentation

/**
 * =============================================================================
 * Mean Molecular Weight Calculator for Interstellar Gas
 * Header File
 * =============================================================================
 * 
 * This library calculates the mean molecular weight (mu) of interstellar gas
 * by solving the coupled Saha ionization equations for H and He simultaneously.
 * 
 * The coupled equations are reduced to a 4th-order polynomial and solved
 * using the bisection method for numerical stability.
 * 
 * Optimized for galaxy formation simulations with density range:
 *   10^-31 to 10^-18 g/cm^3 (IGM to star-forming cores)
 * 
 * Usage:
 *   1. Call init_mu_grid() to initialize the lookup tables
 *   2. Call get_mean_molecular_weight(rho, T, Z_metal) for interpolated values
 *      or calculate_mu(rho, T, Z_metal) for direct calculation
 *   3. Call free_mu_grid() when done to release memory
 * 
 * =============================================================================
 */

#ifndef MEAN_MOLECULAR_WEIGHT_H
#define MEAN_MOLECULAR_WEIGHT_H

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 *                         PHYSICAL CONSTANTS (CGS)
 * ============================================================================ */

#define MMW_C_LIGHT       2.99792458e10    /* Speed of light [cm/s] */
#define MMW_K_BOLTZMANN   1.380649e-16     /* Boltzmann constant [erg/K] */
#define MMW_H_PLANCK      6.62607015e-27   /* Planck constant [erg*s] */
#define MMW_M_ELECTRON    9.1093837e-28    /* Electron mass [g] */
#define MMW_M_HYDROGEN    1.6735575e-24    /* Hydrogen atom mass [g] */
#define MMW_EV_TO_ERG     1.602176634e-12  /* 1 eV in erg */

/* Ionization energies [eV] */
#define MMW_CHI_H_I_EV    13.598           /* H ionization */
#define MMW_CHI_HE_I_EV   24.587           /* He I ionization */
#define MMW_CHI_HE_II_EV  54.418           /* He II ionization */

/* Primordial abundances (mass fraction) */
#define MMW_X_PRIMORDIAL  0.76             /* Hydrogen mass fraction */
#define MMW_Y_PRIMORDIAL  0.24             /* Helium mass fraction */
#define MMW_Z_SOLAR       0.02             /* Solar metallicity */

/* Grid bounds (optimized for galaxy formation simulations) */
#define MMW_LOG_RHO_MIN   -31.0            /* log10(rho_min) [g/cm^3]: IGM */
#define MMW_LOG_RHO_MAX   -18.0            /* log10(rho_max) [g/cm^3]: SF cores */
#define MMW_TEMP_MIN       1000.0          /* Minimum temperature [K] */
#define MMW_TEMP_MAX       10000.0         /* Maximum temperature [K] */
#define MMW_Z_MIN          0.0             /* Minimum metallicity [Z_solar] */
#define MMW_Z_MAX          3.0             /* Maximum metallicity [Z_solar] */

/* ============================================================================
 *                           GRID FUNCTIONS
 * ============================================================================ */

/**
 * Initialize the mean molecular weight lookup tables
 * 
 * This pre-computes mu values on a 3D grid of (density, temperature, metallicity)
 * and builds an inverse lookup table T(T/mu, rho, Z).
 * 
 * Density range is optimized for galaxy formation simulations:
 *   - IGM:  10^-31 to 10^-28 g/cm^3
 *   - CGM:  10^-28 to 10^-26 g/cm^3
 *   - ISM:  10^-26 to 10^-24 g/cm^3
 *   - GMC:  10^-24 to 10^-21 g/cm^3
 *   - SF:   10^-21 to 10^-18 g/cm^3
 * 
 * Must be called before using get_mean_molecular_weight() or 
 * get_temperature_from_Tmu().
 * 
 * @return 0 on success, -1 on failure
 */
int init_mu_grid(void);

/**
 * Free all memory allocated for the lookup tables
 */
void free_mu_grid(void);

/**
 * Print information about the current grid
 */
void print_grid_info(void);

/**
 * Save the grid to a binary file
 * 
 * @param filename  Output filename
 * @return          0 on success, -1 on failure
 */
int save_grid_to_file(const char *filename);

/**
 * Load the grid from a binary file
 * 
 * @param filename  Input filename
 * @return          0 on success, -1 on failure
 */
int load_grid_from_file(const char *filename);

/* ============================================================================
 *                    MEAN MOLECULAR WEIGHT FUNCTIONS
 * ============================================================================ */

/**
 * Calculate mean molecular weight by solving coupled Saha equations
 * 
 * This function directly solves the 4th-order polynomial derived from the
 * coupled H + He Saha equations using bisection.
 * 
 * The metallicity affects the H and He mass fractions:
 *   X = X_primordial * (1 - Z_metal * Z_solar)
 *   Y = Y_primordial * (1 - Z_metal * Z_solar)
 * 
 * @param rho       Physical density [g/cm^3]
 * @param T         Temperature [K]
 * @param Z_metal   Metallicity in solar units (Z_solar = 0.02)
 * @return          Mean molecular weight [m_H units]
 */
double calculate_mu(double rho, double T, double Z_metal);

/**
 * Get mean molecular weight from pre-computed lookup table
 * 
 * This is faster than calculate_mu() but requires init_mu_grid() to be called.
 * Uses trilinear interpolation on the pre-computed grid.
 * 
 * @param rho       Physical density [g/cm^3]
 * @param T         Temperature [K]
 * @param Z_metal   Metallicity in solar units (Z_solar ~ 0.02)
 * @return          Mean molecular weight [m_H units]
 */
double get_mean_molecular_weight(double rho, double T, double Z_metal);

/**
 * Get temperature from T/mu and density (inverse lookup)
 * 
 * This is useful when internal energy is known but temperature is not:
 *   u = (3/2) * k_B * T / (mu * m_H)
 *   T/mu = 2 * u * m_H / (3 * k_B)
 * 
 * @param Tmu       Temperature divided by mean molecular weight [K]
 * @param rho       Physical density [g/cm^3]
 * @param Z_metal   Metallicity in solar units
 * @return          Temperature [K]
 */
double get_temperature_from_Tmu(double Tmu, double rho, double Z_metal);

#ifdef __cplusplus
}
#endif

#endif /* MEAN_MOLECULAR_WEIGHT_H */

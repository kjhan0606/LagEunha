/*
 * Stellar Cluster Supernova Library
 * 
 * Calculates supernova rates, energy release, mass ejection, and
 * chemical yields for stellar clusters based on:
 * - Cluster mass and age
 * - Initial Mass Function (IMF)
 * - Metallicity
 * 
 * Author: Based on Keegans+ 2023, Limongi & Chieffi 2018
 */

#ifndef CLUSTER_SUPERNOVA_H
#define CLUSTER_SUPERNOVA_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Physical constants */
#define SOLAR_MASS 1.989e33      // grams
#define YEAR_IN_SECONDS 3.156e7  // seconds
#define FOE 1.0e51               // 10^51 erg (typical SN energy)

/* Stellar mass limits */
#define MIN_MASS_STARS 0.08      // Minimum mass for hydrogen burning (M_sun)
#define MAX_MASS_STARS 120.0     // Maximum stellar mass (M_sun) - normal Pop I/II
#define MIN_MASS_SNII 8.0        // Minimum mass for core-collapse SN (M_sun)
#define MAX_MASS_SNII 40.0       // Maximum mass for successful SN (beyond this: direct BH)
#define CHANDRASEKHAR_MASS 1.4   // White dwarf mass for Type Ia (M_sun)

/* Pop III / PISN limits */
#define MIN_MASS_PISN 140.0      // Minimum mass for pair-instability SN (M_sun)
#define MAX_MASS_PISN 260.0      // Maximum mass for PISN (beyond this: direct BH)
#define MAX_MASS_POPIII 300.0    // Maximum Pop III stellar mass (M_sun)
#define N_PISN_MODELS 5          // Number of PISN yield models

/* Number of elements tracked */
#define N_ELEMENTS 11

/* Element enumeration */
typedef enum {
    H_ELEM = 0,
    HE_ELEM,
    C_ELEM,
    N_ELEM,
    O_ELEM,
    NE_ELEM,
    MG_ELEM,
    SI_ELEM,
    S_ELEM,
    CA_ELEM,
    FE_ELEM
} ElementType;

/* IMF type */
typedef enum {
    IMF_SALPETER,         // Salpeter (1955): dN/dM # M^(-2.35)
    IMF_KROUPA,           // Kroupa (2001): Broken power law
    IMF_CHABRIER,         // Chabrier (2003): Log-normal + power law
    IMF_POPIII_TOPHEAVY   // Pop III top-heavy (Z=0): Peak ~100 M_sun, M^(-1.0)
} IMFType;

/* Star formation mode */
typedef enum {
    SF_INSTANTANEOUS,    // All stars formed at once (single burst)
    SF_CONTINUOUS        // Constant star formation rate
} StarFormationMode;

/* Pop III mode */
typedef enum {
    POPIII_DISABLED = 0,     // No Pop III stars
    POPIII_ENABLED = 1       // Include Pop III stars (PISN)
} PopIIIMode;

/* Structure for cluster properties */
typedef struct {
    double total_mass;              // Total cluster mass (M_sun)
    double age;                     // Cluster age (years)
    double metallicity;             // Metallicity Z (mass fraction)
    IMFType imf_type;               // Type of IMF to use
    StarFormationMode sf_mode;      // Star formation mode
    double star_formation_rate;     // SFR (M_sun/year) - only for SF_CONTINUOUS
    
    /* Pop III options (thread-safe: read-only per calculation) */
    PopIIIMode popiii_mode;         // Enable/disable Pop III stars
    double popiii_imf_max;          // Max mass for Pop III IMF (default: 260 M_sun)
    double popiii_mass_fraction;    // Fraction of mass in Pop III (default: 0.10 = 10%)
} ClusterProperties;

/* Structure for supernova rates */
typedef struct {
    double sn_ia_rate;      // Type Ia SN rate (per year)
    double sn_ii_rate;      // Type II SN rate (per year)
    double pisn_rate;       // Pair-instability SN rate (per year) - Pop III
    double total_rate;      // Total SN rate (per year)
} SNRates;

/* Structure for energy release */
typedef struct {
    double sn_ia_energy;    // Energy from Type Ia (erg/year)
    double sn_ii_energy;    // Energy from Type II (erg/year)
    double pisn_energy;     // Energy from PISN (erg/year) - Pop III
    double total_energy;    // Total energy (erg/year)
} EnergyRelease;

/* Structure for mass ejection */
typedef struct {
    double sn_ia_mass;      // Mass from Type Ia (M_sun/year)
    double sn_ii_mass;      // Mass from Type II (M_sun/year)
    double pisn_mass;       // Mass from PISN (M_sun/year) - Pop III
    double total_mass;      // Total ejected mass (M_sun/year)
} MassEjection;

/* Structure for elemental yields */
typedef struct {
    double element_yields[N_ELEMENTS];  // Yields per element (M_sun/year)
    double total_metals;                 // Total metals (M_sun/year)
    double mean_metallicity;             // Mean Z of ejected material
} ElementalYields;

/* Structure for complete cluster output */
typedef struct {
    SNRates rates;
    EnergyRelease energy;
    MassEjection mass;
    ElementalYields yields;
} ClusterOutput;

/* ===== Core Functions ===== */

/**
 * Calculate all supernova properties for a stellar cluster
 * 
 * @param cluster: Cluster properties (mass, age, Z, IMF)
 * @param output: Pointer to output structure (will be filled)
 * @return: 0 on success, -1 on error
 */
int calculate_cluster_supernovae(const ClusterProperties *cluster, 
                                  ClusterOutput *output);

/**
 * Calculate supernova rates
 */
int calculate_sn_rates(const ClusterProperties *cluster, SNRates *rates);

/**
 * Calculate energy release
 */
int calculate_energy_release(const SNRates *rates, EnergyRelease *energy);

/**
 * Calculate mass ejection
 */
int calculate_mass_ejection(const ClusterProperties *cluster, 
                             const SNRates *rates, 
                             MassEjection *mass);

/**
 * Calculate elemental yields
 */
int calculate_elemental_yields(const ClusterProperties *cluster,
                                const SNRates *rates,
                                ElementalYields *yields);

/* ===== Helper Functions ===== */

/**
 * Calculate stellar lifetime as function of mass and metallicity
 * Based on Hurley et al. (2000)
 */
double stellar_lifetime(double mass, double metallicity);

/**
 * Evaluate IMF at given mass
 */
double evaluate_imf(double mass, IMFType imf_type);

/**
 * Calculate IMF normalization constant
 */
double calculate_imf_normalization(IMFType imf_type);

/**
 * Type Ia delay time distribution (DTD)
 * Based on Maoz & Graur (2017): dN/dt # t^(-1)
 */
double type_ia_dtd(double time);

/**
 * Get element name from enum
 */
const char* get_element_name(ElementType elem);

/**
 * Get metallicity-dependent yield for specific element
 * (Uses interpolation from our yield tables)
 */
double get_element_yield_value(ElementType elem, 
                                int sn_type,  // 0=Type Ia, 1=Type II
                                double mass,  // Progenitor mass
                                double metallicity);

/* ===== Utility Functions ===== */

/**
 * Print cluster output in formatted table
 */
void print_cluster_output(const ClusterOutput *output);

/**
 * Print cluster properties
 */
void print_cluster_properties(const ClusterProperties *cluster);

/**
 * Validate cluster properties
 */
int validate_cluster_properties(const ClusterProperties *cluster);

#endif /* CLUSTER_SUPERNOVA_H */

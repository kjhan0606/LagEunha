/*
 * ============================================================================
 * STELLAR WIND YIELD LIBRARY
 * ============================================================================
 * 
 * A library for calculating metal yields from stellar winds ejected into the
 * interstellar medium (ISM) for galaxy formation simulations.
 * 
 * This library provides functions to calculate:
 * - IMF-integrated yields for stellar populations
 * - Time-dependent yields for star clusters of given mass and age
 * - Mass-dependent yields for individual stars
 * 
 * Yield tables based on:
 *   - Kobayashi et al. (2020) ApJ 900, 179
 *   - Limongi & Chieffi (2018) ApJS 237, 13
 *   - Higgins et al. (2023) MNRAS 526, 534 (VMS winds)
 *   - Karakas & Lugaro (2016) ApJ 825, 26 (AGB yields)
 * 
 * Mass-loss prescriptions:
 *   - Vink et al. (2001, 2022)
 *   - Sabhahit et al. (2022) for VMS
 * 
 * Author: Generated for galaxy simulation research
 * Version: 2.0
 * Date: 2025
 * 
 * ============================================================================
 */

#ifndef STELLAR_WIND_LIB_H
#define STELLAR_WIND_LIB_H

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * PHYSICAL CONSTANTS
 * ============================================================================ */

#define SWL_M_SUN       1.989e33      /* Solar mass in grams */
#define SWL_L_SUN       3.828e33      /* Solar luminosity in erg/s */
#define SWL_Z_SUN       0.0134        /* Solar metallicity (Asplund et al. 2009) */
#define SWL_PI          3.14159265359

/* Mass limits for stellar populations */
#define SWL_M_MIN       0.08          /* Minimum stellar mass (M_sun) */
#define SWL_M_MAX       300.0         /* Maximum stellar mass (M_sun) */
#define SWL_M_AGB_MIN   0.8           /* Minimum AGB progenitor mass (M_sun) */
#define SWL_M_AGB_MAX   8.0           /* Maximum AGB progenitor mass (M_sun) */
#define SWL_M_MASSIVE   8.0           /* Minimum massive star mass (M_sun) */
#define SWL_M_VMS       100.0         /* Very Massive Star threshold (M_sun) */

/* Number of elements tracked */
#define SWL_N_ELEMENTS  12

/* ============================================================================
 * ENUMERATIONS
 * ============================================================================ */

/**
 * @brief Element indices for yield arrays
 */
typedef enum {
    SWL_EL_H = 0,      /**< Hydrogen */
    SWL_EL_HE,         /**< Helium */
    SWL_EL_C,          /**< Carbon */
    SWL_EL_N,          /**< Nitrogen */
    SWL_EL_O,          /**< Oxygen */
    SWL_EL_NE,         /**< Neon */
    SWL_EL_MG,         /**< Magnesium */
    SWL_EL_SI,         /**< Silicon */
    SWL_EL_S,          /**< Sulfur */
    SWL_EL_CA,         /**< Calcium */
    SWL_EL_FE,         /**< Iron */
    SWL_EL_Z_TOTAL     /**< Total metals */
} SWL_Element;

/**
 * @brief Available Initial Mass Function (IMF) models
 */
typedef enum {
    SWL_IMF_SALPETER = 0,        /**< Salpeter (1955): phi(m) ~ m^-2.35 */
    SWL_IMF_KROUPA,              /**< Kroupa (2001): three-segment power law */
    SWL_IMF_CHABRIER,            /**< Chabrier (2003): log-normal + power law */
    SWL_IMF_BALDRY_GLAZEBROOK,   /**< Baldry & Glazebrook (2003) */
    SWL_N_IMF_TYPES
} SWL_IMFType;

/**
 * @brief Stellar evolutionary phases
 */
typedef enum {
    SWL_PHASE_MS = 0,     /**< Main Sequence */
    SWL_PHASE_RGB,        /**< Red Giant Branch */
    SWL_PHASE_AGB,        /**< Asymptotic Giant Branch */
    SWL_PHASE_WR,         /**< Wolf-Rayet */
    SWL_PHASE_REMNANT     /**< Stellar remnant (WD, NS, BH) */
} SWL_StellarPhase;

/**
 * @brief Error codes returned by library functions
 */
typedef enum {
    SWL_SUCCESS = 0,
    SWL_ERR_INVALID_IMF,
    SWL_ERR_INVALID_MASS,
    SWL_ERR_INVALID_METALLICITY,
    SWL_ERR_INVALID_AGE,
    SWL_ERR_NULL_POINTER,
    SWL_ERR_INTEGRATION_FAILED
} SWL_ErrorCode;

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================ */

/**
 * @brief Structure holding yields for individual elements
 * 
 * All yields are in units of solar masses (M_sun).
 * Negative values indicate net consumption of that element.
 */
typedef struct {
    double yields[SWL_N_ELEMENTS];   /**< Net yield for each element (M_sun) */
} SWL_ElementYields;

/**
 * @brief Result structure for single star wind yields
 */
typedef struct {
    double initial_mass;             /**< Initial stellar mass (M_sun) */
    double metallicity;              /**< Initial metallicity Z */
    SWL_ElementYields element_yields; /**< Yields for each element */
    double total_wind_mass;          /**< Total mass lost in winds (M_sun) */
    double wind_velocity;            /**< Terminal wind velocity (km/s) */
} SWL_StarWindYield;

/**
 * @brief Result structure for time-integrated cluster wind yields
 * 
 * This structure contains the cumulative yields from stellar winds
 * for a star cluster of given mass and age.
 * 
 * IMPORTANT: Two types of yields are provided:
 * 
 * 1. net_yields (element_yields): 
 *    Net production = Ejected - Initial
 *    Can be NEGATIVE for H (consumed as fuel)
 *    Use for: Chemical evolution, metallicity tracking
 * 
 * 2. ejected_mass:
 *    Actual mass ejected into ISM via stellar winds
 *    Always POSITIVE
 *    Use for: ISM enrichment, feedback calculations, hydrodynamic simulations
 */
typedef struct {
    /* Input parameters */
    SWL_IMFType imf_type;            /**< IMF model used */
    double metallicity;              /**< Initial metallicity Z */
    double cluster_mass;             /**< Initial cluster mass (M_sun) */
    double cluster_age;              /**< Cluster age (years) */
    
    /* Yield results - TWO DIFFERENT QUANTITIES */
    SWL_ElementYields element_yields; /**< NET yields: Ejected - Initial (M_sun) */
    SWL_ElementYields ejected_mass;   /**< ACTUAL mass ejected to ISM (M_sun) */
    double total_wind_mass;          /**< Total mass ejected in winds (M_sun) */
    
    /* Mass budget */
    double mass_in_living_stars;     /**< Mass still in living stars (M_sun) */
    double mass_in_remnants;         /**< Mass locked in remnants (M_sun) */
    double mass_returned_total;      /**< Total mass returned to ISM (M_sun) */
    
    /* Stellar population info */
    double turnoff_mass;             /**< Current main sequence turnoff mass (M_sun) */
    double n_stars_total;            /**< Total number of stars formed */
    double n_stars_dead;             /**< Number of stars that have died */
    
    /* Rates (instantaneous at cluster_age) */
    double wind_mass_rate;           /**< Current wind mass-loss rate (M_sun/yr) */
    SWL_ElementYields yield_rate;    /**< Current NET yield rate (M_sun/yr) */
    SWL_ElementYields ejection_rate; /**< Current ejection rate to ISM (M_sun/yr) */
    
} SWL_ClusterWindYield;

/**
 * @brief Configuration options for yield calculations
 */
typedef struct {
    double m_min;                    /**< Minimum stellar mass to consider (M_sun) */
    double m_max;                    /**< Maximum stellar mass to consider (M_sun) */
    int n_integration_steps;         /**< Number of integration steps */
    int include_agb_winds;           /**< Include AGB star winds (0/1) */
    int include_massive_winds;       /**< Include massive star winds (0/1) */
    int include_vms_enhancement;     /**< Include VMS wind enhancement (0/1) */
    double metallicity_floor;        /**< Minimum metallicity for scaling */
} SWL_Config;

/* ============================================================================
 * LIBRARY INITIALIZATION AND CONFIGURATION
 * ============================================================================ */

/**
 * @brief Initialize the stellar wind library
 * 
 * Must be called before using any other library functions.
 * Sets up internal tables and default configuration.
 * 
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_init(void);

/**
 * @brief Clean up library resources
 * 
 * Should be called when done using the library.
 */
void swl_cleanup(void);

/**
 * @brief Get default configuration
 * 
 * @param config Pointer to configuration structure to fill
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_get_default_config(SWL_Config* config);

/**
 * @brief Set library configuration
 * 
 * @param config Pointer to configuration structure
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_set_config(const SWL_Config* config);

/* ============================================================================
 * MAIN API FUNCTIONS
 * ============================================================================ */

/**
 * @brief Calculate wind yields for a star cluster at a given age
 * 
 * This is the main function for calculating time-dependent stellar wind
 * yields from a star cluster. It integrates over the IMF and stellar
 * lifetimes to compute the cumulative wind yields up to the given age.
 * 
 * @param imf_type      IMF model to use (Salpeter, Kroupa, Chabrier, etc.)
 * @param metallicity   Initial metallicity Z of the cluster
 * @param cluster_mass  Total initial mass of the cluster (M_sun)
 * @param cluster_age   Age of the cluster (years)
 * @param result        Pointer to result structure to fill
 * 
 * @return SWL_SUCCESS on success, error code otherwise
 * 
 * @note The returned yields represent the NET yields (ejected - initial).
 *       Negative values mean the element was consumed (e.g., H, He).
 * 
 * Example usage:
 * @code
 *     SWL_ClusterWindYield result;
 *     swl_init();
 *     
 *     // Calculate yields for a 10^6 M_sun cluster at 10 Myr
 *     swl_calculate_cluster_wind_yield(
 *         SWL_IMF_CHABRIER,    // IMF model
 *         0.02,                // Solar metallicity
 *         1.0e6,               // 10^6 M_sun cluster
 *         1.0e7,               // 10 Myr age
 *         &result
 *     );
 *     
 *     printf("Total wind mass: %e M_sun\n", result.total_wind_mass);
 *     printf("Carbon yield: %e M_sun\n", result.element_yields.yields[SWL_EL_C]);
 *     
 *     swl_cleanup();
 * @endcode
 */
SWL_ErrorCode swl_calculate_cluster_wind_yield(
    SWL_IMFType imf_type,
    double metallicity,
    double cluster_mass,
    double cluster_age,
    SWL_ClusterWindYield* result
);

/**
 * @brief Calculate wind yields for a single star
 * 
 * Returns the total wind yields integrated over the star's lifetime
 * (or up to the current stellar age if star is still alive).
 * 
 * @param initial_mass  Initial stellar mass (M_sun)
 * @param metallicity   Initial metallicity Z
 * @param stellar_age   Age of the star (years), or -1 for lifetime-integrated
 * @param result        Pointer to result structure to fill
 * 
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_calculate_star_wind_yield(
    double initial_mass,
    double metallicity,
    double stellar_age,
    SWL_StarWindYield* result
);

/**
 * @brief Calculate IMF-integrated yields per unit stellar mass formed
 * 
 * This function returns the yields normalized to 1 M_sun of stars formed,
 * useful for chemical evolution models.
 * 
 * @param imf_type      IMF model to use
 * @param metallicity   Initial metallicity Z
 * @param yields        Pointer to yield structure to fill
 * @param wind_mass     Pointer to store total wind mass per M_sun formed
 * 
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_calculate_imf_integrated_yields(
    SWL_IMFType imf_type,
    double metallicity,
    SWL_ElementYields* yields,
    double* wind_mass
);

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================ */

/**
 * @brief Get stellar lifetime for a given mass and metallicity
 * 
 * @param mass          Stellar mass (M_sun)
 * @param metallicity   Metallicity Z
 * 
 * @return Stellar lifetime in years
 */
double swl_stellar_lifetime(double mass, double metallicity);

/**
 * @brief Get turnoff mass for a given age and metallicity
 * 
 * Returns the mass of stars currently at the main sequence turnoff.
 * 
 * @param age           Age (years)
 * @param metallicity   Metallicity Z
 * 
 * @return Turnoff mass in M_sun
 */
double swl_turnoff_mass(double age, double metallicity);

/**
 * @brief Get mass-loss rate for a star
 * 
 * @param mass          Stellar mass (M_sun)
 * @param luminosity    Stellar luminosity (L_sun)
 * @param T_eff         Effective temperature (K)
 * @param metallicity   Metallicity Z
 * @param phase         Stellar evolutionary phase
 * 
 * @return Mass-loss rate in M_sun/yr
 */
double swl_mass_loss_rate(
    double mass, 
    double luminosity, 
    double T_eff,
    double metallicity,
    SWL_StellarPhase phase
);

/**
 * @brief Get IMF value at a given mass
 * 
 * Returns the normalized IMF phi(m) such that:
 * integral of m * phi(m) dm from m_min to m_max = 1 M_sun
 * 
 * @param mass          Stellar mass (M_sun)
 * @param imf_type      IMF model
 * 
 * @return IMF value (normalized)
 */
double swl_imf(double mass, SWL_IMFType imf_type);

/**
 * @brief Get remnant mass for a given initial mass
 * 
 * @param initial_mass  Initial stellar mass (M_sun)
 * @param metallicity   Metallicity Z
 * 
 * @return Remnant mass in M_sun (WD, NS, or BH)
 */
double swl_remnant_mass(double initial_mass, double metallicity);

/**
 * @brief Get element name string
 * 
 * @param element       Element index
 * 
 * @return Pointer to element name string
 */
const char* swl_get_element_name(SWL_Element element);

/**
 * @brief Get IMF name string
 * 
 * @param imf_type      IMF type
 * 
 * @return Pointer to IMF name string
 */
const char* swl_get_imf_name(SWL_IMFType imf_type);

/**
 * @brief Get solar abundance for an element
 * 
 * @param element       Element index
 * 
 * @return Solar mass fraction (Asplund et al. 2009)
 */
double swl_get_solar_abundance(SWL_Element element);

/**
 * @brief Get error message for error code
 * 
 * @param error         Error code
 * 
 * @return Pointer to error message string
 */
const char* swl_get_error_message(SWL_ErrorCode error);

/* ============================================================================
 * ADVANCED FUNCTIONS
 * ============================================================================ */

/**
 * @brief Calculate instantaneous wind yield rate
 * 
 * Returns the current rate of element ejection from stellar winds
 * at the given cluster age.
 * 
 * @param imf_type      IMF model
 * @param metallicity   Metallicity Z
 * @param cluster_mass  Cluster mass (M_sun)
 * @param cluster_age   Cluster age (years)
 * @param yield_rate    Pointer to yield rate structure (M_sun/yr)
 * @param mass_rate     Pointer to store total mass-loss rate (M_sun/yr)
 * 
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_calculate_wind_rate(
    SWL_IMFType imf_type,
    double metallicity,
    double cluster_mass,
    double cluster_age,
    SWL_ElementYields* yield_rate,
    double* mass_rate
);

/**
 * @brief Calculate cumulative yields at multiple ages
 * 
 * Efficiently computes yields at an array of ages.
 * 
 * @param imf_type      IMF model
 * @param metallicity   Metallicity Z
 * @param cluster_mass  Cluster mass (M_sun)
 * @param ages          Array of ages (years)
 * @param n_ages        Number of ages
 * @param results       Array of result structures to fill
 * 
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_calculate_yield_history(
    SWL_IMFType imf_type,
    double metallicity,
    double cluster_mass,
    const double* ages,
    int n_ages,
    SWL_ClusterWindYield* results
);

/**
 * @brief Export yields to CSV file
 * 
 * @param result        Pointer to result structure
 * @param filename      Output filename
 * 
 * @return SWL_SUCCESS on success, error code otherwise
 */
SWL_ErrorCode swl_export_to_csv(
    const SWL_ClusterWindYield* result,
    const char* filename
);

#ifdef __cplusplus
}
#endif

#endif /* STELLAR_WIND_LIB_H */

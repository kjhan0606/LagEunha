/*
 * Stellar Cluster Supernova Library - Implementation
 */

#include "cluster_supernova.h"
#include <string.h>

/* ===== Internal Constants ===== */

#define N_METALLICITIES 5
#define SOLAR_Z 0.02

/* Typical supernova energies (in foe = 10^51 erg) */
#define ENERGY_SN_IA 1.0    // Type Ia: ~1 foe
#define ENERGY_SN_II 1.0    // Type II: ~1 foe (range 0.5-2)

/* Type Ia parameters */
//#define DTD_NORMALIZATION 1.3e-3  // Normalization for DTD (per M_sun per year)
#define DTD_NORMALIZATION 8.e-3  // Normalization for DTD (per M_sun per year)
#define MIN_DELAY_SNIA 40e6       // Minimum delay time (40 Myr)
#define MAX_DELAY_SNIA 10e9       // Maximum delay time (10 Gyr)

/* ===== Metallicity Grid and Yield Tables ===== */

typedef struct {
    double Z;
    double log_Z_solar;
} MetallicityPoint;

static const MetallicityPoint Z_grid[N_METALLICITIES] = {
    {0.0,      -999.0},
    {0.0001,   -2.30},
    {0.001,    -1.30},
    {0.02,      0.0},
    {0.10,      0.70}
};

/* ===== Nucleosynthesis Yield Tables ===== */
/*
 * REFERENCES:
 * 
 * Type Ia Supernovae:
 *   Iwamoto et al. (1999), ApJS, 125, 439
 *   Model: W7 (deflagration, 1.4 M_sun WD)
 * 
 * Type II Supernovae:
 *   Limongi & Chieffi (2018), ApJS, 237, 13
 *   Models: Non-rotating, various masses and metallicities
 *   Note: Ca yields calibrated to Cayrel et al. (2004), A&A, 416, 1117
 * 
 * IMF:
 *   Kroupa (2001), MNRAS, 322, 231
 * 
 * For detailed references, see YIELD_TABLE_REFERENCES.md
 */

/* Element yield tables (same as in previous code) */
/* * Corrected Yield Table based on Physics of Mass Loss
 * * Correction 1 (H/He Scaling):
 * - High mass stars (25M) lose more envelope -> Less H/He left (Factor < 1.0).
 * - Low mass stars (15M) lose less envelope -> relative to core, comparable H (Factor ~ 0.8-1.0 due to smaller total mass).
 * * Correction 2 (Mg):
 * - Increased Mg yield slightly to match standard O/Mg ratios (~0.15-0.20).
 *
 * Correction 3 (25M Metal Scaling):
 * - 25M stars produce significantly more O/Mg/C than 20M stars. Factor increased to ~1.8.
 */
/* Element yield tables (same as in previous code) */
typedef struct {
    const char *name;
    double atomic_mass;
    double sn_ia_yield[N_METALLICITIES];
    double sn_ii_yield_20M[N_METALLICITIES];
    double mass_scaling_15M;
    double mass_scaling_25M;
} YieldTable;

static const YieldTable yield_table[N_ELEMENTS] = {
    /* name  mass    Type Ia yields (Z grid)             Type II yields (20 M_sun)         15M    25M */

    // H & He: 25M loses most of its envelope. 15M keeps more relative to 20M, but total mass is less.
    {"H",  1.008, {0.0,    0.0,    0.0,    0.0,    0.0},     {6.5,   6.0,   5.5,   4.8,   3.8},      0.90, 0.60},
    {"He", 4.003, {0.0,    0.0,    0.0,    0.0,    0.0},     {3.8,   3.5,   3.0,   2.5,   1.9},      0.85, 0.70},

    // Carbon: 25M produces much more.
    {"C",  12.01, {0.040,  0.042,  0.045,  0.048,  0.052},   {0.25,  0.22,  0.20,  0.18,  0.16},     0.60, 1.80},

    // Nitrogen: Secondary element, stronger in higher mass/Z.
    {"N",  14.01, {1.0e-5, 1.1e-5, 1.2e-5, 1.2e-5, 1.3e-5},  {0.005, 0.008, 0.010, 0.012, 0.020},    0.60, 1.60},

    // Alpha Elements (O, Mg, Si, S, Ca): Strong scaling with mass (25M >> 20M).
    {"O",  16.00, {0.130,  0.135,  0.138,  0.140,  0.145},   {2.8,   2.5,   2.1,   1.8,   1.5},      0.55, 2.20},
    {"Ne", 20.18, {0.0042, 0.0044, 0.0046, 0.0048, 0.0051},  {0.38,  0.35,  0.32,  0.28,  0.24},     0.55, 2.00},
    // Mg corrected from ~0.11 to ~0.18
    {"Mg", 24.31, {0.0075, 0.0080, 0.0083, 0.0085, 0.0090},  {0.20,  0.18,  0.16,  0.14,  0.12},     0.60, 2.10},
    {"Si", 28.09, {0.150,  0.155,  0.158,  0.160,  0.165},   {0.30,  0.28,  0.26,  0.24,  0.21},     0.65, 1.80},
    {"S",  32.07, {0.080,  0.083,  0.084,  0.085,  0.088},   {0.14,  0.13,  0.12,  0.11,  0.095},    0.70, 1.70},
    {"Ca", 40.08, {0.012,  0.014,  0.015,  0.016,  0.018},   {0.015, 0.013, 0.011, 0.009, 0.008},    0.70, 1.60},

    // Iron: Weak dependence on mass in Core Collapse (determined by mass cut/explosion energy).
    {"Fe", 55.85, {0.710,  0.720,  0.730,  0.740,  0.755},   {0.10,  0.09,  0.088, 0.085, 0.078},    0.80, 1.10}
};

/* ===== Pair-Instability Supernova (PISN) Yield Tables ===== */
/*
 * REFERENCES:
 * 
 * Pair-Instability Supernovae (Pop III):
 *   Heger & Woosley (2002), ApJ, 567, 532
 *   Heger et al. (2003), ApJ, 591, 288
 *   
 * Mass range: 140-260 M_sun
 * Energy: 10-100 foe (10^52 - 10^53 erg)
 * 
 * Key features:
 *   - Complete disruption (no remnant)
 *   - Odd-even effect (low Na, Al, P, etc.)
 *   - Massive 56Ni production (up to 57 M_sun)
 *   - No elements > Zn (no r/s-process)
 */

#define ENERGY_PISN_BASE 10.0  // Base PISN energy in foe (10^52 erg)

typedef struct {
    double mass;           // Progenitor mass (M_sun)
    double he_core;        // Helium core mass (M_sun)
    double ni56;           // 56Ni mass produced (M_sun)
    double energy_foe;     // Explosion energy (foe)
    double yields[N_ELEMENTS];  // Element yields (M_sun) - THREAD-SAFE: const data
} PISNModel;

/* * PISN models based on Heger & Woosley (2002), ApJ, 567, 532
 * Re-calibrated for Mass Conservation (Sum of Yields approx He_core)
 * * Note 1: 'H' is 0.0 assuming we are looking at the Helium Core yields.
 * Note 2: 'Fe(final)' implies decayed 56Ni + other Fe isotopes.
 * Note 3: Pop III stars produce negligible Nitrogen.
 */

static const PISNModel pisn_models[N_PISN_MODELS] = {
    /* M_init   He_core   56Ni    Energy    H      He      C       N        O       Ne      Mg      Si      S       Ca      Fe(final) */

    // Model 1: 150 Msun (Weak explosion, lots of unburned He, low Ni)
    // He_core ~ 72.0 -> Yield Sum ~ 72.0
    {150,       72.0,     0.10,   15.0,    {0.0,   16.0,   0.60,   1.0e-5,  42.0,   3.5,    2.0,    5.0,    1.8,    0.4,    0.10}},

    // Model 2: 175 Msun (Intermediate)
    // He_core ~ 86.0 -> Yield Sum ~ 86.0
    {175,       86.0,     2.50,   25.0,    {0.0,   10.0,   0.90,   1.0e-5,  46.0,   4.5,    3.0,    14.0,   4.0,    1.0,    2.50}},

    // Model 3: 200 Msun (Standard PISN)
    // He_core ~ 100.0 -> Yield Sum ~ 100.0
    {200,       100.0,    16.0,   45.0,    {0.0,   6.0,    1.20,   1.0e-5,  44.0,   5.0,    3.5,    18.0,   5.0,    1.3,    16.0}},

    // Model 4: 225 Msun (High Energy)
    // He_core ~ 114.0 -> Yield Sum ~ 114.0
    {225,       114.0,    30.0,   65.0,    {0.0,   3.5,    1.80,   1.0e-5,  40.0,   6.0,    3.8,    22.0,   6.0,    1.5,    30.0}},

    // Model 5: 250 Msun (Extreme PISN, almost full burning)
    // He_core ~ 128.0 -> Yield Sum ~ 128.0
    // Very little He left. Most mass is O, Si, Fe(Ni).
    {250,       128.0,    50.0,   85.0,    {0.0,   2.0,    2.50,   1.0e-5,  35.0,   7.0,    4.0,    24.0,   8.0,    2.0,    50.0}}
};
/* ===== Internal Helper Functions ===== */

static double interpolate_metallicity(double Z, const double values[N_METALLICITIES]) {
    if (Z <= Z_grid[0].Z) return values[0];
    if (Z >= Z_grid[N_METALLICITIES-1].Z) return values[N_METALLICITIES-1];
    if (Z < 1e-6) return values[0];
    
    int i_lower = 0;
    for (int i = 0; i < N_METALLICITIES - 1; i++) {
        if (Z >= Z_grid[i].Z && Z < Z_grid[i+1].Z) {
            i_lower = i;
            break;
        }
    }
    
    double log_Z = log10(Z);
    double log_Z_lower = log10(Z_grid[i_lower].Z + 1e-10);
    double log_Z_upper = log10(Z_grid[i_lower+1].Z);
    double frac = (log_Z - log_Z_lower) / (log_Z_upper - log_Z_lower);
    
    return values[i_lower] * (1.0 - frac) + values[i_lower+1] * frac;
}

/* ===== IMF Functions ===== */

double evaluate_imf(double mass, IMFType imf_type) {
    /* Allow evaluation beyond MAX_MASS_STARS for Pop III/PISN calculations */
    if (mass < MIN_MASS_STARS) return 0.0;
    
    switch (imf_type) {
        case IMF_SALPETER:
            /* Salpeter (1955): dN/dM ∝ M^(-2.35) */
            return pow(mass, -2.35);
            
        case IMF_KROUPA:
            /* Kroupa (2001): Broken power law */
            if (mass < 0.5) {
                return pow(mass / 0.5, -1.3) * pow(0.5, -2.3);
            } else {
                return pow(mass, -2.3);
            }
            
        case IMF_CHABRIER:
            /* Chabrier (2003): Log-normal for M < 1, power law for M > 1 */
            if (mass < 1.0) {
                double mc = 0.079;  // Characteristic mass
                double sigma = 0.69;
                double log_m = log10(mass);
                double log_mc = log10(mc);
                return (1.0 / mass) * exp(-pow(log_m - log_mc, 2) / (2.0 * sigma * sigma));
            } else {
                return pow(mass, -2.3);
            }
            
        case IMF_POPIII_TOPHEAVY:
            /* 
             * Pop III top-heavy IMF (for Z ~ 0)
             * Based on simulations: Hirano et al. (2014), Susa et al. (2014)
             * 
             * Characteristics:
             * - Peak around 100 M_sun (very different from normal ~0.3 M_sun)
             * - No low-mass stars (M < 10 M_sun)
             * - Flatter slope (M^-1.0 instead of M^-2.3)
             * - Extends to very high masses (up to 300 M_sun)
             * 
             * Form:
             * M < 10 M_sun: ξ = 0 (no low-mass Pop III)
             * 10 < M < 100 M_sun: ξ ∝ (M/100)^0.5 (rising to peak)
             * M > 100 M_sun: ξ ∝ (M/100)^-1.0 (top-heavy tail)
             */
            if (mass < 10.0) {
                return 0.0;  // No low-mass Pop III stars
            } else if (mass < 100.0) {
                /* Rising to peak at 100 M_sun */
                return pow(mass / 100.0, 0.5);
            } else {
                /* Top-heavy tail: much flatter than Salpeter */
                return pow(mass / 100.0, -1.0);
            }
            
        default:
            return pow(mass, -2.35);  // Default to Salpeter
    }
}

double calculate_imf_normalization(IMFType imf_type) {
    /* Calculate normalization so that ∫ M * IMF(M) dM = 1 M_sun */
    double sum = 0.0;
    int n_bins = 1000;
    double log_min = log10(MIN_MASS_STARS);
    double log_max = log10(MAX_MASS_STARS);
    double dlog_m = (log_max - log_min) / n_bins;
    
    for (int i = 0; i < n_bins; i++) {
        double log_m = log_min + (i + 0.5) * dlog_m;
        double m = pow(10.0, log_m);
        double dm = m * log(10.0) * dlog_m;
        sum += m * evaluate_imf(m, imf_type) * dm;
    }
    
    return 1.0 / sum;
}

/* Calculate fraction of stellar mass that becomes Type II SNe */
double calculate_snii_mass_fraction(IMFType imf_type) {
    /*
     * For continuous SF: NUMBER of SN II per M_sun of stars formed
     * NOT mass fraction!
     * 
     * Rate = ∫ ξ(M) dM over 8-40 M_sun
     * where ξ is normalized IMF
     */
    double imf_norm = calculate_imf_normalization(imf_type);
    
    double n_snii = 0.0;  // Number, not mass!
    int n_bins = 500;
    double log_min = log10(MIN_MASS_SNII);
    double log_max = log10(MAX_MASS_SNII);
    double dlog_m = (log_max - log_min) / n_bins;
    
    for (int i = 0; i < n_bins; i++) {
        double log_m = log_min + (i + 0.5) * dlog_m;
        double m = pow(10.0, log_m);
        double dm = m * log(10.0) * dlog_m;
        /* Number of stars: ξ(M) dM, NOT M*ξ(M) dM */
        n_snii += imf_norm * evaluate_imf(m, imf_type) * dm;
    }
    
    return n_snii;
}

/* ===== Stellar Lifetime Function ===== */

double stellar_lifetime(double mass, double metallicity) {
    /*
     * Stellar main sequence lifetime
     * Based on Hurley, Pols, & Tout (2000)
     * 
     * For M > 0.5 M_sun: t ~ M^(-2.5) * f(Z)
     * Adjusted for metallicity effects
     */
    
    if (mass < MIN_MASS_STARS) return 1e12;  // Very long for brown dwarfs
    
    /* Base lifetime (solar metallicity) */
    double t_ms;
    if (mass < 0.5) {
        t_ms = 1.0e10 * pow(mass, -2.0);  // Low mass stars
    } else {
        t_ms = 1.0e10 * pow(mass, -2.5);  // Intermediate and high mass
    }
    
    /* Metallicity correction */
    double Z_ratio = metallicity / SOLAR_Z;
    if (Z_ratio < 0.01) Z_ratio = 0.01;  // Avoid extreme values
    
    /* Higher metallicity -> stronger winds -> shorter lifetime */
    double Z_factor = 1.0 - 0.1 * log10(Z_ratio);
    if (Z_factor < 0.5) Z_factor = 0.5;
    
    return t_ms * Z_factor;
}

/* ===== PISN Functions (THREAD-SAFE) ===== */

/*
 * Get PISN lifetime
 * Very massive Pop III stars have short lifetimes
 * Based on Heger & Woosley (2002)
 * 
 * THREAD-SAFE: Pure function, no global state
 */
static double pisn_lifetime(double mass) {
    /* 
     * For very massive Pop III stars (140-260 M_sun):
     * Approximate lifetime: t ~ 2.5-3.5 Myr
     * 
     * Scaling: t ∝ M^-1.5 (approximately)
     */
    if (mass < MIN_MASS_PISN) return 1e10;  // Not a PISN
    if (mass > MAX_MASS_PISN) return 1e10;  // Direct BH collapse
    
    /* Reference: 200 M_sun -> 2.8 Myr */
    double t_ref = 2.8e6;  // years
    double m_ref = 200.0;  // M_sun
    
    return t_ref * pow(m_ref / mass, 1.5);
}

/*
 * Get PISN yields by interpolation
 * THREAD-SAFE: Only reads static const data
 * 
 * @param mass: Progenitor mass (M_sun)
 * @param element: Element index
 * @return: Yield in M_sun
 */
static double get_pisn_yield(double mass, int element) {
    /* Outside PISN range */
    if (mass < MIN_MASS_PISN || mass > MAX_MASS_PISN) return 0.0;
    if (element < 0 || element >= N_ELEMENTS) return 0.0;
    
    /* Find bracketing models - THREAD-SAFE: read-only */
    int i_lower = 0;
    for (int i = 0; i < N_PISN_MODELS - 1; i++) {
        if (mass >= pisn_models[i].mass && mass < pisn_models[i+1].mass) {
            i_lower = i;
            break;
        }
    }
    
    /* Handle boundaries */
    if (mass <= pisn_models[0].mass) {
        return pisn_models[0].yields[element];
    }
    if (mass >= pisn_models[N_PISN_MODELS-1].mass) {
        return pisn_models[N_PISN_MODELS-1].yields[element];
    }
    
    /* Linear interpolation - THREAD-SAFE: all local variables */
    double m1 = pisn_models[i_lower].mass;
    double m2 = pisn_models[i_lower + 1].mass;
    double y1 = pisn_models[i_lower].yields[element];
    double y2 = pisn_models[i_lower + 1].yields[element];
    
    double frac = (mass - m1) / (m2 - m1);
    return y1 * (1.0 - frac) + y2 * frac;
}

/*
 * Get PISN explosion energy
 * THREAD-SAFE: Pure function
 * 
 * @param mass: Progenitor mass (M_sun)
 * @return: Energy in erg
 */
static double get_pisn_energy(double mass) {
    if (mass < MIN_MASS_PISN || mass > MAX_MASS_PISN) return 0.0;
    
    /* Find bracketing models */
    int i_lower = 0;
    for (int i = 0; i < N_PISN_MODELS - 1; i++) {
        if (mass >= pisn_models[i].mass && mass < pisn_models[i+1].mass) {
            i_lower = i;
            break;
        }
    }
    
    /* Boundaries */
    if (mass <= pisn_models[0].mass) {
        return pisn_models[0].energy_foe * FOE;
    }
    if (mass >= pisn_models[N_PISN_MODELS-1].mass) {
        return pisn_models[N_PISN_MODELS-1].energy_foe * FOE;
    }
    
    /* Interpolate energy */
    double m1 = pisn_models[i_lower].mass;
    double m2 = pisn_models[i_lower + 1].mass;
    double e1 = pisn_models[i_lower].energy_foe;
    double e2 = pisn_models[i_lower + 1].energy_foe;
    
    double frac = (mass - m1) / (m2 - m1);
    double energy_foe = e1 * (1.0 - frac) + e2 * frac;
    
    return energy_foe * FOE;  // Convert to erg
}

/*
 * Calculate PISN contribution to cluster
 * THREAD-SAFE: No global state, all calculations local
 * 
 * @param cluster: Cluster properties (const - read-only)
 * @param pisn_rate: Output - PISN rate (/year)
 * @param pisn_energy: Output - Energy release (erg/year)
 * @param pisn_mass: Output - Total ejected mass (M_sun/year)
 * @param pisn_yields: Output - Element yields (M_sun/year)
 */
static void calculate_pisn_contribution(
    const ClusterProperties *cluster,
    double *pisn_rate,
    double *pisn_energy,
    double *pisn_mass,
    double pisn_yields[N_ELEMENTS]
) {
    /* Initialize outputs - THREAD-SAFE: local pointers */
    *pisn_rate = 0.0;
    *pisn_energy = 0.0;
    *pisn_mass = 0.0;
    for (int i = 0; i < N_ELEMENTS; i++) {
        pisn_yields[i] = 0.0;
    }
    
    /* Check if Pop III is enabled */
    if (cluster->popiii_mode != POPIII_ENABLED) return;
    
    /* PISN only for very low metallicity */
    if (cluster->metallicity > 0.001) return;  // Z > 0.001 -> no PISN
    
    /* Determine mass range */
    double m_min = MIN_MASS_PISN;  // 140 M_sun
    double m_max = cluster->popiii_imf_max;
    if (m_max > MAX_MASS_PISN) m_max = MAX_MASS_PISN;  // 260 M_sun
    if (m_max <= m_min) return;
    
    /*
     * SIMPLIFIED PISN CALCULATION
     * 
     * For a cluster of total_mass M_cluster formed with IMF ξ(M):
     * 
     * Number of stars in mass range [M, M+dM]:
     *   dN = M_cluster * ξ(M) dM / M
     * 
     * where ξ is normalized: ∫ M * ξ(M) dM = 1 over standard range
     * 
     * For PISN range (140-260 M_sun), we use same ξ but account for
     * the fact that this is ADDITIONAL mass beyond the standard 0.08-120 range.
     */
    
    /* Get IMF normalization (for 0.08-120 M_sun range) */
    double imf_norm = calculate_imf_normalization(cluster->imf_type);
    
    /* Pop III mass fraction from cluster properties (default: 10%) */
    double popiii_mass_fraction = cluster->popiii_mass_fraction;
    
    /* Validate and apply default if needed */
    if (popiii_mass_fraction <= 0.0 || popiii_mass_fraction > 1.0) {
        popiii_mass_fraction = 0.10;  // Default to 10%
    }
    
    double mass_in_popiii = cluster->total_mass * popiii_mass_fraction;
    
    #ifdef DEBUG_PISN
    fprintf(stderr, "DEBUG PISN: IMF norm = %.6e\n", imf_norm);
    fprintf(stderr, "DEBUG PISN: Pop III mass = %.6e M_sun (%.1f%% of total)\n",
            mass_in_popiii, popiii_mass_fraction * 100.0);
    fprintf(stderr, "DEBUG PISN: Age = %.3f Myr\n", cluster->age / 1e6);
    #endif
    
    /* Integration over PISN mass range */
    const int n_bins = 60;
    double log_min = log10(m_min);
    double log_max = log10(m_max);
    double dlog_m = (log_max - log_min) / n_bins;
    
    int n_bins_matched = 0;  // Debug counter
    
    /* THREAD-SAFE: All loop variables are local */
    for (int i = 0; i < n_bins; i++) {
        double log_m = log_min + (i + 0.5) * dlog_m;
        double m = pow(10.0, log_m);
        double dm = m * log(10.0) * dlog_m;
        
        /* IMF value at this mass */
        double xi_m = evaluate_imf(m, cluster->imf_type);
        
        /* 
         * Number of stars in this bin:
         * 
         * IMF normalization: ∫ M * ξ(M) dM = 1
         * So: imf_norm * ∫ M * IMF(M) dM = 1
         * 
         * Number of stars in bin [M, M+dM]:
         * dN = (available_mass) * imf_norm * ξ(M) dM
         * 
         * NOTE: No division by M! That's already in the normalization.
         */
        double n_stars = mass_in_popiii * imf_norm * xi_m * dm;
        
        #ifdef DEBUG_PISN
        if (i < 3 || (i >= 28 && i <= 35)) {
            fprintf(stderr, "Bin %d: M=%.1f, xi=%.3e, imf*xi*dm=%.3e, n_stars=%.3e\n",
                    i, m, xi_m, imf_norm * xi_m * dm, n_stars);
        }
        #endif
        
        /* PISN lifetime for this mass */
        double lifetime = pisn_lifetime(m);
        
        /* Calculate rate contribution */
        double rate_contribution = 0.0;
        
        if (cluster->sf_mode == SF_INSTANTANEOUS) {
            /*
             * Instantaneous burst at t=0
             * Stars of mass M explode at t = lifetime(M)
             * 
             * PISN: lifetime ~ 2-3.5 Myr for 260-140 M_sun
             * All PISN explode in a narrow time window
             */
            
            /* Time window for "exploding now" */
            double time_window = 0.3e6;  // 0.3 Myr
            
            /* Is current age within window of this mass's lifetime? */
            if (fabs(cluster->age - lifetime) < time_window) {
                /* Yes - these stars are exploding now */
                rate_contribution = n_stars / (2.0 * time_window);
                n_bins_matched++;
                
                #ifdef DEBUG_PISN
                fprintf(stderr, "DEBUG PISN: Bin %d matched!\n", i);
                fprintf(stderr, "  Mass = %.1f M_sun, lifetime = %.3f Myr\n",
                        m, lifetime / 1e6);
                fprintf(stderr, "  n_stars = %.6e, rate = %.6e /yr\n",
                        n_stars, rate_contribution);
                #endif
            }
        } else {
            /* Continuous star formation */
            if (cluster->age >= lifetime) {
                /* Stars formed continuously */
                double sfr_popiii = cluster->star_formation_rate * popiii_mass_fraction;
                double n_stars_per_year = sfr_popiii * imf_norm * xi_m * dm;
                rate_contribution = n_stars_per_year;
            }
        }
        
        /* Accumulate outputs - THREAD-SAFE: local variables */
        *pisn_rate += rate_contribution;
        
        /* Energy from these PISN */
        double energy_per_sn = get_pisn_energy(m);
        *pisn_energy += rate_contribution * energy_per_sn;
        
        /* Elemental yields */
        double total_ejecta = 0.0;
        for (int elem = 0; elem < N_ELEMENTS; elem++) {
            double yield_elem = get_pisn_yield(m, elem);
            pisn_yields[elem] += rate_contribution * yield_elem;
            total_ejecta += yield_elem;
        }
        
        *pisn_mass += rate_contribution * total_ejecta;
    }
    
    #ifdef DEBUG_PISN
    fprintf(stderr, "DEBUG PISN: Matched %d / %d bins\n", n_bins_matched, n_bins);
    fprintf(stderr, "DEBUG PISN: Total PISN rate = %.6e /yr\n", *pisn_rate);
    fprintf(stderr, "DEBUG PISN: Total PISN energy = %.6e erg/yr\n", *pisn_energy);
    fprintf(stderr, "\n");
    #endif
}

/* ===== Type Ia DTD ===== */

double type_ia_dtd(double time) {
    /*
     * Delay Time Distribution for Type Ia supernovae
     * Based on Maoz & Graur (2017): dN/dt ∝ t^(-1)
     * 
     * Returns: probability per unit time per unit stellar mass
     */
    
    if (time < MIN_DELAY_SNIA || time > MAX_DELAY_SNIA) {
        return 0.0;
    }
    
    /* Power law: t^(-1) */
    double dtd = DTD_NORMALIZATION / time;
    
    return dtd;
}

/* ===== Main Calculation Functions ===== */

int validate_cluster_properties(const ClusterProperties *cluster) {
    if (cluster == NULL) return -1;
    if (cluster->total_mass <= 0.0) return -1;
    if (cluster->age < 0.0) return -1;
    if (cluster->metallicity < 0.0 || cluster->metallicity > 0.15) return -1;
    return 0;
}

int calculate_sn_rates(const ClusterProperties *cluster, SNRates *rates) {
    if (validate_cluster_properties(cluster) != 0) return -1;
    if (rates == NULL) return -1;
    
    double age = cluster->age;
    double total_mass = cluster->total_mass;
    double Z = cluster->metallicity;
    IMFType imf = cluster->imf_type;
    StarFormationMode sf_mode = cluster->sf_mode;
    
    /* Calculate IMF normalization */
    double imf_norm = calculate_imf_normalization(imf);
    
    /* ===== CONTINUOUS STAR FORMATION MODE ===== */
    if (sf_mode == SF_CONTINUOUS) {
        double sfr = cluster->star_formation_rate;
        
        if (sfr <= 0.0) {
            fprintf(stderr, "Warning: SFR must be positive for SF_CONTINUOUS mode\n");
            rates->sn_ia_rate = 0.0;
            rates->sn_ii_rate = 0.0;
            rates->total_rate = 0.0;
            return 0;
        }
        
        /* Type II rate: Recent star formation */
        /* SN II progenitors (8-40 M_sun) have lifetimes 3-40 Myr */
        /* For constant SFR, rate = SFR × (mass fraction in 8-40 M_sun) */
        double f_snii = calculate_snii_mass_fraction(imf);
        
        /* Time delay: average lifetime of SN II progenitors ~20 Myr */
        /* For ages < 40 Myr, need to account for startup */
        double max_snii_lifetime = stellar_lifetime(MIN_MASS_SNII, Z);
        
        if (age < max_snii_lifetime) {
            /* Partial contribution: only stars old enough to explode */
            double active_fraction = age / max_snii_lifetime;
            rates->sn_ii_rate = sfr * f_snii * active_fraction;
        } else {
            /* Steady state: full SN II rate */
            rates->sn_ii_rate = sfr * f_snii;
        }
        
        /* Type Ia rate: Integral of past star formation with DTD */
        /* For constant SFR: ∫ SFR × DTD(t-t') dt' from 0 to age */
        rates->sn_ia_rate = 0.0;
        
        if (age > MIN_DELAY_SNIA) {
            /* Integrate DTD */
            int n_time_bins = 100;
            double t_min = MIN_DELAY_SNIA;
            double t_max = (age < MAX_DELAY_SNIA) ? age : MAX_DELAY_SNIA;
            double dt = (t_max - t_min) / n_time_bins;
            
            for (int i = 0; i < n_time_bins; i++) {
                double t = t_min + (i + 0.5) * dt;
                rates->sn_ia_rate += sfr * type_ia_dtd(t) * dt;
            }
        }
        
        rates->total_rate = rates->sn_ia_rate + rates->sn_ii_rate;
        return 0;
    }
    
    /* ===== INSTANTANEOUS BURST MODE ===== */
    
    /* Type II supernova rate */
    /* 
     * For instantaneous burst: Stars in mass bin dm die over time dτ
     * Rate = (number of stars in dm) / (time to exhaust dm)
     * 
     * Use window method: count stars dying within Δt of current age
     * and normalize by window size
     */
    double sn_ii_rate = 0.0;
    
    /* Use adaptive window based on cluster age */
    /* Younger clusters: narrower window (faster evolution) */
    /* Older clusters: wider window (slower evolution) */
    double window_years;
    if (age < 30e6) {
        window_years = 5e6;  // 5 Myr for young clusters
    } else if (age < 100e6) {
        window_years = 20e6;  // 20 Myr for intermediate age
    } else {
        window_years = 50e6;  // 50 Myr for old clusters
    }
    
    /* Count stars dying within window */
    double n_stars_in_window = 0.0;
    
    int n_mass_bins = 200;
    double log_m_min = log10(MIN_MASS_SNII);
    double log_m_max = log10(MAX_MASS_SNII);
    double dlog_m = (log_m_max - log_m_min) / n_mass_bins;
    
    for (int i = 0; i < n_mass_bins; i++) {
        double log_m = log_m_min + (i + 0.5) * dlog_m;
        double m = pow(10.0, log_m);
        double dm = m * log(10.0) * dlog_m;
        
        /* Get lifetime for this mass */
        double lifetime = stellar_lifetime(m, Z);
        
        /* Check if this mass dies within window */
        if (fabs(age - lifetime) < window_years / 2.0) {
            /* Number of stars in this mass bin */
            double n_stars = total_mass * imf_norm * evaluate_imf(m, imf) * dm;
            n_stars_in_window += n_stars;
        }
    }
    
    /* Calculate rate: stars per window / years per window */
    sn_ii_rate = n_stars_in_window / window_years;
    
    /* Type Ia supernova rate */
    /* Use delay time distribution */
    double sn_ia_rate = 0.0;
    
    if (age > MIN_DELAY_SNIA) {
        /* Binary fraction ~ 0.5 for solar-type stars */
        double binary_fraction = 0.5;
        
        /* Mass range for Type Ia progenitors (A stars, ~1-3 M_sun) */
        double m_min_ia = 1.0;
        double m_max_ia = 8.0;
        
        /* Calculate total stellar mass in this range */
        double stellar_mass_ia = 0.0;
        int n_ia_bins = 50;
        double log_m_ia_min = log10(m_min_ia);
        double log_m_ia_max = log10(m_max_ia);
        double dlog_m_ia = (log_m_ia_max - log_m_ia_min) / n_ia_bins;
        
        for (int i = 0; i < n_ia_bins; i++) {
            double log_m = log_m_ia_min + (i + 0.5) * dlog_m_ia;
            double m = pow(10.0, log_m);
            double dm = m * log(10.0) * dlog_m_ia;
            stellar_mass_ia += total_mass * imf_norm * evaluate_imf(m, imf) * dm *m;
        }
        
        /* Apply DTD */
        sn_ia_rate = binary_fraction * stellar_mass_ia * type_ia_dtd(age);
    }
    
    /* Fill output structure */
    rates->sn_ia_rate = sn_ia_rate;
    rates->sn_ii_rate = sn_ii_rate;
    rates->pisn_rate = 0.0;  // Will be filled by PISN calculation if enabled
    rates->total_rate = sn_ia_rate + sn_ii_rate;
    
    return 0;
}

int calculate_energy_release(const SNRates *rates, EnergyRelease *energy) {
    if (rates == NULL || energy == NULL) return -1;
    
    /* Energy in erg per year */
    energy->sn_ia_energy = rates->sn_ia_rate * ENERGY_SN_IA * FOE;
    energy->sn_ii_energy = rates->sn_ii_rate * ENERGY_SN_II * FOE;
    energy->pisn_energy = 0.0;  // Will be filled by PISN calculation if enabled
    energy->total_energy = energy->sn_ia_energy + energy->sn_ii_energy;
    
    return 0;
}

int calculate_mass_ejection(const ClusterProperties *cluster,
                             const SNRates *rates,
                             MassEjection *mass) {
    if (cluster == NULL || rates == NULL || mass == NULL) return -1;
    
    /* Typical ejecta masses */
    /* Type Ia: ~1.4 M_sun (entire WD) */
    double ejecta_ia = 1.4;
    
    /* Type II: depends on mass, typically ~10-20 M_sun */
    /* Use average of ~15 M_sun for typical 20 M_sun progenitor */
    double ejecta_ii = 15.0;
    
    /* For continuous SF, also account for stellar winds and AGB mass loss */
    if (cluster->sf_mode == SF_CONTINUOUS) {
        /* In continuous SF, assume ~30% of formed mass is returned */
        /* This includes: stellar winds (massive stars), planetary nebulae (AGB), SNe */
        double sfr = cluster->star_formation_rate;
        double total_returned_mass = sfr * 0.30;  // 30% return fraction
        
        /* SN contribution */
        double sn_mass = rates->sn_ia_rate * ejecta_ia + rates->sn_ii_rate * ejecta_ii;
        
        /* Additional mass from winds/AGB */
        double wind_mass = total_returned_mass - sn_mass;
        if (wind_mass < 0) wind_mass = 0;
        
        mass->sn_ia_mass = rates->sn_ia_rate * ejecta_ia;
        mass->sn_ii_mass = rates->sn_ii_rate * ejecta_ii + wind_mass;
        mass->pisn_mass = 0.0;  // Will be filled by PISN calculation if enabled
        mass->total_mass = mass->sn_ia_mass + mass->sn_ii_mass;
    } else {
        /* Instantaneous burst: only SNe contribute */
        mass->sn_ia_mass = rates->sn_ia_rate * ejecta_ia;
        mass->sn_ii_mass = rates->sn_ii_rate * ejecta_ii;
        mass->pisn_mass = 0.0;  // Will be filled by PISN calculation if enabled
        mass->total_mass = mass->sn_ia_mass + mass->sn_ii_mass;
    }
    
    return 0;
}

double get_element_yield_value(ElementType elem,
                                int sn_type,
                                double mass,
                                double metallicity) {
    if (elem < 0 || elem >= N_ELEMENTS) return 0.0;
    
    double base_yield;
    
    if (sn_type == 0) {  /* Type Ia */
        base_yield = interpolate_metallicity(metallicity, yield_table[elem].sn_ia_yield);
    } else {  /* Type II */
        double yield_20M = interpolate_metallicity(metallicity, 
                                                     yield_table[elem].sn_ii_yield_20M);
        
        /* Scale for mass */
        if (mass <= 15.0) {
            base_yield = yield_20M * yield_table[elem].mass_scaling_15M;
        } else if (mass <= 20.0) {
            double frac = (mass - 15.0) / 5.0;
            double yield_15M = yield_20M * yield_table[elem].mass_scaling_15M;
            base_yield = yield_15M * (1.0 - frac) + yield_20M * frac;
        } else if (mass <= 25.0) {
            double frac = (mass - 20.0) / 5.0;
            double yield_25M = yield_20M * yield_table[elem].mass_scaling_25M;
            base_yield = yield_20M * (1.0 - frac) + yield_25M * frac;
        } else {
            double yield_25M = yield_20M * yield_table[elem].mass_scaling_25M;
            base_yield = yield_25M * (mass / 25.0);
        }
    }
    
    return base_yield;
}

int calculate_elemental_yields(const ClusterProperties *cluster,
                                const SNRates *rates,
                                ElementalYields *yields) {
    if (cluster == NULL || rates == NULL || yields == NULL) return -1;
    
    /* Initialize */
    for (int i = 0; i < N_ELEMENTS; i++) {
        yields->element_yields[i] = 0.0;
    }
    yields->total_metals = 0.0;
    
    /* Type Ia contribution (use 1.4 M_sun progenitor) */
    for (int i = 0; i < N_ELEMENTS; i++) {
        double yield_per_sn = get_element_yield_value(i, 0, CHANDRASEKHAR_MASS, 
                                                       cluster->metallicity);
        yields->element_yields[i] += rates->sn_ia_rate * yield_per_sn;
    }
    
    /* Type II contribution - IMF-weighted average over mass range */
    /* 
     * CRITICAL: Type II yields depend on progenitor mass!
     * Must integrate over 8-40 M_sun with IMF weighting
     * 
     * Yield_element = ∫ N(M) × Y(M,Z) dM / ∫ N(M) dM
     * where N(M) is number of stars from IMF
     */
    
    double imf_norm = calculate_imf_normalization(cluster->imf_type);
    
    /* Mass bins for Type II progenitors */
    int n_mass_bins = 50;
    double log_m_min = log10(MIN_MASS_SNII);
    double log_m_max = log10(MAX_MASS_SNII);
    double dlog_m = (log_m_max - log_m_min) / n_mass_bins;
    
    /* Calculate total number of Type II SNe for normalization */
    double total_snii_number = 0.0;
    for (int j = 0; j < n_mass_bins; j++) {
        double log_m = log_m_min + (j + 0.5) * dlog_m;
        double m = pow(10.0, log_m);
        double dm = m * log(10.0) * dlog_m;
        total_snii_number += imf_norm * evaluate_imf(m, cluster->imf_type) * dm;
    }
    
    /* Calculate IMF-weighted yields for each element */
    for (int i = 0; i < N_ELEMENTS; i++) {
        double weighted_yield = 0.0;
        
        /* Integrate over mass range */
        for (int j = 0; j < n_mass_bins; j++) {
            double log_m = log_m_min + (j + 0.5) * dlog_m;
            double m = pow(10.0, log_m);
            double dm = m * log(10.0) * dlog_m;
            
            /* Number of stars in this mass bin */
            double n_stars = imf_norm * evaluate_imf(m, cluster->imf_type) * dm;
            
            /* Yield for this mass */
            double yield_this_mass = get_element_yield_value(i, 1, m, 
                                                              cluster->metallicity);
            
            /* Weighted contribution */
            weighted_yield += n_stars * yield_this_mass;
        }
        
        /* Average yield per SN */
        double avg_yield_per_sn = weighted_yield / total_snii_number;
        
        /* Total yield rate */
        yields->element_yields[i] += rates->sn_ii_rate * avg_yield_per_sn;
    }
    
    /* Calculate total metals (everything except H and He) */
    for (int i = C_ELEM; i < N_ELEMENTS; i++) {
        yields->total_metals += yields->element_yields[i];
    }
    
    /* Calculate mean metallicity of ejected material */
    double total_ejected = 0.0;
    for (int i = 0; i < N_ELEMENTS; i++) {
        total_ejected += yields->element_yields[i];
    }
    
    if (total_ejected > 0.0) {
        yields->mean_metallicity = yields->total_metals / total_ejected;
    } else {
        yields->mean_metallicity = 0.0;
    }
    
    return 0;
}

int calculate_cluster_supernovae(const ClusterProperties *cluster,
                                  ClusterOutput *output) {
    if (cluster == NULL || output == NULL) return -1;
    
    /* Validate input */
    if (validate_cluster_properties(cluster) != 0) {
        fprintf(stderr, "Error: Invalid cluster properties\n");
        return -1;
    }
    
    /* Calculate all components */
    if (calculate_sn_rates(cluster, &output->rates) != 0) {
        fprintf(stderr, "Error: Failed to calculate SN rates\n");
        return -1;
    }
    
    if (calculate_energy_release(&output->rates, &output->energy) != 0) {
        fprintf(stderr, "Error: Failed to calculate energy release\n");
        return -1;
    }
    
    if (calculate_mass_ejection(cluster, &output->rates, &output->mass) != 0) {
        fprintf(stderr, "Error: Failed to calculate mass ejection\n");
        return -1;
    }
    
    if (calculate_elemental_yields(cluster, &output->rates, &output->yields) != 0) {
        fprintf(stderr, "Error: Failed to calculate elemental yields\n");
        return -1;
    }
    
    /* ===== Pop III / PISN Contribution (THREAD-SAFE) ===== */
    if (cluster->popiii_mode == POPIII_ENABLED) {
        /* THREAD-SAFE: All variables are local to this scope */
        double pisn_rate = 0.0;
        double pisn_energy = 0.0;
        double pisn_mass = 0.0;
        double pisn_yields[N_ELEMENTS];
        
        /* Calculate PISN contribution - THREAD-SAFE function */
        calculate_pisn_contribution(cluster, &pisn_rate, &pisn_energy, 
                                   &pisn_mass, pisn_yields);
        
        /* Add PISN to existing outputs - THREAD-SAFE: local accumulation */
        output->rates.pisn_rate = pisn_rate;
        output->rates.total_rate += pisn_rate;
        
        output->energy.pisn_energy = pisn_energy;
        output->energy.total_energy += pisn_energy;
        
        output->mass.pisn_mass = pisn_mass;
        output->mass.total_mass += pisn_mass;
        
        /* Add PISN yields to total */
        for (int i = 0; i < N_ELEMENTS; i++) {
            output->yields.element_yields[i] += pisn_yields[i];
        }
        
        /* Recalculate total metals and mean metallicity */
        output->yields.total_metals = 0.0;
        for (int i = C_ELEM; i < N_ELEMENTS; i++) {  /* C and heavier */
            output->yields.total_metals += output->yields.element_yields[i];
        }
        
        double total_ejecta = 0.0;
        for (int i = 0; i < N_ELEMENTS; i++) {
            total_ejecta += output->yields.element_yields[i];
        }
        
        if (total_ejecta > 0) {
            output->yields.mean_metallicity = output->yields.total_metals / total_ejecta;
        }
    }
    
    return 0;
}

/* ===== Utility Functions ===== */

const char* get_element_name(ElementType elem) {
    if (elem < 0 || elem >= N_ELEMENTS) return "Unknown";
    return yield_table[elem].name;
}

void print_cluster_properties(const ClusterProperties *cluster) {
    if (cluster == NULL) return;
    
    printf("=== Cluster Properties ===\n");
    printf("Total Mass:    %.2e M_sun\n", cluster->total_mass);
    printf("Age:           %.2e years (%.2f Myr)\n", cluster->age, cluster->age / 1e6);
    printf("Metallicity:   %.4f (%.2f Z_sun)\n", 
           cluster->metallicity, cluster->metallicity / SOLAR_Z);
    printf("IMF Type:      ");
    switch (cluster->imf_type) {
        case IMF_SALPETER: printf("Salpeter\n"); break;
        case IMF_KROUPA:   printf("Kroupa\n"); break;
        case IMF_CHABRIER: printf("Chabrier\n"); break;
        case IMF_POPIII_TOPHEAVY: printf("Pop III Top-Heavy\n"); break;
        default:           printf("Unknown\n");
    }
    printf("SF Mode:       ");
    switch (cluster->sf_mode) {
        case SF_INSTANTANEOUS: printf("Instantaneous burst\n"); break;
        case SF_CONTINUOUS:    
            printf("Continuous (SFR = %.2e M_sun/yr)\n", 
                   cluster->star_formation_rate); 
            break;
        default: printf("Unknown\n");
    }
    
    /* Pop III info */
    printf("Pop III:       ");
    if (cluster->popiii_mode == POPIII_ENABLED) {
        printf("Enabled (M_max = %.0f M_sun, mass fraction = %.1f%%)\n", 
               cluster->popiii_imf_max, 
               cluster->popiii_mass_fraction * 100.0);
    } else {
        printf("Disabled\n");
    }
    printf("\n");
}

void print_cluster_output(const ClusterOutput *output) {
    if (output == NULL) return;
    
    printf("=== Supernova Rates (per year) ===\n");
    printf("Type Ia:       %.6e\n", output->rates.sn_ia_rate);
    printf("Type II:       %.6e\n", output->rates.sn_ii_rate);
    if (output->rates.pisn_rate > 0) {
        printf("PISN (Pop III):%.6e\n", output->rates.pisn_rate);
    }
    printf("Total:         %.6e\n", output->rates.total_rate);
    printf("\n");
    
    printf("=== Energy Release (per year) ===\n");
    printf("Type Ia:       %.6e erg/yr (%.3g foe/yr)\n", 
           output->energy.sn_ia_energy, output->energy.sn_ia_energy / FOE);
    printf("Type II:       %.6e erg/yr (%.3g foe/yr)\n",
           output->energy.sn_ii_energy, output->energy.sn_ii_energy / FOE);
    if (output->energy.pisn_energy > 0) {
        printf("PISN (Pop III):%.6e erg/yr (%.3g foe/yr)\n",
               output->energy.pisn_energy, output->energy.pisn_energy / FOE);
    }
    printf("Total:         %.6e erg/yr (%.3g foe/yr)\n",
           output->energy.total_energy, output->energy.total_energy / FOE);
    printf("\n");
    
    printf("=== Mass Ejection (per year) ===\n");
    printf("Type Ia:       %.6e M_sun/yr\n", output->mass.sn_ia_mass);
    printf("Type II:       %.6e M_sun/yr\n", output->mass.sn_ii_mass);
    if (output->mass.pisn_mass > 0) {
        printf("PISN (Pop III):%.6e M_sun/yr\n", output->mass.pisn_mass);
    }
    printf("Total:         %.6e M_sun/yr\n", output->mass.total_mass);
    printf("\n");
    
    printf("=== Elemental Yields (per year) ===\n");
    printf("%-10s %15s\n", "Element", "Yield (M_sun/yr)");
    printf("----------------------------------------\n");
    for (int i = 0; i < N_ELEMENTS; i++) {
        printf("%-10s %15.6e\n", 
               get_element_name(i), output->yields.element_yields[i]);
    }
    printf("----------------------------------------\n");
    printf("Total Metals:  %.6e M_sun/yr\n", output->yields.total_metals);
    printf("Mean Z (ejecta): %.6f\n", output->yields.mean_metallicity);
    printf("\n");
}

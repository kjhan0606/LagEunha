/**
 * star_formation.h
 * 
 * Star formation module for galaxy formation simulations
 * Implements modern star formation models and IMF sampling
 * 
 * References:
 * - Kennicutt (1998): The Global Schmidt Law in Star-forming Galaxies, ApJ 498, 541
 * - Krumholz & McKee (2005): A General Theory of Turbulence-regulated Star Formation, ApJ 630, 250
 * - Padoan et al. (2012): The Star Formation Rate of Supersonic MHD Turbulence, ApJ 759, L27
 * - Chabrier (2003): Galactic Stellar and Substellar Initial Mass Function, PASP 115, 763
 * - Kroupa (2001): On the variation of the initial mass function, MNRAS 322, 231
 * - Kim et al. (2024): TIGRESS simulations (arXiv:2502.13244)
 * - Ostriker & Kim (2022): Pressure-Regulated, Feedback-Modulated star formation (ApJ 936, 137)
 */

#ifndef STAR_FORMATION_H
#define STAR_FORMATION_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

/* Define M_PI if not already defined */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Physical constants in CGS units */
#define GRAVITATIONAL_CONSTANT 6.674e-8     /* cm^3 g^-1 s^-2 */
#define BOLTZMANN_CONSTANT 1.381e-16        /* erg K^-1 */
#define PROTON_MASS 1.673e-24               /* g */
#define SOLAR_MASS 1.989e33                 /* g */
#define PARSEC 3.086e18                     /* cm */
#define YEAR 3.156e7                        /* s */
#define METALLICITY_SOLAR 0.0134            /* Solar metallicity (mass fraction) */

/* Star formation model parameters */
#define SF_EFFICIENCY_PER_FREEFALL 0.01     /* ε_ff ~ 1% (Krumholz & McKee 2005) */
#define KENNICUTT_SCHMIDT_INDEX 1.4         /* N in Σ_SFR ∝ Σ_gas^N (Kennicutt 1998) */
#define KENNICUTT_SCHMIDT_COEFF 2.5e-4      /* Normalization in M_sun yr^-1 kpc^-2 */

/* Density and temperature thresholds */
#define DENSITY_THRESHOLD_CGS 1.0e-22       /* g cm^-3, ~100 H atoms/cm^3 */
#define TEMPERATURE_THRESHOLD 1.0e4         /* K, threshold for cold gas */

/* IMF types */
typedef enum {
    IMF_CHABRIER_2003,
    IMF_KROUPA_2001,
    IMF_SALPETER_1955
} IMF_Type;

/* Structure to hold gas element properties */
typedef struct {
    double density;          /* Gas density in g cm^-3 */
    double temperature;      /* Gas temperature in K */
    double metallicity;      /* Metal mass fraction (Z/Z_sun) */
    double velocity_dispersion; /* 3D velocity dispersion in cm s^-1 */
    double volume;           /* Volume of gas element in cm^3 */
} GasElement;

/* Structure to hold star particle properties */
typedef struct {
    double mass;             /* Total stellar mass in M_sun */
    int num_stars;           /* Number of individual stars sampled */
    double *stellar_masses;  /* Array of individual stellar masses in M_sun */
    double velocity[3];      /* Velocity components in cm s^-1 */
    double metallicity;      /* Inherited metallicity */
    double age;              /* Age in years (typically 0 for new stars) */
} StarParticle;

/**
 * Calculate free-fall time for given density
 * 
 * t_ff = sqrt(3π / (32Gρ))
 * 
 * Reference: Eq. 2 in Krumholz & McKee (2005)
 * 
 * @param density Gas density in g cm^-3
 * @return Free-fall time in years
 */
double calculate_freefall_time(double density) {
    double t_ff_seconds = sqrt(3.0 * M_PI / (32.0 * GRAVITATIONAL_CONSTANT * density));
    return t_ff_seconds / YEAR;
}

/**
 * Calculate Jeans mass for given density and temperature
 * 
 * M_J = (5kT / GμmH)^(3/2) * (3 / 4πρ)^(1/2)
 * 
 * Reference: Standard ISM physics
 * 
 * @param density Gas density in g cm^-3
 * @param temperature Gas temperature in K
 * @param mean_molecular_weight Mean molecular weight (typically 1.3 for neutral gas)
 * @return Jeans mass in M_sun
 */
double calculate_jeans_mass(double density, double temperature, double mean_molecular_weight) {
    double numerator = pow(5.0 * BOLTZMANN_CONSTANT * temperature / 
                          (GRAVITATIONAL_CONSTANT * mean_molecular_weight * PROTON_MASS), 1.5);
    double denominator = sqrt(4.0 * M_PI * density / 3.0);
    double jeans_mass_grams = numerator / denominator;
    return jeans_mass_grams / SOLAR_MASS;
}

/**
 * Calculate turbulent star formation efficiency
 * Based on Padoan et al. (2012) model
 * 
 * ε_ff depends on virial parameter, Mach number, and density distribution
 * 
 * Reference: Padoan et al. (2012), ApJ 759, L27
 * 
 * @param mach_number Mach number (σ_v / c_s)
 * @param virial_parameter α_vir = 5σ_v^2 R / (GM)
 * @return Star formation efficiency per free-fall time
 */
double calculate_turbulent_sf_efficiency(double mach_number, double virial_parameter) {
    /* Simplified model - more complex dependencies can be added */
    /* Higher Mach number -> more turbulent support -> lower efficiency */
    /* Based on Padoan et al. (2012) multi-freefall model */
    
    double base_efficiency = SF_EFFICIENCY_PER_FREEFALL;
    double mach_factor = exp(-0.4 * mach_number);
    double virial_factor = exp(-0.5 * virial_parameter);
    
    return base_efficiency * mach_factor * virial_factor;
}

/**
 * Calculate star formation rate using Kennicutt-Schmidt law
 * 
 * Σ_SFR = A * Σ_gas^N
 * 
 * Reference: Kennicutt (1998), ApJ 498, 541
 *           Kennicutt & Evans (2012), ARA&A 50, 531
 * 
 * @param gas Gas element properties
 * @param timestep Time interval in years
 * @return Star formation rate in M_sun yr^-1
 */
double calculate_sfr_kennicutt_schmidt(const GasElement *gas, double timestep) {
    (void)timestep; /* Unused in this function */
    
    /* Only form stars in cold, dense gas */
    if (gas->temperature > TEMPERATURE_THRESHOLD || 
        gas->density < DENSITY_THRESHOLD_CGS) {
        return 0.0;
    }
    
    /* Convert volume density to surface density (approximation) */
    /* Assume characteristic length scale ~ (volume)^(1/3) */
    double length_scale_pc = pow(gas->volume, 1.0/3.0) / PARSEC;
    double gas_mass_msun = (gas->density * gas->volume) / SOLAR_MASS;
    double surface_density = gas_mass_msun / (length_scale_pc * length_scale_pc); /* M_sun pc^-2 */
    
    /* Kennicutt-Schmidt law: Σ_SFR = A * Σ_gas^N */
    /* Converting to kpc^-2 for standard units */
    double surface_density_kpc2 = surface_density * 1.0e6; /* M_sun kpc^-2 */
    double sfr_surface_density = KENNICUTT_SCHMIDT_COEFF * 
                                 pow(surface_density_kpc2, KENNICUTT_SCHMIDT_INDEX);
    
    /* Convert back to total SFR */
    double area_kpc2 = (length_scale_pc * length_scale_pc) / 1.0e6;
    double sfr = sfr_surface_density * area_kpc2;
    
    return sfr;
}

/**
 * Calculate star formation rate using volumetric approach
 * Based on free-fall time and star formation efficiency
 * 
 * SFR = ε_ff * M_gas / t_ff
 * 
 * Reference: Krumholz & McKee (2005), ApJ 630, 250
 * 
 * @param gas Gas element properties
 * @param timestep Time interval in years
 * @return Star formation rate in M_sun yr^-1
 */
double calculate_sfr_volumetric(const GasElement *gas, double timestep) {
    (void)timestep; /* Unused in this function */
    
    /* Only form stars in cold, dense gas */
    if (gas->temperature > TEMPERATURE_THRESHOLD || 
        gas->density < DENSITY_THRESHOLD_CGS) {
        return 0.0;
    }
    
    /* Calculate Mach number and virial parameter for turbulent efficiency */
    double sound_speed = sqrt(BOLTZMANN_CONSTANT * gas->temperature / 
                             (1.3 * PROTON_MASS)); /* cm s^-1 */
    double mach_number = gas->velocity_dispersion / sound_speed;
    
    /* Rough estimate of virial parameter */
    double length_scale = pow(gas->volume, 1.0/3.0);
    double gas_mass = gas->density * gas->volume;
    double virial_parameter = 5.0 * gas->velocity_dispersion * gas->velocity_dispersion * 
                             length_scale / (GRAVITATIONAL_CONSTANT * gas_mass);
    
    /* Calculate efficiency including turbulence effects */
    double epsilon_ff = calculate_turbulent_sf_efficiency(mach_number, virial_parameter);
    
    /* Calculate free-fall time */
    double t_ff = calculate_freefall_time(gas->density);
    
    /* Calculate gas mass */
    double gas_mass_msun = gas_mass / SOLAR_MASS;
    
    /* Star formation rate */
    double sfr = epsilon_ff * gas_mass_msun / t_ff;
    
    return sfr;
}

/**
 * Sample stellar mass from Chabrier (2003) IMF
 * 
 * ξ(m) = 0.158/m * exp[-(log10(m) - log10(0.08))^2 / (2*0.69^2)] for m < 1 M_sun
 * ξ(m) ∝ m^(-2.3) for m >= 1 M_sun
 * 
 * Reference: Chabrier (2003), PASP 115, 763
 * 
 * @return Stellar mass in M_sun
 */
double sample_chabrier_imf(void) {
    double u = (double)rand() / RAND_MAX;
    
    /* Approximate normalization: ~70% of stars are < 1 M_sun */
    if (u < 0.7) {
        /* Log-normal part for m < 1 M_sun */
        /* Using Box-Muller transform for Gaussian sampling */
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        
        double log_m = log10(0.08) + 0.69 * z;
        double m = pow(10.0, log_m);
        
        /* Clip to physical range */
        if (m < 0.01) m = 0.01;
        if (m > 1.0) m = 1.0;
        
        return m;
    } else {
        /* Power-law part for m >= 1 M_sun */
        /* m^(-2.3) -> cumulative: m^(-1.3) */
        double u_scaled = (double)rand() / RAND_MAX;
        double m_min = 1.0;
        double m_max = 100.0; /* Maximum stellar mass */
        
        /* Inverse transform sampling for power law */
        double alpha = 2.3;
        double m = pow(pow(m_min, 1.0 - alpha) + u_scaled * 
                      (pow(m_max, 1.0 - alpha) - pow(m_min, 1.0 - alpha)), 
                      1.0 / (1.0 - alpha));
        
        return m;
    }
}

/**
 * Sample stellar mass from Kroupa (2001) IMF
 * 
 * ξ(m) ∝ m^(-0.3) for 0.01 < m < 0.08 M_sun
 * ξ(m) ∝ m^(-1.3) for 0.08 < m < 0.5 M_sun
 * ξ(m) ∝ m^(-2.3) for m > 0.5 M_sun
 * 
 * Reference: Kroupa (2001), MNRAS 322, 231
 * 
 * @return Stellar mass in M_sun
 */
double sample_kroupa_imf(void) {
    double u = (double)rand() / RAND_MAX;
    
    /* Mass ranges with different slopes */
    double m_min = 0.01, m1 = 0.08, m2 = 0.5, m_max = 100.0;
    double alpha1 = 0.3, alpha2 = 1.3, alpha3 = 2.3;
    
    /* Compute cumulative weights (approximate) */
    double w1 = 0.15;  /* ~15% brown dwarfs */
    double w2 = 0.50;  /* ~50% low-mass stars */
    double w3 = 0.35;  /* ~35% higher-mass stars */
    
    double u_scaled = (double)rand() / RAND_MAX;
    double m;
    
    if (u < w1) {
        /* Brown dwarf regime */
        m = pow(pow(m_min, 1.0 - alpha1) + (u_scaled / w1) * 
               (pow(m1, 1.0 - alpha1) - pow(m_min, 1.0 - alpha1)), 
               1.0 / (1.0 - alpha1));
    } else if (u < w1 + w2) {
        /* Low-mass star regime */
        m = pow(pow(m1, 1.0 - alpha2) + ((u_scaled - w1) / w2) * 
               (pow(m2, 1.0 - alpha2) - pow(m1, 1.0 - alpha2)), 
               1.0 / (1.0 - alpha2));
    } else {
        /* High-mass star regime */
        m = pow(pow(m2, 1.0 - alpha3) + ((u_scaled - w1 - w2) / w3) * 
               (pow(m_max, 1.0 - alpha3) - pow(m2, 1.0 - alpha3)), 
               1.0 / (1.0 - alpha3));
    }
    
    return m;
}

/**
 * Sample stellar mass from Salpeter (1955) IMF
 * 
 * ξ(m) ∝ m^(-2.35)
 * 
 * Reference: Salpeter (1955), ApJ 121, 161
 * 
 * @return Stellar mass in M_sun
 */
double sample_salpeter_imf(void) {
    double u = (double)rand() / RAND_MAX;
    double m_min = 0.1;   /* Minimum mass */
    double m_max = 100.0; /* Maximum mass */
    double alpha = 2.35;
    
    /* Inverse transform sampling */
    double m = pow(pow(m_min, 1.0 - alpha) + u * 
                  (pow(m_max, 1.0 - alpha) - pow(m_min, 1.0 - alpha)), 
                  1.0 / (1.0 - alpha));
    
    return m;
}

/**
 * Create a star particle from a gas element
 * 
 * This function:
 * 1. Calculates the star formation rate
 * 2. Determines the stellar mass to form
 * 3. Samples individual stellar masses from the IMF
 * 4. Assigns velocities and metallicity
 * 
 * @param gas Input gas element
 * @param timestep Simulation timestep in years
 * @param imf_type Type of IMF to use
 * @param star_particle Output star particle (allocated by caller)
 * @return true if star formation occurred, false otherwise
 */
bool form_stars(const GasElement *gas, double timestep, IMF_Type imf_type, 
                StarParticle *star_particle) {
    /* Calculate SFR using preferred model */
    /* Can switch between Kennicutt-Schmidt and volumetric approaches */
    double sfr = calculate_sfr_volumetric(gas, timestep);
    
    /* Alternative: use Kennicutt-Schmidt law */
    /* double sfr = calculate_sfr_kennicutt_schmidt(gas, timestep); */
    
    if (sfr <= 0.0) {
        star_particle->mass = 0.0;
        star_particle->num_stars = 0;
        return false;
    }
    
    /* Calculate total stellar mass formed in this timestep */
    double stellar_mass_formed = sfr * timestep;
    
    /* Ensure we don't form more stars than available gas mass */
    double gas_mass_msun = (gas->density * gas->volume) / SOLAR_MASS;
    if (stellar_mass_formed > gas_mass_msun * 0.5) {
        stellar_mass_formed = gas_mass_msun * 0.5; /* Max 50% conversion per timestep */
    }
    
    /* Sample individual stellar masses from IMF */
    int max_stars = 10000; /* Maximum number of stars to sample */
    star_particle->stellar_masses = (double *)malloc(max_stars * sizeof(double));
    
    double total_sampled_mass = 0.0;
    int num_stars = 0;
    
    while (total_sampled_mass < stellar_mass_formed && num_stars < max_stars) {
        double stellar_mass;
        
        switch (imf_type) {
            case IMF_CHABRIER_2003:
                stellar_mass = sample_chabrier_imf();
                break;
            case IMF_KROUPA_2001:
                stellar_mass = sample_kroupa_imf();
                break;
            case IMF_SALPETER_1955:
                stellar_mass = sample_salpeter_imf();
                break;
            default:
                stellar_mass = sample_chabrier_imf();
        }
        
        /* Check if adding this star would exceed total mass */
        if (total_sampled_mass + stellar_mass > stellar_mass_formed * 1.1) {
            break;
        }
        
        star_particle->stellar_masses[num_stars] = stellar_mass;
        total_sampled_mass += stellar_mass;
        num_stars++;
    }
    
    /* Normalize to match exact total mass */
    double normalization = stellar_mass_formed / total_sampled_mass;
    for (int i = 0; i < num_stars; i++) {
        star_particle->stellar_masses[i] *= normalization;
    }
    
    /* Set star particle properties */
    star_particle->mass = stellar_mass_formed;
    star_particle->num_stars = num_stars;
    
    /* Inherit velocity from gas (with small random perturbation) */
    /* For simplicity, we assume gas has bulk velocity ~0 and only velocity dispersion */
    for (int i = 0; i < 3; i++) {
        double random_gaussian = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        star_particle->velocity[i] = gas->velocity_dispersion * random_gaussian / sqrt(3.0);
    }
    
    /* Inherit metallicity */
    star_particle->metallicity = gas->metallicity;
    
    /* New stars have age = 0 */
    star_particle->age = 0.0;
    
    return true;
}

/**
 * Free memory allocated for star particle
 * 
 * @param star_particle Star particle to free
 */
void free_star_particle(StarParticle *star_particle) {
    if (star_particle->stellar_masses != NULL) {
        free(star_particle->stellar_masses);
        star_particle->stellar_masses = NULL;
    }
}

/**
 * Calculate mean stellar mass for a given IMF
 * Useful for statistics and mass budgeting
 * 
 * @param imf_type Type of IMF
 * @return Mean stellar mass in M_sun
 */
double calculate_mean_stellar_mass(IMF_Type imf_type) {
    int n_samples = 100000;
    double total_mass = 0.0;
    
    for (int i = 0; i < n_samples; i++) {
        double m;
        switch (imf_type) {
            case IMF_CHABRIER_2003:
                m = sample_chabrier_imf();
                break;
            case IMF_KROUPA_2001:
                m = sample_kroupa_imf();
                break;
            case IMF_SALPETER_1955:
                m = sample_salpeter_imf();
                break;
            default:
                m = sample_chabrier_imf();
        }
        total_mass += m;
    }
    
    return total_mass / n_samples;
}

#endif /* STAR_FORMATION_H */

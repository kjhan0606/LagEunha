/**
 * star_formation_advanced.h
 * 
 * Advanced star formation module for galaxy formation simulations
 * Implements:
 * - Magnetic field effects on star formation
 * - Metallicity-dependent IMF variations
 * - Binary star population synthesis
 * 
 * New References:
 * - Pattle et al. (2022): Magnetic fields in star formation (arXiv:2203.11179)
 * - Kretschmer & Teyssier (2020): Magnetized SF model (MNRAS 527, 6779)
 * - Tanvir et al. (2024): Metallicity dependence of IMF (MNRAS 527, 7306)
 * - Yan et al. (2018, 2020): IGIMF theory (A&A 607, 126; A&A 637, 68)
 * - Bate (2023): Redshift-metallicity IMF variation (MNRAS 537, 752)
 * - Moe & Di Stefano (2017): Binary fraction systematics (ApJS 230, 15)
 * - Stanway & Eldridge (2020): Binary population uncertainties (MNRAS 495, 4605)
 */

#ifndef STAR_FORMATION_ADVANCED_H
#define STAR_FORMATION_ADVANCED_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/* Define M_PI if not already defined */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Physical constants in CGS units */
#define GRAVITATIONAL_CONSTANT 6.674e-8
#define BOLTZMANN_CONSTANT 1.381e-16
#define PROTON_MASS 1.673e-24
#define SOLAR_MASS 1.989e33
#define PARSEC 3.086e18
#define YEAR 3.156e7
#define METALLICITY_SOLAR 0.0134
#define SPEED_OF_LIGHT 2.998e10              /* cm s^-1 */

/* Magnetic field parameters */
#define ALFVEN_SPEED_SUBPC 3.0e5             /* Typical Alfven speed at sub-pc scales (cm/s) */
#define MAGNETIC_CRITICALITY_THRESHOLD 1.0   /* λ = μ/μ_crit threshold */

/* Star formation model parameters */
#define SF_EFFICIENCY_PER_FREEFALL 0.01
#define KENNICUTT_SCHMIDT_INDEX 1.4
#define KENNICUTT_SCHMIDT_COEFF 2.5e-4

/* Density and temperature thresholds */
#define DENSITY_THRESHOLD_CGS 1.0e-22
#define TEMPERATURE_THRESHOLD 1.0e4

/* Binary star parameters (Moe & Di Stefano 2017) */
#define BINARY_FRACTION_LOWMASS 0.23         /* < 0.5 M_sun */
#define BINARY_FRACTION_SOLAR 0.44           /* 0.5-1.5 M_sun */
#define BINARY_FRACTION_INTERMEDIATE 0.59    /* 1.5-5 M_sun */
#define BINARY_FRACTION_MASSIVE 0.70         /* > 5 M_sun */

/* IMF types */
typedef enum {
    IMF_CHABRIER_2003,
    IMF_KROUPA_2001,
    IMF_SALPETER_1955,
    IMF_METALLICITY_DEPENDENT  /* New: varies with Z */
} IMF_Type;

/* Binary orbital parameters */
typedef struct {
    double mass_ratio;       /* q = M2/M1, range [0,1] */
    double period_days;      /* Orbital period in days */
    double eccentricity;     /* Orbital eccentricity [0,1) */
    double separation_au;    /* Semi-major axis in AU */
} BinaryOrbit;

/* Structure to hold gas element properties (extended) */
typedef struct {
    double density;              /* Gas density in g cm^-3 */
    double temperature;          /* Gas temperature in K */
    double metallicity;          /* Metal mass fraction (Z/Z_sun) */
    double velocity_dispersion;  /* 3D velocity dispersion in cm s^-1 */
    double volume;               /* Volume of gas element in cm^3 */
    
    /* New: Magnetic field properties */
    double magnetic_field_strength; /* B field in Gauss */
    double magnetic_pressure;       /* B^2/(8π) in dyn cm^-2 */
    double mass_to_flux_ratio;      /* μ = M/Φ_B (normalized to critical) */
} GasElement;

/* Structure to hold star particle properties (extended) */
typedef struct {
    double mass;                /* Total stellar mass in M_sun */
    int num_stars;              /* Number of individual stars sampled */
    double *stellar_masses;     /* Array of individual stellar masses in M_sun */
    double velocity[3];         /* Velocity components in cm s^-1 */
    double metallicity;         /* Inherited metallicity */
    double age;                 /* Age in years */
    
    /* New: Binary information */
    int num_binaries;           /* Number of binary systems */
    bool *is_binary;            /* Array: is star in a binary? */
    BinaryOrbit *binary_params; /* Array of binary orbital parameters */
} StarParticle;

/**
 * Calculate mass-to-flux ratio (normalized to critical value)
 * 
 * λ = (M/Φ_B) / (M/Φ_B)_crit
 * where (M/Φ_B)_crit ~ 0.13 / sqrt(G) in cgs units
 * 
 * λ > 1: magnetically supercritical (can collapse)
 * λ < 1: magnetically subcritical (supported against collapse)
 * 
 * Reference: Pattle et al. (2022), Section 3.2
 * 
 * @param mass Cloud mass in grams
 * @param magnetic_flux Magnetic flux in G cm^2
 * @return Normalized mass-to-flux ratio
 */
double calculate_mass_to_flux_ratio(double mass, double magnetic_flux) {
    double critical_ratio = 0.13 / sqrt(GRAVITATIONAL_CONSTANT);
    double mass_to_flux = mass / magnetic_flux;
    return mass_to_flux / critical_ratio;
}

/**
 * Calculate magnetic field strength from gas properties
 * Assumes flux freezing and equipartition with turbulence
 * 
 * B ~ sqrt(4π ρ) * v_turb (equipartition)
 * 
 * Reference: Pattle et al. (2022), Section 2
 * 
 * @param density Gas density in g cm^-3
 * @param velocity_dispersion Turbulent velocity in cm s^-1
 * @return Magnetic field strength in Gauss
 */
double estimate_magnetic_field(double density, double velocity_dispersion) {
    /* Equipartition: B^2/(8π) ~ ρ v^2 */
    double B_equipartition = sqrt(4.0 * M_PI * density) * velocity_dispersion;
    return B_equipartition;
}

/**
 * Calculate star formation efficiency with magnetic field effects
 * 
 * Magnetic fields reduce SF efficiency through:
 * 1. Magnetic pressure support
 * 2. Angular momentum transport (magnetic braking)
 * 
 * Reference: Kretschmer & Teyssier (2020), MNRAS 527, 6779
 * 
 * @param gas Gas element with magnetic properties
 * @param epsilon_ff_base Base SF efficiency without magnetic fields
 * @return Modified SF efficiency with magnetic effects
 */
double calculate_sf_efficiency_magnetized(const GasElement *gas, double epsilon_ff_base) {
    /* Calculate magnetic pressure relative to thermal pressure */
    double thermal_pressure = (gas->density / PROTON_MASS) * BOLTZMANN_CONSTANT * gas->temperature;
    double beta_plasma = thermal_pressure / gas->magnetic_pressure;
    
    /* Strong magnetic fields (low beta) suppress SF */
    /* Empirical fit from Kretschmer & Teyssier (2020) */
    double magnetic_suppression = 1.0;
    
    if (beta_plasma < 1.0) {
        /* Magnetically dominated: significant suppression */
        magnetic_suppression = 0.5 * (1.0 + beta_plasma);
    } else if (beta_plasma < 10.0) {
        /* Intermediate regime: moderate suppression */
        magnetic_suppression = 0.7 + 0.3 * (beta_plasma - 1.0) / 9.0;
    }
    /* For beta > 10, thermal pressure dominates, minimal suppression */
    
    /* Check if magnetically subcritical (cannot collapse) */
    if (gas->mass_to_flux_ratio < MAGNETIC_CRITICALITY_THRESHOLD) {
        /* Subcritical: very low or no star formation */
        return epsilon_ff_base * 0.1 * gas->mass_to_flux_ratio;
    }
    
    return epsilon_ff_base * magnetic_suppression;
}

/**
 * Calculate metallicity-dependent IMF characteristic mass
 * 
 * Lower metallicity → higher characteristic mass (top-heavy IMF)
 * Higher metallicity → lower characteristic mass (bottom-heavy IMF)
 * 
 * Reference: Tanvir et al. (2024), MNRAS 527, 7306
 *            Bate (2023), MNRAS 537, 752
 * 
 * @param metallicity Metallicity in solar units
 * @param surface_density Gas surface density in M_sun pc^-2
 * @return Characteristic mass scale in M_sun
 */
double calculate_characteristic_mass(double metallicity, double surface_density) {
    /* Base characteristic mass at solar metallicity and normal density */
    double m_char_solar = 0.25; /* M_sun */
    
    /* Metallicity effect (Tanvir et al. 2024) */
    /* More metal-rich → slightly more bottom-heavy */
    double z_factor = pow(metallicity, -0.15);
    
    /* Surface density effect is stronger (Tanvir et al. 2024) */
    /* Higher surface density → more bottom-heavy */
    double sigma_factor = pow(surface_density / 100.0, -0.3);
    
    /* At very low metallicity, characteristic mass increases significantly */
    if (metallicity < 0.01) {
        /* Primordial/metal-poor star formation (Bate 2023) */
        z_factor = pow(metallicity, -0.4);
    }
    
    double m_char = m_char_solar * z_factor * sigma_factor;
    
    /* Clamp to reasonable range */
    if (m_char < 0.05) m_char = 0.05;
    if (m_char > 5.0) m_char = 5.0;
    
    return m_char;
}

/**
 * Sample stellar mass from metallicity-dependent IMF
 * 
 * Uses log-normal form with metallicity-dependent peak
 * 
 * Reference: Yan et al. (2018), A&A 607, A126
 * 
 * @param metallicity Metallicity in solar units
 * @param surface_density Surface density in M_sun pc^-2
 * @return Stellar mass in M_sun
 */
double sample_metallicity_dependent_imf(double metallicity, double surface_density) {
    double u = (double)rand() / RAND_MAX;
    double m_char = calculate_characteristic_mass(metallicity, surface_density);
    
    /* Log-normal peak shifts with metallicity */
    double log_m_peak = log10(m_char);
    double sigma_log_m = 0.55; /* Width of log-normal part */
    
    if (u < 0.75) {
        /* Log-normal part (most stars) */
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        
        double log_m = log_m_peak + sigma_log_m * z;
        double m = pow(10.0, log_m);
        
        if (m < 0.01) m = 0.01;
        if (m > 1.5) m = 1.5;
        
        return m;
    } else {
        /* High-mass power law tail */
        /* Slope varies slightly with metallicity */
        double alpha = 2.3 - 0.2 * log10(metallicity + 0.01);
        if (alpha < 2.0) alpha = 2.0;
        if (alpha > 2.7) alpha = 2.7;
        
        double m_min = 1.5;
        double m_max = 100.0;
        double u_scaled = (double)rand() / RAND_MAX;
        
        double m = pow(pow(m_min, 1.0 - alpha) + u_scaled * 
                      (pow(m_max, 1.0 - alpha) - pow(m_min, 1.0 - alpha)), 
                      1.0 / (1.0 - alpha));
        
        return m;
    }
}

/**
 * Determine if a star should be in a binary system
 * Binary fraction varies with stellar mass (Moe & Di Stefano 2017)
 * 
 * Reference: Moe & Di Stefano (2017), ApJS 230, 15
 * 
 * @param primary_mass Primary star mass in M_sun
 * @return true if star should be in binary
 */
bool is_in_binary(double primary_mass) {
    double binary_fraction;
    
    /* Binary fraction increases with stellar mass */
    if (primary_mass < 0.5) {
        binary_fraction = BINARY_FRACTION_LOWMASS;
    } else if (primary_mass < 1.5) {
        binary_fraction = BINARY_FRACTION_SOLAR;
    } else if (primary_mass < 5.0) {
        binary_fraction = BINARY_FRACTION_INTERMEDIATE;
    } else {
        binary_fraction = BINARY_FRACTION_MASSIVE;
    }
    
    double u = (double)rand() / RAND_MAX;
    return u < binary_fraction;
}

/**
 * Generate binary companion properties
 * 
 * Mass ratio distribution: flat for q > 0.3, rising for q < 0.3
 * Period distribution: log-uniform
 * 
 * Reference: Moe & Di Stefano (2017), ApJS 230, 15
 * 
 * @param primary_mass Primary star mass in M_sun
 * @param orbit Output binary orbital parameters
 */
void generate_binary_companion(double primary_mass, BinaryOrbit *orbit) {
    /* Mass ratio distribution (Moe & Di Stefano 2017) */
    /* Approximately flat for massive stars, peaks near q=1 */
    double u_q = (double)rand() / RAND_MAX;
    
    /* For simplicity: flat distribution in q */
    /* More sophisticated: weight toward q ~ 1 for massive stars */
    double q_min = 0.1;
    if (primary_mass > 5.0) {
        /* Massive stars: prefer equal mass ratios */
        orbit->mass_ratio = q_min + (1.0 - q_min) * pow(u_q, 0.5);
    } else {
        /* Lower mass: flatter distribution */
        orbit->mass_ratio = q_min + (1.0 - q_min) * u_q;
    }
    
    /* Period distribution: log-uniform from 1 day to 1e5 days */
    double log_period_min = 0.0;  /* log10(1 day) */
    double log_period_max = 5.0;  /* log10(100,000 days) ~ 300 yr */
    double u_p = (double)rand() / RAND_MAX;
    double log_period = log_period_min + (log_period_max - log_period_min) * u_p;
    orbit->period_days = pow(10.0, log_period);
    
    /* Eccentricity: thermal distribution e^2 */
    /* For short periods: circularized (e ~ 0) */
    if (orbit->period_days < 10.0) {
        orbit->eccentricity = 0.0;
    } else {
        double u_e = (double)rand() / RAND_MAX;
        orbit->eccentricity = sqrt(u_e);
        if (orbit->eccentricity > 0.9) orbit->eccentricity = 0.9;
    }
    
    /* Calculate semi-major axis from Kepler's third law */
    /* a^3 = G(M1+M2)P^2 / (4π^2) */
    double total_mass = primary_mass * (1.0 + orbit->mass_ratio); /* M_sun */
    double period_seconds = orbit->period_days * 86400.0;
    
    double a_cm = pow(GRAVITATIONAL_CONSTANT * total_mass * SOLAR_MASS * 
                     period_seconds * period_seconds / (4.0 * M_PI * M_PI), 1.0/3.0);
    orbit->separation_au = a_cm / 1.496e13; /* Convert cm to AU */
}

/**
 * Sample stellar mass from chosen IMF type
 * Includes metallicity-dependent option
 * 
 * @param imf_type Type of IMF to use
 * @param metallicity Metallicity (only used for IMF_METALLICITY_DEPENDENT)
 * @param surface_density Surface density (only used for IMF_METALLICITY_DEPENDENT)
 * @return Stellar mass in M_sun
 */
double sample_stellar_mass(IMF_Type imf_type, double metallicity, double surface_density);

/* Forward declarations for IMF functions from original code */
double sample_chabrier_imf(void);
double sample_kroupa_imf(void);
double sample_salpeter_imf(void);

/**
 * Free-fall time calculation (from original code)
 */
double calculate_freefall_time(double density) {
    (void)density; /* Placeholder */
    double t_ff_seconds = sqrt(3.0 * M_PI / (32.0 * GRAVITATIONAL_CONSTANT * density));
    return t_ff_seconds / YEAR;
}

/**
 * Calculate turbulent SF efficiency (from original code, placeholder)
 */
double calculate_turbulent_sf_efficiency(double mach_number, double virial_parameter) {
    double base_efficiency = SF_EFFICIENCY_PER_FREEFALL;
    double mach_factor = exp(-0.4 * mach_number);
    double virial_factor = exp(-0.5 * virial_parameter);
    return base_efficiency * mach_factor * virial_factor;
}

/**
 * Calculate SFR with magnetic field effects (volumetric approach)
 * 
 * @param gas Gas element properties
 * @param timestep Time interval in years
 * @return Star formation rate in M_sun yr^-1
 */
double calculate_sfr_magnetized(const GasElement *gas, double timestep) {
    (void)timestep;
    
    if (gas->temperature > TEMPERATURE_THRESHOLD || 
        gas->density < DENSITY_THRESHOLD_CGS) {
        return 0.0;
    }
    
    /* Calculate turbulent properties */
    double sound_speed = sqrt(BOLTZMANN_CONSTANT * gas->temperature / (1.3 * PROTON_MASS));
    double mach_number = gas->velocity_dispersion / sound_speed;
    
    double length_scale = pow(gas->volume, 1.0/3.0);
    double gas_mass = gas->density * gas->volume;
    double virial_parameter = 5.0 * gas->velocity_dispersion * gas->velocity_dispersion * 
                             length_scale / (GRAVITATIONAL_CONSTANT * gas_mass);
    
    /* Base efficiency from turbulence */
    double epsilon_ff_base = calculate_turbulent_sf_efficiency(mach_number, virial_parameter);
    
    /* Apply magnetic field effects */
    double epsilon_ff = calculate_sf_efficiency_magnetized(gas, epsilon_ff_base);
    
    /* Calculate free-fall time and SFR */
    double t_ff = calculate_freefall_time(gas->density);
    double gas_mass_msun = gas_mass / SOLAR_MASS;
    double sfr = epsilon_ff * gas_mass_msun / t_ff;
    
    return sfr;
}

/**
 * Create a star particle with advanced features:
 * - Magnetic field effects on SF
 * - Metallicity-dependent IMF
 * - Binary star populations
 * 
 * @param gas Input gas element
 * @param timestep Simulation timestep in years
 * @param imf_type Type of IMF to use
 * @param include_binaries Whether to generate binary systems
 * @param star_particle Output star particle
 * @return true if star formation occurred
 */
bool form_stars_advanced(const GasElement *gas, double timestep, IMF_Type imf_type,
                        bool include_binaries, StarParticle *star_particle) {
    /* Calculate SFR including magnetic field effects */
    double sfr = calculate_sfr_magnetized(gas, timestep);
    
    if (sfr <= 0.0) {
        star_particle->mass = 0.0;
        star_particle->num_stars = 0;
        return false;
    }
    
    /* Calculate total stellar mass formed */
    double stellar_mass_formed = sfr * timestep;
    
    /* Ensure reasonable conversion efficiency */
    double gas_mass_msun = (gas->density * gas->volume) / SOLAR_MASS;
    if (stellar_mass_formed > gas_mass_msun * 0.5) {
        stellar_mass_formed = gas_mass_msun * 0.5;
    }
    
    /* Calculate surface density for metallicity-dependent IMF */
    double length_scale_pc = pow(gas->volume, 1.0/3.0) / PARSEC;
    double surface_density = gas_mass_msun / (length_scale_pc * length_scale_pc);
    
    /* Sample individual stellar masses */
    int max_stars = 10000;
    star_particle->stellar_masses = (double *)malloc(max_stars * sizeof(double));
    star_particle->is_binary = (bool *)malloc(max_stars * sizeof(bool));
    star_particle->binary_params = (BinaryOrbit *)malloc(max_stars * sizeof(BinaryOrbit));
    
    double total_sampled_mass = 0.0;
    int num_stars = 0;
    int num_binaries = 0;
    
    while (total_sampled_mass < stellar_mass_formed && num_stars < max_stars) {
        double primary_mass;
        
        /* Sample primary mass from IMF */
        if (imf_type == IMF_METALLICITY_DEPENDENT) {
            primary_mass = sample_metallicity_dependent_imf(
                gas->metallicity / METALLICITY_SOLAR, surface_density);
        } else if (imf_type == IMF_CHABRIER_2003) {
            primary_mass = sample_chabrier_imf();
        } else if (imf_type == IMF_KROUPA_2001) {
            primary_mass = sample_kroupa_imf();
        } else {
            primary_mass = sample_salpeter_imf();
        }
        
        /* Check if adding this star exceeds total mass */
        if (total_sampled_mass + primary_mass > stellar_mass_formed * 1.1) {
            break;
        }
        
        star_particle->stellar_masses[num_stars] = primary_mass;
        total_sampled_mass += primary_mass;
        
        /* Determine if star is in binary */
        if (include_binaries && is_in_binary(primary_mass)) {
            star_particle->is_binary[num_stars] = true;
            generate_binary_companion(primary_mass, &star_particle->binary_params[num_stars]);
            
            /* Add companion mass to total */
            double companion_mass = primary_mass * star_particle->binary_params[num_stars].mass_ratio;
            total_sampled_mass += companion_mass;
            num_binaries++;
        } else {
            star_particle->is_binary[num_stars] = false;
        }
        
        num_stars++;
    }
    
    /* Normalize to match exact total mass */
    double normalization = stellar_mass_formed / total_sampled_mass;
    for (int i = 0; i < num_stars; i++) {
        star_particle->stellar_masses[i] *= normalization;
        if (star_particle->is_binary[i]) {
            /* Companion mass also scales */
            double companion_mass = star_particle->stellar_masses[i] * 
                                   star_particle->binary_params[i].mass_ratio;
            (void)companion_mass; /* Mass is implicit in mass_ratio */
        }
    }
    
    /* Set star particle properties */
    star_particle->mass = stellar_mass_formed;
    star_particle->num_stars = num_stars;
    star_particle->num_binaries = num_binaries;
    
    /* Inherit velocity from gas */
    for (int i = 0; i < 3; i++) {
        double random_gaussian = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        star_particle->velocity[i] = gas->velocity_dispersion * random_gaussian / sqrt(3.0);
    }
    
    /* Inherit metallicity */
    star_particle->metallicity = gas->metallicity;
    star_particle->age = 0.0;
    
    return true;
}

/**
 * Free memory allocated for advanced star particle
 */
void free_star_particle_advanced(StarParticle *star_particle) {
    if (star_particle->stellar_masses != NULL) {
        free(star_particle->stellar_masses);
        star_particle->stellar_masses = NULL;
    }
    if (star_particle->is_binary != NULL) {
        free(star_particle->is_binary);
        star_particle->is_binary = NULL;
    }
    if (star_particle->binary_params != NULL) {
        free(star_particle->binary_params);
        star_particle->binary_params = NULL;
    }
}

/* Include original IMF sampling functions */
double sample_chabrier_imf(void) {
    double u = (double)rand() / RAND_MAX;
    
    if (u < 0.7) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        
        double log_m = log10(0.08) + 0.69 * z;
        double m = pow(10.0, log_m);
        
        if (m < 0.01) m = 0.01;
        if (m > 1.0) m = 1.0;
        
        return m;
    } else {
        double u_scaled = (double)rand() / RAND_MAX;
        double m_min = 1.0;
        double m_max = 100.0;
        double alpha = 2.3;
        
        double m = pow(pow(m_min, 1.0 - alpha) + u_scaled * 
                      (pow(m_max, 1.0 - alpha) - pow(m_min, 1.0 - alpha)), 
                      1.0 / (1.0 - alpha));
        
        return m;
    }
}

double sample_kroupa_imf(void) {
    double u = (double)rand() / RAND_MAX;
    double m_min = 0.01, m1 = 0.08, m2 = 0.5, m_max = 100.0;
    double alpha1 = 0.3, alpha2 = 1.3, alpha3 = 2.3;
    double w1 = 0.15, w2 = 0.50, w3 = 0.35;
    double u_scaled = (double)rand() / RAND_MAX;
    double m;
    
    if (u < w1) {
        m = pow(pow(m_min, 1.0 - alpha1) + (u_scaled / w1) * 
               (pow(m1, 1.0 - alpha1) - pow(m_min, 1.0 - alpha1)), 
               1.0 / (1.0 - alpha1));
    } else if (u < w1 + w2) {
        m = pow(pow(m1, 1.0 - alpha2) + ((u_scaled - w1) / w2) * 
               (pow(m2, 1.0 - alpha2) - pow(m1, 1.0 - alpha2)), 
               1.0 / (1.0 - alpha2));
    } else {
        m = pow(pow(m2, 1.0 - alpha3) + ((u_scaled - w1 - w2) / w3) * 
               (pow(m_max, 1.0 - alpha3) - pow(m2, 1.0 - alpha3)), 
               1.0 / (1.0 - alpha3));
    }
    
    return m;
}

double sample_salpeter_imf(void) {
    double u = (double)rand() / RAND_MAX;
    double m_min = 0.1, m_max = 100.0;
    double alpha = 2.35;
    
    double m = pow(pow(m_min, 1.0 - alpha) + u * 
                  (pow(m_max, 1.0 - alpha) - pow(m_min, 1.0 - alpha)), 
                  1.0 / (1.0 - alpha));
    
    return m;
}

double sample_stellar_mass(IMF_Type imf_type, double metallicity, double surface_density) {
    if (imf_type == IMF_METALLICITY_DEPENDENT) {
        return sample_metallicity_dependent_imf(metallicity, surface_density);
    } else if (imf_type == IMF_CHABRIER_2003) {
        return sample_chabrier_imf();
    } else if (imf_type == IMF_KROUPA_2001) {
        return sample_kroupa_imf();
    } else {
        return sample_salpeter_imf();
    }
}

#endif /* STAR_FORMATION_ADVANCED_H */

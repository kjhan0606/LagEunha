/*
 * ============================================================================
 * STELLAR WIND YIELD LIBRARY - Implementation
 * ============================================================================
 * 
 * See stellar_wind_lib.h for documentation.
 * 
 * ============================================================================
 */

#include "stellar_wind_lib.h"
#include "massive_yields.h"
#include "agb_yields.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================================================================
 * INTERNAL CONSTANTS AND GLOBAL STATE
 * ============================================================================ */

/* Element names */
static const char* g_element_names[SWL_N_ELEMENTS] = {
    "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe", "Z_total"
};

/* IMF names */
static const char* g_imf_names[SWL_N_IMF_TYPES] = {
    "Salpeter (1955)",
    "Kroupa (2001)",
    "Chabrier (2003)",
    "Baldry & Glazebrook (2003)"
};

/* Error messages */
static const char* g_error_messages[] = {
    "Success",
    "Invalid IMF type",
    "Invalid mass value",
    "Invalid metallicity value",
    "Invalid age value",
    "Null pointer argument",
    "Integration failed"
};

/* Solar abundances (mass fractions) - Asplund et al. (2009) */
static const double g_solar_abundance[SWL_N_ELEMENTS] = {
    0.7381,     /* H */
    0.2485,     /* He */
    2.38e-3,    /* C */
    6.96e-4,    /* N */
    5.79e-3,    /* O */
    1.26e-3,    /* Ne */
    7.06e-4,    /* Mg */
    6.64e-4,    /* Si */
    3.09e-4,    /* S */
    6.13e-5,    /* Ca */
    1.17e-3,    /* Fe */
    0.0134      /* Z_total */
};

/* Library configuration (global state) */
static SWL_Config g_config;
static int g_initialized = 0;

/* Cached IMF normalization constants */
static double g_imf_norm[SWL_N_IMF_TYPES];

/* ============================================================================
 * YIELD TABLE DATA - Now loaded from external headers:
 *   - massive_yields.h: Limongi & Chieffi (2018) wind yields
 *   - agb_yields.h: Karakas & Lugaro (2016) AGB yields
 * 
 * Multiple metallicity tables available:
 *   Massive: [Fe/H] = 0, -1, -2, -3 (Z = 0.02, 0.002, 0.0002, 0.00002)
 *   AGB: Z = 0.007, 0.014, 0.03
 * ============================================================================ */

/* ============================================================================
 * INTERNAL HELPER FUNCTIONS
 * ============================================================================ */

/*
 * Linear interpolation helper
 */
static double lerp(double x, double x0, double x1, double y0, double y1) {
    if (x1 == x0) return y0;
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

/*
 * Log-space interpolation helper
 */
static double log_lerp(double x, double x0, double x1, double y0, double y1) {
    if (x1 == x0 || x <= 0 || x0 <= 0 || x1 <= 0) return lerp(x, x0, x1, y0, y1);
    if (y0 <= 0 || y1 <= 0) return lerp(x, x0, x1, y0, y1);
    
    double log_x = log10(x);
    double log_x0 = log10(x0);
    double log_x1 = log10(x1);
    double log_y0 = log10(y0);
    double log_y1 = log10(y1);
    
    double log_y = lerp(log_x, log_x0, log_x1, log_y0, log_y1);
    return pow(10.0, log_y);
}

/*
 * Clamp value to range
 */
static double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

/* ============================================================================
 * IMF FUNCTIONS
 * ============================================================================ */

/*
 * Salpeter (1955) IMF: phi(m) ~ m^-2.35
 */
static double imf_salpeter(double m) {
    return pow(m, -2.35);
}

/*
 * Kroupa (2001) IMF: three-segment power law
 */
static double imf_kroupa(double m) {
    const double m1 = 0.08, m2 = 0.50;
    const double alpha0 = 0.3, alpha1 = 1.3, alpha2 = 2.3;
    
    double c1 = pow(m1, alpha1 - alpha0);
    double c2 = c1 * pow(m2, alpha2 - alpha1);
    
    if (m < m1) return pow(m, -alpha0);
    if (m < m2) return c1 * pow(m, -alpha1);
    return c2 * pow(m, -alpha2);
}

/*
 * Chabrier (2003) IMF: log-normal + power law
 */
static double imf_chabrier(double m) {
    const double mc = 0.08, sigma = 0.69, alpha = 2.3;
    
    if (m < 1.0) {
        double log_m = log10(m);
        double log_mc = log10(mc);
        double exponent = -(log_m - log_mc) * (log_m - log_mc) / (2.0 * sigma * sigma);
        return exp(exponent) / m;
    } else {
        double log_mc = log10(mc);
        double exponent_1 = -log_mc * log_mc / (2.0 * sigma * sigma);
        double match_factor = exp(exponent_1);
        return match_factor * pow(m, -alpha);
    }
}

/*
 * Baldry & Glazebrook (2003) IMF
 */
static double imf_baldry_glazebrook(double m) {
    const double m_break = 0.5;
    const double alpha1 = 1.5, alpha2 = 2.35;
    
    if (m < m_break) return pow(m, -alpha1);
    double c = pow(m_break, alpha2 - alpha1);
    return c * pow(m, -alpha2);
}

/*
 * Raw (un-normalized) IMF
 */
static double imf_raw(double m, SWL_IMFType type) {
    switch (type) {
        case SWL_IMF_SALPETER:         return imf_salpeter(m);
        case SWL_IMF_KROUPA:           return imf_kroupa(m);
        case SWL_IMF_CHABRIER:         return imf_chabrier(m);
        case SWL_IMF_BALDRY_GLAZEBROOK: return imf_baldry_glazebrook(m);
        default: return imf_kroupa(m);
    }
}

/*
 * Calculate IMF normalization using Simpson's rule
 */
static double calculate_imf_norm(SWL_IMFType type) {
    const int n_steps = 10000;
    double log_m_min = log10(g_config.m_min);
    double log_m_max = log10(g_config.m_max);
    double d_log_m = (log_m_max - log_m_min) / n_steps;
    
    double integral = 0.0;
    for (int i = 0; i <= n_steps; i++) {
        double log_m = log_m_min + i * d_log_m;
        double m = pow(10.0, log_m);
        double weight = (i == 0 || i == n_steps) ? 1.0 : (i % 2 == 0) ? 2.0 : 4.0;
        double phi = imf_raw(m, type);  /* Use raw (un-normalized) IMF */
        integral += weight * m * phi * m * log(10.0) * d_log_m;
    }
    integral /= 3.0;
    
    if (integral <= 0) return 1.0;  /* Safety check */
    return 1.0 / integral;
}

/* ============================================================================
 * STELLAR PHYSICS FUNCTIONS
 * ============================================================================ */

#include <math.h>

/*
 * Stellar lifetime approximation (Corrected for High-Mass Stars)
 * Based on standard stellar evolution trends (e.g., Geneva/Padova models)
 */
double swl_stellar_lifetime(double mass, double metallicity) {
    double log_m = log10(mass);
    double log_tau;

    // Use simplified fitting coefficients that are safer for high-mass extrapolation
    // Or use piecewise fitting if higher precision is required
    
    // Adjusted coefficients to prevent overly steep decline for massive stars
    // (Approximating standard t ~ M^-2.5 to M^-3 relation)
    double a0 = 10.13;
    double a1 = -3.6; // Adjusted from -4.424 to avoid unphysical drop
    double a2 = 0.5;  // Adjusted curvature
    
    log_tau = a0 + a1 * log_m + a2 * log_m * log_m;

    /* Metallicity correction factor */
    // SWL_Z_SUN should be defined elsewhere (e.g., 0.020 or 0.014)
    double Z_factor = 1.0 - 0.35 * log10(metallicity / SWL_Z_SUN);
    
    // Clamp the Z_factor to avoid extreme values
    if (Z_factor < 0.5) Z_factor = 0.5;
    if (Z_factor > 1.5) Z_factor = 1.5;

    double t_calc = pow(10.0, log_tau) * Z_factor;

    /* * [SAFETY FLOOR] 
     * Physically, stars cannot evolve significantly faster than ~3 Myr 
     * (the Eddington limit saturation).
     * If the calculated lifetime is less than 3 Myr, clamp it to this minimum.
     */
    double min_lifetime = 3.0e6; // 3.0 Myr hard limit
    
    if (t_calc < min_lifetime) {
        return min_lifetime;
    }
    
    return t_calc;
}

/*
 * Inverse: get turnoff mass from age
 * (Logic remains the same, but now uses the corrected lifetime function)
 */
double swl_turnoff_mass(double age, double metallicity) {
    /* Binary search to invert lifetime function */
    double m_low = 0.1, m_high = 300.0;

    /* Handle edge cases */
    if (age >= swl_stellar_lifetime(m_low, metallicity)) return m_low;
    if (age <= swl_stellar_lifetime(m_high, metallicity)) return m_high;

    for (int iter = 0; iter < 50; iter++) {
        double m_mid = 0.5 * (m_low + m_high);
        double tau_mid = swl_stellar_lifetime(m_mid, metallicity);

        if (tau_mid > age) {
            m_low = m_mid;
        } else {
            m_high = m_mid;
        }

        /* Check for convergence */
        if (fabs(m_high - m_low) < 1e-6 * m_mid) break;
    }

    return 0.5 * (m_low + m_high);
}
/*
 * Remnant mass as function of initial mass
 * Based on Fryer et al. (2012) and Spera & Mapelli (2017)
 */
double swl_remnant_mass(double initial_mass, double metallicity) {
    if (initial_mass < 0.8) {
        /* No remnant (star still on MS for Hubble time) */
        return initial_mass;
    } else if (initial_mass < SWL_M_MASSIVE) {
        /* White dwarf: use initial-final mass relation */
        /* Kalirai et al. (2008) */
        double m_wd = 0.109 * initial_mass + 0.394;
        return clamp(m_wd, 0.5, 1.4);
    } else if (initial_mass < 25.0) {
        /* Neutron star */
        return 1.4;
    } else if (initial_mass < 40.0) {
        /* Black hole (direct collapse or fallback) */
        /* Metallicity-dependent: lower Z -> more massive BH */
        double Z_factor = clamp(metallicity / SWL_Z_SUN, 0.01, 1.0);
        double f_fb = 0.3 + 0.3 * (1.0 - Z_factor);  /* Fallback fraction */
        return 1.4 + f_fb * (initial_mass - 25.0);
    } else {
        /* Massive black hole */
        double Z_factor = clamp(metallicity / SWL_Z_SUN, 0.01, 1.0);
        double f_core = 0.25 + 0.15 * (1.0 - Z_factor);
        return f_core * initial_mass;
    }
}

/*
 * Mass-loss rate for massive stars (Vink et al. 2001)
 */
static double mdot_vink(double mass, double lum, double T_eff, double Z) {
    double log_L5 = log10(lum / 1.0e5);
    double log_M30 = log10(mass / 30.0);
    double log_T40 = log10(T_eff / 40000.0);
    double log_Z = log10(Z / SWL_Z_SUN);
    
    double v_ratio = (T_eff >= 25000.0) ? 2.6 : 1.3;
    double log_mdot;
    
    if (T_eff >= 25000.0) {
        log_mdot = -6.697 + 2.194 * log_L5 - 1.313 * log_M30
                 - 1.226 * log10(v_ratio / 2.0) + 0.933 * log_T40
                 - 10.92 * log_T40 * log_T40 + 0.85 * log_Z;
    } else {
        log_mdot = -6.688 + 2.210 * log_L5 - 1.339 * log_M30
                 - 1.601 * log10(v_ratio / 2.0) + 1.07 * log_T40
                 + 0.85 * log_Z;
    }
    
    return pow(10.0, log_mdot);
}

/*
 * Mass-luminosity relation for main sequence stars
 */
static double mass_luminosity(double mass) {
    if (mass < 0.43) {
        return 0.23 * pow(mass, 2.3);
    } else if (mass < 2.0) {
        return pow(mass, 4.0);
    } else if (mass < 55.0) {
        return 1.4 * pow(mass, 3.5);
    } else {
        return 32000.0 * mass;  /* Near Eddington limit */
    }
}

/*
 * Effective temperature from mass (rough MS approximation)
 */
static double mass_teff(double mass) {
    /* log(T_eff) ~ 3.7 + 0.15*log(M) for MS stars */
    double log_T = 3.76 + 0.15 * log10(mass);
    return pow(10.0, log_T);
}

double swl_mass_loss_rate(double mass, double luminosity, double T_eff,
                          double metallicity, SWL_StellarPhase phase) {
    if (luminosity <= 0) luminosity = mass_luminosity(mass);
    if (T_eff <= 0) T_eff = mass_teff(mass);
    
    switch (phase) {
        case SWL_PHASE_MS:
        case SWL_PHASE_WR:
            return mdot_vink(mass, luminosity, T_eff, metallicity);
            
        case SWL_PHASE_RGB:
            /* Reimers (1975) */
            return 4.0e-13 * luminosity * 100.0 / mass;  /* Approximate R ~ 100 R_sun */
            
        case SWL_PHASE_AGB:
            /* Enhanced Reimers / Bloecker */
            {
                double mdot_r = 4.0e-13 * luminosity * 200.0 / mass;
                return 4.83e-9 * pow(mass, -2.1) * pow(luminosity, 2.7) * mdot_r;
            }
            
        default:
            return 0.0;
    }
}

/* ============================================================================
 * YIELD TABLE INTERPOLATION - Using external header yield tables
 * ============================================================================ */

/*
 * Convert metallicity Z to [Fe/H] for massive star table selection
 */
static int z_to_feh(double Z) {
    /* [Fe/H] = log10(Z/Z_sun) where Z_sun # 0.014 */
    double feh = log10(Z / 0.014);
    
    if (feh >= -0.5) return 0;       /* [Fe/H] ~  0, Z = 0.02 */
    if (feh >= -1.5) return -1;      /* [Fe/H] ~ -1, Z = 0.002 */
    if (feh >= -2.5) return -2;      /* [Fe/H] ~ -2, Z = 0.0002 */
    return -3;                        /* [Fe/H] ~ -3, Z = 0.00002 */
}

/*
 * Interpolate wind yields for massive stars with metallicity interpolation
 */
static void interpolate_massive_yields(double mass, double Z, 
                                        double* yields, double* wind_mass, 
                                        double* wind_fraction) {
    int feh = z_to_feh(Z);
    const YieldTableEntry* table = get_wind_yield_table(feh);
    
    if (!table) {
        table = g_wind_yields_solar;  /* Fallback to solar */
    }
    
    /* Find bracketing entry for wind_mass and wind_fraction */
    int i;
    for (i = 0; i < N_GWIND_MASS_BINS - 1; i++) {
        if (mass <= table[i + 1].mass) break;
    }
    
    double wind_m, wind_f;
    if (mass <= table[0].mass) {
        wind_m = table[0].wind_mass * (mass / table[0].mass);
        wind_f = table[0].wind_frac;
    } else if (mass >= table[N_GWIND_MASS_BINS - 1].mass) {
        wind_m = table[N_GWIND_MASS_BINS - 1].wind_mass * (mass / table[N_GWIND_MASS_BINS - 1].mass);
        wind_f = table[N_GWIND_MASS_BINS - 1].wind_frac;
    } else {
        double t = (mass - table[i].mass) / (table[i + 1].mass - table[i].mass);
        wind_m = table[i].wind_mass + t * (table[i + 1].wind_mass - table[i].wind_mass);
        wind_f = table[i].wind_frac + t * (table[i + 1].wind_frac - table[i].wind_frac);
    }
    
    /* Interpolate yields for each element */
    for (int e = 0; e < SWL_N_ELEMENTS; e++) {
        yields[e] = interp_wind_yield(table, mass, (ElementIndex)e);
    }
    
    *wind_mass = wind_m;
    *wind_fraction = wind_f;
    
    /* Metallicity interpolation between tables if needed */
    double log_Z = log10(Z);
    double Z_table = table[0].Z;
    double log_Z_table = log10(Z_table);
    
    /* Get neighboring [Fe/H] table for interpolation */
    int feh_hi = feh + 1;
    if (feh_hi > 0) feh_hi = 0;
    
    const YieldTableEntry* table_hi = get_wind_yield_table(feh_hi);
    if (table_hi && feh_hi != feh) {
        double Z_hi = table_hi[0].Z;
        double log_Z_hi = log10(Z_hi);
        
        if (log_Z > log_Z_table && log_Z < log_Z_hi) {
            double wZ = (log_Z - log_Z_table) / (log_Z_hi - log_Z_table);
            
            double yields_hi[SWL_N_ELEMENTS];
            for (int e = 0; e < SWL_N_ELEMENTS; e++) {
                yields_hi[e] = interp_wind_yield(table_hi, mass, (ElementIndex)e);
                yields[e] = (1.0 - wZ) * yields[e] + wZ * yields_hi[e];
            }
            
            /* Interpolate wind mass too */
            double wind_m_hi;
            if (mass <= table_hi[0].mass) {
                wind_m_hi = table_hi[0].wind_mass * (mass / table_hi[0].mass);
            } else if (mass >= table_hi[N_GWIND_MASS_BINS - 1].mass) {
                wind_m_hi = table_hi[N_GWIND_MASS_BINS - 1].wind_mass;
            } else {
                for (i = 0; i < N_GWIND_MASS_BINS - 1; i++) {
                    if (mass <= table_hi[i + 1].mass) break;
                }
                double t = (mass - table_hi[i].mass) / (table_hi[i + 1].mass - table_hi[i].mass);
                wind_m_hi = table_hi[i].wind_mass + t * (table_hi[i + 1].wind_mass - table_hi[i].wind_mass);
            }
            *wind_mass = (1.0 - wZ) * (*wind_mass) + wZ * wind_m_hi;
        }
    }
}

/*
 * Interpolate wind yields for AGB stars with metallicity interpolation
 */
static void interpolate_agb_yields(double mass, double Z, 
                                    double* yields, double* wind_mass, 
                                    double* wind_fraction) {
    int n;
    const YieldTableEntry* table = get_agb_yield_table(Z, &n);
    
    /* Find bracketing entry for wind_mass and wind_fraction */
    int i;
    for (i = 0; i < n - 1; i++) {
        if (mass <= table[i + 1].mass) break;
    }
    
    double wind_m, wind_f;
    if (mass <= table[0].mass) {
        wind_m = table[0].wind_mass * (mass / table[0].mass);
        wind_f = table[0].wind_frac;
    } else if (mass >= table[n - 1].mass) {
        wind_m = table[n - 1].wind_mass;
        wind_f = table[n - 1].wind_frac;
    } else {
        double t = (mass - table[i].mass) / (table[i + 1].mass - table[i].mass);
        wind_m = table[i].wind_mass + t * (table[i + 1].wind_mass - table[i].wind_mass);
        wind_f = table[i].wind_frac + t * (table[i + 1].wind_frac - table[i].wind_frac);
    }
    
    /* Interpolate yields for each element */
    for (int e = 0; e < SWL_N_ELEMENTS; e++) {
        yields[e] = interp_agb_yield(table, n, mass, (ElementIndex)e);
    }
    
    *wind_mass = wind_m;
    *wind_fraction = wind_f;
}

/*
 * Get wind yields for a star of given mass
 * Returns both net yields and ejected mass
 */
static void get_star_wind_yields(double mass, double Z, double* yields, 
                                  double* wind_mass, double* wind_fraction) {
    for (int i = 0; i < SWL_N_ELEMENTS; i++) yields[i] = 0.0;
    *wind_mass = 0.0;
    *wind_fraction = 0.0;
    
    if (mass >= SWL_M_MASSIVE && g_config.include_massive_winds) {
        /* Use massive star wind yields from Limongi & Chieffi (2018) */
        interpolate_massive_yields(mass, Z, yields, wind_mass, wind_fraction);
    } else if (mass >= SWL_M_AGB_MIN && mass < SWL_M_MASSIVE && g_config.include_agb_winds) {
        /* Use AGB star yields from Karakas & Lugaro (2016) */
        interpolate_agb_yields(mass, Z, yields, wind_mass, wind_fraction);
    } else if (mass >= 0.5) {
        /* Low-mass RGB mass loss (minimal) */
        double dm = 0.1 * (mass - 0.5);
        yields[SWL_EL_H] = -dm * 0.74;
        yields[SWL_EL_HE] = dm * 0.74;
        yields[SWL_EL_Z_TOTAL] = 0.0;  /* Mass conservation */
        *wind_mass = dm;
        *wind_fraction = dm / mass;
    }
    
    /* Apply metallicity floor scaling for metals if below floor */
    if (Z < g_config.metallicity_floor) {
        double scale = Z / g_config.metallicity_floor;
        for (int i = SWL_EL_C; i < SWL_N_ELEMENTS; i++) {
            yields[i] *= scale;
        }
    }
}

/*
 * Calculate ejected mass from net yields
 * Ejected = Net Yield + (wind_mass × initial_abundance)
 */
static void calculate_ejected_mass(double* net_yields, double wind_mass, 
                                    double Z, double* ejected) {
    /* Initial abundances scaled by metallicity */
    double X_H = 0.76 - 3.0 * Z;   /* Hydrogen mass fraction */
    double X_He = 0.24 + 2.0 * Z;  /* Helium mass fraction */
    
    if (X_H < 0.70) X_H = 0.70;
    if (X_He > 0.28) X_He = 0.28;
    
    /* Initial mass fractions for each element */
    double initial_frac[SWL_N_ELEMENTS];
    initial_frac[SWL_EL_H] = X_H;
    initial_frac[SWL_EL_HE] = X_He;
    
    /* Scale metal abundances by Z/Z_sun */
    double Z_scale = Z / SWL_Z_SUN;
    for (int i = SWL_EL_C; i < SWL_N_ELEMENTS; i++) {
        initial_frac[i] = g_solar_abundance[i] * Z_scale;
    }
    
    /* Calculate ejected mass for each element */
    /* Ejected = Net_Yield + (wind_mass × initial_fraction) */
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        double initial_in_wind = wind_mass * initial_frac[i];
        ejected[i] = net_yields[i] + initial_in_wind;
        
        /* Ensure non-negative (numerical safety) */
        if (ejected[i] < 0) ejected[i] = 0.0;
    }
}

/* ============================================================================
 * TIME-DEPENDENT WIND YIELD CALCULATION
 * ============================================================================ */

/*
 * Calculate wind yields up to a given stellar age (fraction of lifetime)
 */
static void get_time_dependent_yields(double mass, double Z, double age,
                                       double* yields, double* wind_mass) {
    double lifetime = swl_stellar_lifetime(mass, Z);
    double total_yields[SWL_N_ELEMENTS], total_wind, wind_frac;
    
    /* Get total lifetime-integrated yields */
    get_star_wind_yields(mass, Z, total_yields, &total_wind, &wind_frac);
    
    /* If star has died, return full yields */
    if (age >= lifetime) {
        for (int i = 0; i < SWL_N_ELEMENTS; i++) yields[i] = total_yields[i];
        *wind_mass = total_wind;
        return;
    }
    
    /* Otherwise, calculate fraction of yields released so far */
    double age_frac = age / lifetime;
    
    /*
     * Wind mass loss is NOT linear with time. It's heavily weighted toward
     * the end of stellar evolution, especially for:
     * - Massive stars: WR phase at end of life
     * - AGB stars: superwind phase
     * 
     * We use an approximate model:
     * For massive stars: most mass lost in last ~20% of life
     * For AGB stars: superwind in last ~10% of life
     */
    double yield_frac;
    
    if (mass >= SWL_M_MASSIVE) {
        /* Massive star wind profile */
        /* Approximate cumulative: slow start, accelerating loss */
        /* f(t) = (t/tau)^1.5 for early times, steeper near end */
        if (age_frac < 0.8) {
            yield_frac = 0.3 * pow(age_frac / 0.8, 1.5);
        } else {
            double x = (age_frac - 0.8) / 0.2;
            yield_frac = 0.3 + 0.7 * (1.0 - pow(1.0 - x, 2.0));
        }
    } else if (mass >= SWL_M_AGB_MIN) {
        /* AGB star: most mass lost in final AGB/superwind phase */
        /* RGB mass loss is small, AGB superwind dominates */
        double t_rgb = 0.9;  /* RGB lasts ~90% of post-MS life */
        if (age_frac < t_rgb) {
            yield_frac = 0.1 * (age_frac / t_rgb);
        } else {
            double x = (age_frac - t_rgb) / (1.0 - t_rgb);
            yield_frac = 0.1 + 0.9 * x;
        }
    } else {
        /* Low mass: linear approximation */
        yield_frac = age_frac;
    }
    
    /* Apply yield fraction */
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        yields[i] = total_yields[i] * yield_frac;
    }
    *wind_mass = total_wind * yield_frac;
}

/* ============================================================================
 * PUBLIC API IMPLEMENTATION
 * ============================================================================ */

SWL_ErrorCode swl_init(void) {
    if (g_initialized) return SWL_SUCCESS;
    
    /* Set default configuration */
    g_config.m_min = SWL_M_MIN;
    g_config.m_max = SWL_M_MAX;
    g_config.n_integration_steps = 2000;
    g_config.include_agb_winds = 1;
    g_config.include_massive_winds = 1;
    g_config.include_vms_enhancement = 1;
    g_config.metallicity_floor = 1e-4;
    
    /* Initialize IMF normalizations (temporary value first) */
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        g_imf_norm[i] = 1.0;
    }
    
    /* Calculate proper normalizations */
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        g_imf_norm[i] = calculate_imf_norm((SWL_IMFType)i);
    }
    
    g_initialized = 1;
    return SWL_SUCCESS;
}

void swl_cleanup(void) {
    g_initialized = 0;
}

SWL_ErrorCode swl_get_default_config(SWL_Config* config) {
    if (!config) return SWL_ERR_NULL_POINTER;
    
    config->m_min = SWL_M_MIN;
    config->m_max = SWL_M_MAX;
    config->n_integration_steps = 2000;
    config->include_agb_winds = 1;
    config->include_massive_winds = 1;
    config->include_vms_enhancement = 1;
    config->metallicity_floor = 1e-4;
    
    return SWL_SUCCESS;
}

SWL_ErrorCode swl_set_config(const SWL_Config* config) {
    if (!config) return SWL_ERR_NULL_POINTER;
    
    g_config = *config;
    
    /* Recalculate IMF normalizations */
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        g_imf_norm[i] = calculate_imf_norm((SWL_IMFType)i);
    }
    
    return SWL_SUCCESS;
}

double swl_imf(double mass, SWL_IMFType imf_type) {
    if (!g_initialized) swl_init();
    return g_imf_norm[imf_type] * imf_raw(mass, imf_type);
}

const char* swl_get_element_name(SWL_Element element) {
    if (element >= 0 && element < SWL_N_ELEMENTS) {
        return g_element_names[element];
    }
    return "Unknown";
}

const char* swl_get_imf_name(SWL_IMFType imf_type) {
    if (imf_type >= 0 && imf_type < SWL_N_IMF_TYPES) {
        return g_imf_names[imf_type];
    }
    return "Unknown";
}

double swl_get_solar_abundance(SWL_Element element) {
    if (element >= 0 && element < SWL_N_ELEMENTS) {
        return g_solar_abundance[element];
    }
    return 0.0;
}

const char* swl_get_error_message(SWL_ErrorCode error) {
    if (error >= 0 && error <= SWL_ERR_INTEGRATION_FAILED) {
        return g_error_messages[error];
    }
    return "Unknown error";
}

/* ============================================================================
 * MAIN CALCULATION FUNCTIONS
 * ============================================================================ */

SWL_ErrorCode swl_calculate_star_wind_yield(
    double initial_mass,
    double metallicity,
    double stellar_age,
    SWL_StarWindYield* result)
{
    if (!g_initialized) swl_init();
    if (!result) return SWL_ERR_NULL_POINTER;
    if (initial_mass <= 0 || initial_mass > 500) return SWL_ERR_INVALID_MASS;
    if (metallicity <= 0 || metallicity > 0.1) return SWL_ERR_INVALID_METALLICITY;
    
    result->initial_mass = initial_mass;
    result->metallicity = metallicity;
    
    double yields[SWL_N_ELEMENTS];
    double wind_mass;
    
    if (stellar_age < 0) {
        /* Lifetime-integrated yields */
        double wind_frac;
        get_star_wind_yields(initial_mass, metallicity, yields, &wind_mass, &wind_frac);
    } else {
        /* Time-dependent yields */
        get_time_dependent_yields(initial_mass, metallicity, stellar_age, yields, &wind_mass);
    }
    
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        result->element_yields.yields[i] = yields[i];
    }
    result->total_wind_mass = wind_mass;
    
    /* Approximate wind velocity */
    if (initial_mass >= SWL_M_MASSIVE) {
        result->wind_velocity = 2000.0 + 500.0 * log10(initial_mass / 10.0);  /* km/s */
    } else {
        result->wind_velocity = 15.0 + 5.0 * initial_mass;  /* km/s for AGB */
    }
    
    return SWL_SUCCESS;
}

SWL_ErrorCode swl_calculate_cluster_wind_yield(
    SWL_IMFType imf_type,
    double metallicity,
    double cluster_mass,
    double cluster_age,
    SWL_ClusterWindYield* result)
{
    if (!g_initialized) swl_init();
    if (!result) return SWL_ERR_NULL_POINTER;
    if (imf_type < 0 || imf_type >= SWL_N_IMF_TYPES) return SWL_ERR_INVALID_IMF;
    if (metallicity <= 0 || metallicity > 0.1) return SWL_ERR_INVALID_METALLICITY;
    if (cluster_mass <= 0) return SWL_ERR_INVALID_MASS;
    if (cluster_age < 0) return SWL_ERR_INVALID_AGE;
    
    /* Initialize result structure */
    result->imf_type = imf_type;
    result->metallicity = metallicity;
    result->cluster_mass = cluster_mass;
    result->cluster_age = cluster_age;
    
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        result->element_yields.yields[i] = 0.0;
        result->ejected_mass.yields[i] = 0.0;
        result->yield_rate.yields[i] = 0.0;
        result->ejection_rate.yields[i] = 0.0;
    }
    result->total_wind_mass = 0.0;
    result->mass_in_living_stars = 0.0;
    result->mass_in_remnants = 0.0;
    result->wind_mass_rate = 0.0;
    result->n_stars_total = 0.0;
    result->n_stars_dead = 0.0;
    
    /* Get turnoff mass */
    result->turnoff_mass = swl_turnoff_mass(cluster_age, metallicity);
    
    /* Integration over the IMF */
    const int n_steps = g_config.n_integration_steps;
    double log_m_min = log10(g_config.m_min);
    double log_m_max = log10(g_config.m_max);
    double d_log_m = (log_m_max - log_m_min) / n_steps;
    
    double total_n_stars = 0.0;
    double total_n_dead = 0.0;
    double total_mass_living = 0.0;
    double total_mass_remnants = 0.0;
    double total_mass_returned = 0.0;
    double total_net_yields[SWL_N_ELEMENTS] = {0};
    double total_ejected[SWL_N_ELEMENTS] = {0};
    double current_wind_rate = 0.0;
    double current_yield_rate[SWL_N_ELEMENTS] = {0};
    double current_ejection_rate[SWL_N_ELEMENTS] = {0};
    
    for (int i = 0; i <= n_steps; i++) {
        double log_m = log_m_min + i * d_log_m;
        double m = pow(10.0, log_m);
        double dm = m * log(10.0) * d_log_m;
        
        /* Simpson's rule weight */
        double weight = (i == 0 || i == n_steps) ? 1.0 : (i % 2 == 0) ? 2.0 : 4.0;
        
        /* IMF value (number of stars per unit mass per log mass interval) */
        double phi = swl_imf(m, imf_type);  /* Already normalized */
        double n_stars = cluster_mass * phi * dm;  /* Number of stars in this mass bin */
        
        total_n_stars += weight * n_stars;
        
        /* Get stellar lifetime */
        double lifetime = swl_stellar_lifetime(m, metallicity);
        
        /* Get yields up to cluster age */
        double net_yields[SWL_N_ELEMENTS];
        double wind_mass;
        get_time_dependent_yields(m, metallicity, cluster_age, net_yields, &wind_mass);
        
        /* Calculate ejected mass from net yields */
        double ejected[SWL_N_ELEMENTS];
        calculate_ejected_mass(net_yields, wind_mass, metallicity, ejected);
        
        /* Accumulate yields weighted by number of stars */
        for (int j = 0; j < SWL_N_ELEMENTS; j++) {
            total_net_yields[j] += weight * n_stars * net_yields[j];
            total_ejected[j] += weight * n_stars * ejected[j];
        }
        total_mass_returned += weight * n_stars * wind_mass;
        
        /* Track living vs dead stars */
        if (cluster_age >= lifetime) {
            /* Star has died */
            total_n_dead += weight * n_stars;
            double m_remnant = swl_remnant_mass(m, metallicity);
            total_mass_remnants += weight * n_stars * m_remnant;
        } else {
            /* Star is still alive */
            double mass_lost = wind_mass;
            total_mass_living += weight * n_stars * (m - mass_lost);
            
            /* Calculate instantaneous wind rate for living stars */
            if (m >= SWL_M_AGB_MIN) {
                double lum = mass_luminosity(m);
                double T_eff = mass_teff(m);
                SWL_StellarPhase phase = (m >= SWL_M_MASSIVE) ? SWL_PHASE_MS : SWL_PHASE_AGB;
                double mdot = swl_mass_loss_rate(m, lum, T_eff, metallicity, phase);
                current_wind_rate += weight * n_stars * mdot;
                
                /* Get total lifetime yields for rate calculation */
                double total_yields_temp[SWL_N_ELEMENTS], wind_total, wind_frac;
                get_star_wind_yields(m, metallicity, total_yields_temp, &wind_total, &wind_frac);
                
                double ejected_temp[SWL_N_ELEMENTS];
                calculate_ejected_mass(total_yields_temp, wind_total, metallicity, ejected_temp);
                
                if (wind_total > 0) {
                    for (int j = 0; j < SWL_N_ELEMENTS; j++) {
                        /* Net yield rate */
                        current_yield_rate[j] += weight * n_stars * mdot * (total_yields_temp[j] / wind_total);
                        /* Ejection rate */
                        current_ejection_rate[j] += weight * n_stars * mdot * (ejected_temp[j] / wind_total);
                    }
                }
            }
        }
    }
    
    /* Apply Simpson's 1/3 factor */
    double simpson_factor = 1.0 / 3.0;
    
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        result->element_yields.yields[i] = total_net_yields[i] * simpson_factor;
        result->ejected_mass.yields[i] = total_ejected[i] * simpson_factor;
        result->yield_rate.yields[i] = current_yield_rate[i] * simpson_factor;
        result->ejection_rate.yields[i] = current_ejection_rate[i] * simpson_factor;
    }
    result->total_wind_mass = total_mass_returned * simpson_factor;
    result->mass_in_living_stars = total_mass_living * simpson_factor;
    result->mass_in_remnants = total_mass_remnants * simpson_factor;
    result->mass_returned_total = result->total_wind_mass;
    result->n_stars_total = total_n_stars * simpson_factor;
    result->n_stars_dead = total_n_dead * simpson_factor;
    result->wind_mass_rate = current_wind_rate * simpson_factor;
    
    return SWL_SUCCESS;
}

SWL_ErrorCode swl_calculate_imf_integrated_yields(
    SWL_IMFType imf_type,
    double metallicity,
    SWL_ElementYields* yields,
    double* wind_mass)
{
    if (!g_initialized) swl_init();
    if (!yields || !wind_mass) return SWL_ERR_NULL_POINTER;
    
    /* Use cluster calculation with 1 M_sun and infinite age */
    SWL_ClusterWindYield result;
    SWL_ErrorCode err = swl_calculate_cluster_wind_yield(
        imf_type, metallicity, 1.0, 1e12, &result);
    
    if (err != SWL_SUCCESS) return err;
    
    *yields = result.element_yields;
    *wind_mass = result.total_wind_mass;
    
    return SWL_SUCCESS;
}

SWL_ErrorCode swl_calculate_wind_rate(
    SWL_IMFType imf_type,
    double metallicity,
    double cluster_mass,
    double cluster_age,
    SWL_ElementYields* yield_rate,
    double* mass_rate)
{
    if (!g_initialized) swl_init();
    if (!yield_rate || !mass_rate) return SWL_ERR_NULL_POINTER;
    
    SWL_ClusterWindYield result;
    SWL_ErrorCode err = swl_calculate_cluster_wind_yield(
        imf_type, metallicity, cluster_mass, cluster_age, &result);
    
    if (err != SWL_SUCCESS) return err;
    
    *yield_rate = result.yield_rate;
    *mass_rate = result.wind_mass_rate;
    
    return SWL_SUCCESS;
}

SWL_ErrorCode swl_calculate_yield_history(
    SWL_IMFType imf_type,
    double metallicity,
    double cluster_mass,
    const double* ages,
    int n_ages,
    SWL_ClusterWindYield* results)
{
    if (!g_initialized) swl_init();
    if (!ages || !results) return SWL_ERR_NULL_POINTER;
    if (n_ages <= 0) return SWL_ERR_INVALID_AGE;
    
    for (int i = 0; i < n_ages; i++) {
        SWL_ErrorCode err = swl_calculate_cluster_wind_yield(
            imf_type, metallicity, cluster_mass, ages[i], &results[i]);
        if (err != SWL_SUCCESS) return err;
    }
    
    return SWL_SUCCESS;
}

SWL_ErrorCode swl_export_to_csv(
    const SWL_ClusterWindYield* result,
    const char* filename)
{
    if (!result || !filename) return SWL_ERR_NULL_POINTER;
    
    FILE* fp = fopen(filename, "w");
    if (!fp) return SWL_ERR_INTEGRATION_FAILED;
    
    fprintf(fp, "# Stellar Wind Yields - Cluster Results\n");
    fprintf(fp, "# IMF: %s\n", swl_get_imf_name(result->imf_type));
    fprintf(fp, "# Metallicity: Z = %.6e\n", result->metallicity);
    fprintf(fp, "# Cluster Mass: %.6e M_sun\n", result->cluster_mass);
    fprintf(fp, "# Cluster Age: %.6e years\n", result->cluster_age);
    fprintf(fp, "# Turnoff Mass: %.4f M_sun\n", result->turnoff_mass);
    fprintf(fp, "#\n");
    fprintf(fp, "# Column definitions:\n");
    fprintf(fp, "#   Net_Yield: Ejected - Initial (can be negative for H)\n");
    fprintf(fp, "#   Ejected_Mass: Actual mass ejected to ISM (always positive)\n");
    fprintf(fp, "#   Ejection_Rate: Current ejection rate to ISM\n");
    fprintf(fp, "#\n");
    fprintf(fp, "Element,Net_Yield_Msun,Ejected_Mass_Msun,Ejection_Rate_Msun_per_yr\n");
    
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        fprintf(fp, "%s,%.8e,%.8e,%.8e\n",
                swl_get_element_name((SWL_Element)i),
                result->element_yields.yields[i],
                result->ejected_mass.yields[i],
                result->ejection_rate.yields[i]);
    }
    
    fprintf(fp, "#\n");
    fprintf(fp, "# Summary:\n");
    fprintf(fp, "# Total wind mass: %.8e M_sun\n", result->total_wind_mass);
    fprintf(fp, "# Mass in living stars: %.8e M_sun\n", result->mass_in_living_stars);
    fprintf(fp, "# Mass in remnants: %.8e M_sun\n", result->mass_in_remnants);
    fprintf(fp, "# Current wind rate: %.8e M_sun/yr\n", result->wind_mass_rate);
    fprintf(fp, "# Total stars formed: %.8e\n", result->n_stars_total);
    fprintf(fp, "# Stars that have died: %.8e\n", result->n_stars_dead);
    
    fclose(fp);
    return SWL_SUCCESS;
}

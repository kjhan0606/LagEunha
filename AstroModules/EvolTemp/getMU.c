/**
 * =============================================================================
 * Mean Molecular Weight Calculator for Interstellar Gas
 * =============================================================================
 * 
 * This code calculates the mean molecular weight (mu) of interstellar gas
 * by solving the coupled Saha ionization equations for H and He simultaneously.
 * 
 * Physical Model:
 * ---------------
 * We consider a gas composed of hydrogen (H) and helium (He) with three
 * ionization processes:
 * 
 *   H   + e  <-->  H+   + 2e     (ionization energy chi1 = 13.6 eV)
 *   He  + e  <-->  He+  + 2e     (ionization energy chi2 = 24.6 eV)
 *   He+ + e  <-->  He++ + 2e     (ionization energy chi3 = 54.4 eV)
 * 
 * The Saha equations are:
 * 
 *   n_e * n_H+  / n_H0  = f1(T)    ... (1)
 *   n_e * n_He+ / n_He0 = f2(T)    ... (2)
 *   n_e * n_He++/ n_He+ = f3(T)    ... (3)
 * 
 * Combined with conservation equations:
 * 
 *   n_e = n_H+ + n_He+ + 2*n_He++  ... (charge neutrality)
 *   n_H = n_H0 + n_H+              ... (H conservation)
 *   n_He = n_He0 + n_He+ + n_He++  ... (He conservation)
 * 
 * These 6 equations reduce to a 4th-order polynomial in x = n_e / (n_H + n_He).
 * 
 * Metallicity Treatment:
 * ----------------------
 * The metallicity Z_metal (in solar units) affects the H and He mass fractions:
 *   X = X_primordial * (1 - Z)     where Z = Z_metal * Z_solar
 *   Y = Y_primordial * (1 - Z)
 * 
 * Numerical Approach:
 * -------------------
 * The Saha coefficient involves (2*pi*m_e*k_B/h^2)^(3/2) * T^(3/2).
 * Direct computation of h^2 ~ 4.4e-53 causes numerical underflow.
 * 
 * Solution: Factor as (m_e/h) * (k_B/h) to keep intermediate values manageable:
 *   m_e/h ~ 1.37e-1
 *   k_B/h ~ 2.08e10
 * 
 * Grid Parameters (for Galaxy Formation Simulations):
 * ---------------------------------------------------
 * Density range:     10^-31 to 10^-18 g/cm^3
 *   - IGM (Intergalactic Medium):    ~10^-31 to 10^-28 g/cm^3
 *   - CGM (Circumgalactic Medium):   ~10^-28 to 10^-26 g/cm^3
 *   - ISM (Interstellar Medium):     ~10^-26 to 10^-24 g/cm^3
 *   - Molecular Clouds:              ~10^-24 to 10^-21 g/cm^3
 *   - Star-forming Regions:          ~10^-21 to 10^-18 g/cm^3
 * 
 * Temperature range: 10^3 to 10^4 K (where ionization transitions occur)
 * Metallicity range: 0 to 3 (solar units)
 * 
 * =============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846L
#endif

/* ============================================================================
 *                         PHYSICAL CONSTANTS (CGS)
 * ============================================================================
 * Using CODATA recommended values where applicable.
 */

/* Fundamental constants */
#define C_LIGHT       2.99792458e10    /* Speed of light [cm/s] */
#define K_BOLTZMANN   1.380649e-16     /* Boltzmann constant [erg/K] */
#define H_PLANCK      6.62607015e-27   /* Planck constant [erg*s] */
#define M_ELECTRON    9.1093837e-28    /* Electron mass [g] */
#define M_HYDROGEN    1.6735575e-24    /* Hydrogen atom mass [g] (= m_proton approx) */
#define EV_TO_ERG     1.602176634e-12  /* 1 eV in erg */

/* Helium to Hydrogen mass ratio */
#define HE_TO_H_MASS_RATIO  3.9715     /* m_He / m_H (= ~4 * (1 - binding energy)) */

/* Statistical weights for ionization states */
#define G_H0          2                /* H ground state: 2 (spin degeneracy) */
#define G_Hp          1                /* H+ (proton): 1 */
#define G_He0         1                /* He ground state: 1 (singlet) */
#define G_Hep         2                /* He+: 2 (like H) */
#define G_Hepp        1                /* He++: 1 (bare nucleus) */

/**
 * Saha equation statistical weight ratios:
 * For reaction A -> B + e:  g_ratio = (g_B * g_e) / g_A, where g_e = 2
 */
#define SAHA_G1       (G_Hp * 2 / G_H0)     /* = 1 for H -> H+ */
#define SAHA_G2       (G_Hep * 2 / G_He0)   /* = 4 for He -> He+ */
#define SAHA_G3       (G_Hepp * 2 / G_Hep)  /* = 1 for He+ -> He++ */

/* Ionization energies [erg] */
#define CHI_H_I       (13.598 * EV_TO_ERG)  /* H ionization: 13.598 eV */
#define CHI_HE_I      (24.587 * EV_TO_ERG)  /* He I ionization: 24.587 eV */
#define CHI_HE_II     (54.418 * EV_TO_ERG)  /* He II ionization: 54.418 eV */

/* Primordial abundances by mass (Big Bang Nucleosynthesis values) */
#define Y_PRIMORDIAL  0.24             /* Primordial Helium mass fraction */
#define X_PRIMORDIAL  0.76             /* Primordial Hydrogen mass fraction (1 - Y) */
#define Z_SOLAR       0.02             /* Solar metallicity (mass fraction) */

/* ============================================================================
 *                 GRID PARAMETERS (for Galaxy Formation Simulations)
 * ============================================================================
 * 
 * Density range covers typical gas phases in galaxy simulations:
 *   - IGM:  rho ~ 10^-31 to 10^-28 g/cm^3  (cosmic web, voids)
 *   - CGM:  rho ~ 10^-28 to 10^-26 g/cm^3  (halo gas)
 *   - ISM:  rho ~ 10^-26 to 10^-24 g/cm^3  (diffuse interstellar)
 *   - GMC:  rho ~ 10^-24 to 10^-21 g/cm^3  (molecular clouds)
 *   - SF:   rho ~ 10^-21 to 10^-18 g/cm^3  (star-forming cores)
 */

#define N_DENSITY     128              /* Number of density grid points */
#define N_TEMPERATURE 128              /* Number of temperature grid points */
#define N_METALLICITY 32               /* Number of metallicity grid points */

/* Density range: typical for galaxy formation simulations */
#define LOG_RHO_MIN   -31.0            /* log10(rho_min) [g/cm^3]: IGM */
#define LOG_RHO_MAX   -18.0            /* log10(rho_max) [g/cm^3]: star-forming */

/* Temperature range: ionization transition zone */
#define TEMP_MIN       1000.0          /* Minimum temperature [K] */
#define TEMP_MAX       10000.0         /* Maximum temperature [K] */

/* Metallicity range: from primordial to super-solar */
#define Z_MIN          0.0             /* Minimum metallicity [Z_solar] */
#define Z_MAX          3.0             /* Maximum metallicity [Z_solar] */

/* Mean molecular weight bounds (for T/mu grid) */
#define MU_MIN         0.5             /* Minimum mu (fully ionized) */
#define MU_MAX         1.3             /* Maximum mu (neutral) */

/* Bisection solver parameters */
#define BISECTION_MAX_ITER  1000       /* Maximum iterations for root finding */
#define BISECTION_TOLERANCE 1.0e-10    /* Convergence tolerance */

/* helium increase with metallicity dY/dZ (typically it is 1.5 from observation */
#define DELTA_Y_OVER_DELTA_Z  1.5

/* ============================================================================
 *                            DATA STRUCTURES
 * ============================================================================ */

/**
 * Mean molecular weight lookup table structure
 * 
 * Three tables are stored:
 *   aMU[T][rho][Z]    : mu as function of (T, rho, Z)
 *   aTMU[T][rho][Z]   : T/mu as function of (T, rho, Z)  
 *   aT[Tmu][rho][Z]   : T as function of (T/mu, rho, Z) - inverse lookup
 */
typedef struct {
    double *log_rho;           /* Log10 of density grid points */
    double *log_T;             /* Log10 of temperature grid points */
    double *log_Tmu;           /* Log10 of T/mu grid points */
    double *Z;                 /* Metallicity grid points [Z_solar] */
    
    double ***mu_table;        /* mu(T, rho, Z) lookup table */
    double ***Tmu_table;       /* T/mu(T, rho, Z) lookup table */
    double ***T_from_Tmu;      /* T(T/mu, rho, Z) inverse lookup */
    
    int n_rho;                 /* Number of density points */
    int n_T;                   /* Number of temperature points */
    int n_Tmu;                 /* Number of T/mu points */
    int n_Z;                   /* Number of metallicity points */
    
    double log_T_step;         /* Step size in log10(T) */
    double log_Tmu_step;       /* Step size in log10(T/mu) */
    double log_Tmu_min;        /* Minimum log10(T/mu) */
    double log_Tmu_max;        /* Maximum log10(T/mu) */
    
    int initialized;           /* Initialization flag */
} MuGrid;

/* Global grid instance */
static MuGrid g_grid = {0};

/* ============================================================================
 *                    FOURTH-ORDER POLYNOMIAL SOLVER
 * ============================================================================
 * 
 * The coupled Saha equations reduce to a 4th-order polynomial:
 *   a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
 * 
 * where x = n_e / n_total is the electron fraction.
 * 
 * We solve this using the bisection method, which is robust and guaranteed
 * to converge for a bracketed root.
 */

/**
 * Evaluate the 4th-order polynomial: f(x) = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
 * 
 * Uses Horner's method for numerical stability:
 *   f(x) = a0 + x*(a1 + x*(a2 + x*(a3 + x*a4)))
 * 
 * @param a0, a1, a2, a3, a4  Polynomial coefficients
 * @param x                   Point at which to evaluate
 * @return                    Value of polynomial at x
 */
static double evaluate_quartic(double a0, double a1, double a2, 
                               double a3, double a4, double x) {
    return a0 + x * (a1 + x * (a2 + x * (a3 + x * a4)));
}

/**
 * Find root of 4th-order polynomial using bisection method
 * 
 * The physical constraint is 0 < x < 2 (electron fraction cannot exceed
 * maximum possible ionization: all H+ and all He++).
 * 
 * Algorithm:
 * 1. Scan from x=2 downward to bracket the root (find sign change)
 * 2. Use bisection to refine the root
 * 
 * @param a0, a1, a2, a3, a4  Polynomial coefficients (a4 = 1 typically)
 * @return                    Root x (electron fraction)
 */
static double solve_quartic_bisection(double a0, double a1, double a2, 
                                      double a3, double a4) {
    double x_lower, x_upper, x_mid;
    double f_lower, f_upper, f_mid;
    int i;
    
    /* -------------------------------------------------------------------------
     * Step 1: Bracket the root by scanning from high to low x
     * 
     * We use logarithmic spacing to efficiently search the range [1e-20, 2]
     * ------------------------------------------------------------------------- */
    const double x_max = 2.0;
    const double log_x_max = log10(x_max);
    const double log_x_min = -20.0;  /* x_min = 1e-20 */
    const double log_step = (log_x_max - log_x_min) / BISECTION_MAX_ITER;
    
    x_upper = x_max;
    f_upper = evaluate_quartic(a0, a1, a2, a3, a4, x_upper);
    
    /* Scan downward to find sign change */
    for (i = BISECTION_MAX_ITER - 1; i >= 0; i--) {
        x_lower = pow(10.0, log_step * i + log_x_min);
        f_lower = evaluate_quartic(a0, a1, a2, a3, a4, x_lower);
        
        /* Check for sign change: f_lower * f_upper < 0 */
        if (f_lower * f_upper <= 0.0) {
            break;  /* Root is bracketed between x_lower and x_upper */
        }
        
        /* Move upper bound down */
        x_upper = x_lower;
        f_upper = f_lower;
    }
    
    /* If no bracket found, return 0 (fully neutral gas) */
    if (i <= 0) {
        return 0.0;
    }
    
    /* -------------------------------------------------------------------------
     * Step 2: Bisection refinement
     * 
     * Repeatedly halve the interval until convergence
     * ------------------------------------------------------------------------- */
    double log_x_lower = log10(x_lower);
    double log_x_upper = log10(x_upper);
    
    while (fabs(log_x_upper - log_x_lower) > BISECTION_TOLERANCE) {
        /* Midpoint in log space for better behavior at small x */
        x_mid = pow(10.0, 0.5 * (log_x_upper + log_x_lower));
        f_mid = evaluate_quartic(a0, a1, a2, a3, a4, x_mid);
        
        if (f_mid * f_upper > 0.0) {
            /* Root is in lower half */
            log_x_upper = log10(x_mid);
            f_upper = f_mid;
        } else if (f_mid * f_upper < 0.0) {
            /* Root is in upper half */
            log_x_lower = log10(x_mid);
            f_lower = f_mid;
        } else {
            /* Exact root found */
            return x_mid;
        }
    }
    
    /* Return midpoint of final interval */
    return pow(10.0, 0.5 * (log_x_upper + log_x_lower));
}

/* ============================================================================
 *                    SAHA EQUATION COEFFICIENT CALCULATOR
 * ============================================================================
 * 
 * The Saha equation for ionization A -> B + e is:
 * 
 *   n_B * n_e     g_B * g_e   (2*pi*m_e*k_B*T)^(3/2)
 *   ---------  =  --------- * ----------------------- * exp(-chi / k_B*T)
 *     n_A           g_A              h^3
 * 
 * We define the Saha coefficient f(T) as the RHS.
 * 
 * Numerical trick: Instead of computing h^3 ~ 2.9e-79, we factor as:
 *   
 *   (2*pi*m_e*k_B/h^2)^(3/2) = [2*pi * (m_e/h) * (k_B/h)]^(3/2)
 * 
 * where m_e/h ~ 1.37e-1 and k_B/h ~ 2.08e10 are manageable numbers.
 */

/**
 * Compute Saha coefficient f(T) for a given ionization reaction
 * 
 * f(T) = g_ratio * (CR * T)^(3/2) * exp(-chi / k_B*T)
 * 
 * where CR = 2*pi * (m_e/h) * (k_B/h) is the numerically stable Saha prefactor.
 * 
 * @param T          Temperature [K]
 * @param chi        Ionization energy [erg]
 * @param g_ratio    Statistical weight ratio (g_ion * g_e) / g_neutral
 * @return           Saha coefficient [cm^-3]
 */
static double compute_saha_coefficient(double T, double chi, int g_ratio) {
    /*
     * CR = 2*pi * (m_e/h_P) * (k_B/h_P)
     *    = 2*pi * (9.109e-28 / 6.626e-27) * (1.381e-16 / 6.626e-27)
     *    = 2*pi * 0.1375 * 2.084e10
     *    = 1.800e10  [cm^-2 K^-1]
     * 
     * This avoids computing h^2 ~ 4.4e-53 which causes numerical issues.
     */
    const double CR = 2.0 * M_PI * (M_ELECTRON / H_PLANCK) * (K_BOLTZMANN / H_PLANCK);
    
    /* f = g_ratio * (CR*T)^(3/2) * exp(-chi/kT) */
    double CR_T = CR * T;                           /* [cm^-2] */
    double CR_T_32 = CR_T * sqrt(CR_T);             /* (CR*T)^(3/2) [cm^-3] */
    double boltzmann = exp(-chi / (K_BOLTZMANN * T));
    
    return g_ratio * CR_T_32 * boltzmann;
}

/* ============================================================================
 *                    MEAN MOLECULAR WEIGHT CALCULATOR
 * ============================================================================ */

/**
 * Calculate mean molecular weight by solving coupled Saha equations
 * 
 * This function solves for the ionization state of a H+He gas and returns
 * the mean molecular weight mu, defined as:
 * 
 *   mu = rho / (n_total * m_H)
 * 
 * where n_total = n_H0 + n_H+ + n_He0 + n_He+ + n_He++ + n_e
 * 
 * The coupled Saha equations are reduced to a 4th-order polynomial in the
 * electron fraction x = n_e / (n_H + n_He), then solved by bisection.
 * 
 * Metallicity effect:
 *   Higher metallicity means more metals, which reduce the H and He fractions.
 *   X = X_primordial * (1 - Z_mass)
 *   Y = Y_primordial * (1 - Z_mass)
 *   where Z_mass = Z_metal * Z_solar
 * 
 * @param rho       Physical density [g/cm^3]
 * @param T         Temperature [K]
 * @param Z_metal   Metallicity in solar units (Z_solar = 0.02)
 * @return          Mean molecular weight [m_H units]
 */
double calculate_mu(double rho, double T, double Z_metal) {
    /* -------------------------------------------------------------------------
     * Step 1: Convert metallicity to mass fractions
     * 
     * Z_mass = Z_metal * Z_solar  (metal mass fraction)
     * X = X_primordial * (1 - Z_mass)  (hydrogen mass fraction)
     * Y = Y_primordial * (1 - Z_mass)  (helium mass fraction)
     * 
     * This ensures X + Y + Z = 1 (approximately, for Z << 1)
     * ------------------------------------------------------------------------- */
	double Z_mass = Z_metal * Z_SOLAR; 
	if (Z_mass > 0.5) Z_mass = 0.5; 
	if (Z_mass < 0.0) Z_mass = 0.0;

	double Y_He = Y_PRIMORDIAL + DELTA_Y_OVER_DELTA_Z * Z_mass;
	double X_H = 1.0 - Y_He - Z_mass;

	/* to prevent the overshoot */
	if (X_H < 0.0) { 
		X_H = 0.0; 
		Y_He = 1.0 - Z_mass; 
	}

    
    /* -------------------------------------------------------------------------
     * Step 2: Compute number densities from mass density
     * 
     * Given mass fractions X (hydrogen) and Y (helium):
     *   n_H  = X * rho / m_H
     *   n_He = Y * rho / (4 * m_H)  [He mass ~ 4 * m_H]
     * 
     * Number fraction of He relative to H:
     *   f_He = n_He / n_H = Y / (4 * X)
     * ------------------------------------------------------------------------- */
    double f_He = Y_He / (HE_TO_H_MASS_RATIO * X_H);  /* n_He / n_H */
    
    /* Total hydrogen number density */
    double n_H = X_H * rho / M_HYDROGEN;              /* [cm^-3] */
    
    /* Total helium number density */
    double n_He = n_H * f_He;                         /* [cm^-3] */
    
    /* Total heavy particle density (excluding electrons) */
    double n_total_nuclei = n_H + n_He;               /* [cm^-3] */
    
    /* -------------------------------------------------------------------------
     * Step 3: Handle temperature limits
     * 
     * At very low T: all neutral -> mu = (4*f_He + 1) / (1 + f_He)
     * At very high T: all ionized -> mu = (4*f_He + 1) / (2 + 3*f_He)
     * ------------------------------------------------------------------------- */
    if (T < TEMP_MIN) {
        /* Fully neutral: n_total = n_H + n_He (no free electrons) */
        return (HE_TO_H_MASS_RATIO * f_He + 1.0) / (1.0 + f_He);
    }
    if (T > TEMP_MAX) {
        /* Fully ionized: n_e = n_H + 2*n_He (H+ and He++) */
        return (HE_TO_H_MASS_RATIO * f_He + 1.0) / (2.0 + 3.0 * f_He);
    }
    
    /* -------------------------------------------------------------------------
     * Step 4: Compute Saha coefficients
     * 
     * f1 = Saha coeff for H -> H+    (chi = 13.6 eV, g_ratio = 1)
     * f2 = Saha coeff for He -> He+  (chi = 24.6 eV, g_ratio = 4)
     * f3 = Saha coeff for He+ -> He++ (chi = 54.4 eV, g_ratio = 1)
     * ------------------------------------------------------------------------- */
    double f1 = compute_saha_coefficient(T, CHI_H_I,  SAHA_G1);   /* H ionization */
    double f2 = compute_saha_coefficient(T, CHI_HE_I, SAHA_G2);   /* He I ionization */
    double f3 = compute_saha_coefficient(T, CHI_HE_II, SAHA_G3);  /* He II ionization */
    
    /* -------------------------------------------------------------------------
     * Step 5: Normalize Saha coefficients
     * 
     * Define dimensionless parameters:
     *   p1 = f1 / n_total = Saha coeff relative to total density
     *   p2 = f2 / n_total
     *   p3 = f3 / n_total
     *   pnH  = n_H / n_total  = 1 / (1 + f_He)
     *   pnHe = n_He / n_total = f_He / (1 + f_He)
     * ------------------------------------------------------------------------- */
    double p1 = f1 / n_total_nuclei;
    double p2 = f2 / n_total_nuclei;
    double p3 = f3 / n_total_nuclei;
    double pnH  = n_H / n_total_nuclei;    /* = 1 / (1 + f_He) */
    double pnHe = n_He / n_total_nuclei;   /* = f_He / (1 + f_He) */
    
    /* -------------------------------------------------------------------------
     * Step 6: Construct 4th-order polynomial coefficients
     * 
     * The variable x = n_e / n_total represents the electron fraction.
     * 
     * From the Saha equations and charge neutrality:
     *   n_H+ = n_H0 * f1 / n_e  -->  n_H+ = n_H * f1 / (n_e + f1)
     *   n_He+ = complicated expression involving f2, f3, n_e
     *   n_He++ = n_He+ * f3 / n_e
     *   n_e = n_H+ + n_He+ + 2*n_He++  (charge neutrality)
     * 
     * After algebra, we get a 4th-order polynomial in x:
     *   a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
     * 
     * The coefficients are derived from eliminating n_H0, n_H+, n_He0, n_He+, 
     * n_He++ in favor of x.
     * ------------------------------------------------------------------------- */
    double a4 = 1.0;
    double a3 = p1 + p2;
    double a2 = p1 * p2 + p2 * p3 - p2 * pnHe - p1 * pnH;
    double a1 = p1 * p2 * p3 - p1 * p2 * (pnH + pnHe) - 2.0 * p2 * p3 * pnHe;
    double a0 = -p1 * p2 * p3 * (pnH + 2.0 * pnHe);
    
    /* -------------------------------------------------------------------------
     * Step 7: Solve for electron fraction x
     * ------------------------------------------------------------------------- */
    double x = solve_quartic_bisection(a0, a1, a2, a3, a4);
    
    /* Electron density */
    double n_e = x * n_total_nuclei;
    
    /* -------------------------------------------------------------------------
     * Step 8: Compute individual species densities
     * 
     * From Saha equations:
     *   n_H+ / n_H0 = f1 / n_e  -->  n_H0 = n_H / (1 + f1/n_e)
     *   n_He+ / n_He0 = f2 / n_e
     *   n_He++ / n_He+ = f3 / n_e
     * ------------------------------------------------------------------------- */
    
    /* Avoid division by zero */
    if (n_e < 1.0e-30) n_e = 1.0e-30;
    
    /* Hydrogen species */
    double n_H0 = n_H / (1.0 + f1 / n_e);
    double n_Hp = n_H - n_H0;    /* = n_H0 * f1 / n_e */
    
    /* Helium species 
     * n_He = n_He0 + n_He+ + n_He++
     * n_He+/n_He0 = f2/n_e
     * n_He++/n_He+ = f3/n_e
     * --> n_He0 = n_He / (1 + f2/n_e + f2*f3/n_e^2)
     */
    double ratio2 = f2 / n_e;
    double ratio3 = f3 / n_e;
    double n_He0 = n_He / (1.0 + ratio2 + ratio2 * ratio3);
    double n_Hep = n_He0 * ratio2;
    double n_Hepp = n_Hep * ratio3;
    
    /* -------------------------------------------------------------------------
     * Step 9: Compute ionization fractions and mean molecular weight
     * 
     * Ionization fractions:
     *   x_H   = n_H+ / n_H
     *   x_He  = (n_He+ + n_He++) / n_He  (singly or doubly ionized)
     *   x_He2 = n_He++ / n_He  (doubly ionized only)
     * 
     * Mean molecular weight formula:
     *   mu = (4*f_He + 1) / (1 + x_H + f_He*(1 + x_He + x_He2))
     * 
     * This accounts for:
     * - Numerator: total mass in units of m_H
     * - Denominator: total number of free particles per H nucleus
     * ------------------------------------------------------------------------- */
    double x_H = n_Hp / n_H;
    double x_He = (n_Hep + n_Hepp) / n_He;
    double x_He2 = n_Hepp / n_He;
    
    double mu = (HE_TO_H_MASS_RATIO * f_He + 1.0) / 
                (1.0 + x_H + f_He * (1.0 + x_He + x_He2));
    
    return mu;
}

/* ============================================================================
 *                         3D ARRAY MEMORY MANAGEMENT
 * ============================================================================ */

/**
 * Allocate a 3D array of doubles with dimensions [nx][ny][nz]
 */
static double*** allocate_3d_array(int nx, int ny, int nz) {
    double ***arr = (double ***)malloc(nx * sizeof(double **));
    if (arr == NULL) return NULL;
    
    for (int i = 0; i < nx; i++) {
        arr[i] = (double **)malloc(ny * sizeof(double *));
        if (arr[i] == NULL) return NULL;
        
        for (int j = 0; j < ny; j++) {
            arr[i][j] = (double *)malloc(nz * sizeof(double));
            if (arr[i][j] == NULL) return NULL;
        }
    }
    return arr;
}

/**
 * Free a 3D array allocated by allocate_3d_array
 */
static void free_3d_array(double ***arr, int nx, int ny) {
    if (arr == NULL) return;
    
    for (int i = 0; i < nx; i++) {
        if (arr[i] != NULL) {
            for (int j = 0; j < ny; j++) {
                if (arr[i][j] != NULL) {
                    free(arr[i][j]);
                }
            }
            free(arr[i]);
        }
    }
    free(arr);
}

/* ============================================================================
 *                           GRID INITIALIZATION
 * ============================================================================ */

/**
 * Initialize the mean molecular weight lookup tables
 * 
 * This function pre-computes mu values on a 3D grid of (density, temperature,
 * metallicity) and also builds an inverse lookup table T(T/mu, rho, Z).
 * 
 * The inverse table is useful when energy is known but temperature is not,
 * since internal energy ~ k_B * T / mu.
 * 
 * @return 0 on success, -1 on failure
 */
int init_mu_grid(void) {
    if (g_grid.initialized) {
        printf("Grid already initialized.\n");
        return 0;
    }
    
    printf("Initializing mean molecular weight grid...\n");
    printf("  Density range: 10^%.0f to 10^%.0f g/cm^3 (galaxy simulations)\n",
           LOG_RHO_MIN, LOG_RHO_MAX);
    
    /* -------------------------------------------------------------------------
     * Allocate grid arrays
     * ------------------------------------------------------------------------- */
    g_grid.n_rho = N_DENSITY;
    g_grid.n_T = N_TEMPERATURE;
    g_grid.n_Tmu = N_TEMPERATURE;  /* Same size for inverse table */
    g_grid.n_Z = N_METALLICITY;
    
    g_grid.log_rho = (double *)malloc(N_DENSITY * sizeof(double));
    g_grid.log_T = (double *)malloc(N_TEMPERATURE * sizeof(double));
    g_grid.log_Tmu = (double *)malloc(N_TEMPERATURE * sizeof(double));
    g_grid.Z = (double *)malloc(N_METALLICITY * sizeof(double));
    
    g_grid.mu_table = allocate_3d_array(N_TEMPERATURE, N_DENSITY, N_METALLICITY);
    g_grid.Tmu_table = allocate_3d_array(N_TEMPERATURE, N_DENSITY, N_METALLICITY);
    g_grid.T_from_Tmu = allocate_3d_array(N_TEMPERATURE, N_DENSITY, N_METALLICITY);
    
    if (g_grid.log_rho == NULL || g_grid.log_T == NULL || 
        g_grid.log_Tmu == NULL || g_grid.Z == NULL ||
        g_grid.mu_table == NULL || g_grid.Tmu_table == NULL ||
        g_grid.T_from_Tmu == NULL) {
        fprintf(stderr, "Error: Failed to allocate grid memory.\n");
        return -1;
    }
    
    /* -------------------------------------------------------------------------
     * Set up grid spacing
     * ------------------------------------------------------------------------- */
    double log_T_min = log10(TEMP_MIN);
    double log_T_max = log10(TEMP_MAX);
    g_grid.log_T_step = (log_T_max - log_T_min) / (N_TEMPERATURE - 1);
    
    /* T/mu range: from T_min/mu_max to T_max/mu_min */
    g_grid.log_Tmu_min = log10(TEMP_MIN / MU_MAX);
    g_grid.log_Tmu_max = log10(TEMP_MAX / MU_MIN);
    g_grid.log_Tmu_step = (g_grid.log_Tmu_max - g_grid.log_Tmu_min) / (g_grid.n_Tmu - 1);
    
    double log_rho_step = (LOG_RHO_MAX - LOG_RHO_MIN) / (N_DENSITY - 1);
    double Z_step = (Z_MAX - Z_MIN) / (N_METALLICITY - 1);
    
    /* Fill coordinate arrays */
    for (int i = 0; i < N_DENSITY; i++) {
        g_grid.log_rho[i] = LOG_RHO_MIN + i * log_rho_step;
    }
    for (int i = 0; i < N_TEMPERATURE; i++) {
        g_grid.log_T[i] = log_T_min + i * g_grid.log_T_step;
    }
    for (int i = 0; i < g_grid.n_Tmu; i++) {
        g_grid.log_Tmu[i] = g_grid.log_Tmu_min + i * g_grid.log_Tmu_step;
    }
    for (int k = 0; k < N_METALLICITY; k++) {
        g_grid.Z[k] = Z_MIN + k * Z_step;
    }
    
    /* -------------------------------------------------------------------------
     * Compute mu and T/mu tables
     * ------------------------------------------------------------------------- */
    printf("Computing mu table (%d x %d x %d = %d points)...\n",
           N_TEMPERATURE, N_DENSITY, N_METALLICITY,
           N_TEMPERATURE * N_DENSITY * N_METALLICITY);
    
    for (int j = 0; j < N_DENSITY; j++) {
        double rho = pow(10.0, g_grid.log_rho[j]);
        
        for (int i = 0; i < N_TEMPERATURE; i++) {
            double T = pow(10.0, g_grid.log_T[i]);
            
            for (int k = 0; k < N_METALLICITY; k++) {
                double Z_metal = g_grid.Z[k];
                
                double mu = calculate_mu(rho, T, Z_metal);
                g_grid.mu_table[i][j][k] = mu;
                g_grid.Tmu_table[i][j][k] = T / mu;
            }
        }
        
        /* Progress indicator */
        if ((j + 1) % 20 == 0 || j == N_DENSITY - 1) {
            printf("  Progress: %d/%d density slices\n", j + 1, N_DENSITY);
        }
    }
    
    /* -------------------------------------------------------------------------
     * Build inverse table: T(T/mu, rho, Z)
     * 
     * For each T/mu value, find the corresponding T by interpolating
     * in the forward table.
     * ------------------------------------------------------------------------- */
    printf("Building inverse T(T/mu) table...\n");
    
    for (int j = 0; j < N_DENSITY; j++) {
        for (int k = 0; k < N_METALLICITY; k++) {
            for (int i = 0; i < g_grid.n_Tmu; i++) {
                double Tmu = pow(10.0, g_grid.log_Tmu[i]);
                
                /* Check bounds */
                if (Tmu < g_grid.Tmu_table[0][j][k]) {
                    /* Below minimum T/mu: extrapolate with mu_max */
                    g_grid.T_from_Tmu[i][j][k] = Tmu * MU_MAX;
                }
                else if (Tmu >= g_grid.Tmu_table[N_TEMPERATURE-1][j][k]) {
                    /* Above maximum T/mu: extrapolate with mu_min */
                    g_grid.T_from_Tmu[i][j][k] = Tmu * MU_MIN;
                }
                else {
                    /* Find bracketing interval and interpolate */
                    for (int m = 0; m < N_TEMPERATURE - 1; m++) {
                        if (g_grid.Tmu_table[m][j][k] <= Tmu && 
                            g_grid.Tmu_table[m+1][j][k] > Tmu) {
                            
                            double T_low = pow(10.0, g_grid.log_T[m]);
                            double T_high = pow(10.0, g_grid.log_T[m+1]);
                            double Tmu_low = g_grid.Tmu_table[m][j][k];
                            double Tmu_high = g_grid.Tmu_table[m+1][j][k];
                            
                            /* Linear interpolation */
                            g_grid.T_from_Tmu[i][j][k] = T_low + 
                                (T_high - T_low) / (Tmu_high - Tmu_low) * (Tmu - Tmu_low);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    g_grid.initialized = 1;
    printf("Grid initialization complete.\n");
    
    return 0;
}

/**
 * Free all memory allocated for the grid
 */
void free_mu_grid(void) {
    if (!g_grid.initialized) return;
    
    if (g_grid.log_rho) free(g_grid.log_rho);
    if (g_grid.log_T) free(g_grid.log_T);
    if (g_grid.log_Tmu) free(g_grid.log_Tmu);
    if (g_grid.Z) free(g_grid.Z);
    
    free_3d_array(g_grid.mu_table, N_TEMPERATURE, N_DENSITY);
    free_3d_array(g_grid.Tmu_table, N_TEMPERATURE, N_DENSITY);
    free_3d_array(g_grid.T_from_Tmu, g_grid.n_Tmu, N_DENSITY);
    
    memset(&g_grid, 0, sizeof(MuGrid));
    
    printf("Grid memory freed.\n");
}

/* ============================================================================
 *                         INTERPOLATION FUNCTIONS
 * ============================================================================ */

/**
 * Find index for interpolation using binary search
 * Returns i such that arr[i] <= value < arr[i+1]
 */
static int find_index(const double *arr, int n, double value) {
    if (value <= arr[0]) return 0;
    if (value >= arr[n-1]) return n - 2;
    
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (arr[mid] <= value) lo = mid;
        else hi = mid;
    }
    return lo;
}

/**
 * Get mean molecular weight from lookup table
 * 
 * @param rho     Physical density [g/cm^3]
 * @param T       Temperature [K]
 * @param Z_metal Metallicity [Z_solar]
 * @return        Mean molecular weight [m_H units]
 */
double get_mean_molecular_weight(double rho, double T, double Z_metal) {
    if (!g_grid.initialized) {
        fprintf(stderr, "Error: Grid not initialized. Call init_mu_grid() first.\n");
        return -1.0;
    }
    
    /* Handle out-of-range temperatures analytically */
    double Z_mass = Z_metal * Z_SOLAR;
    if (Z_mass > 0.5) Z_mass = 0.5;
    if (Z_mass < 0.0) Z_mass = 0.0;
    
	double Y_He = Y_PRIMORDIAL + DELTA_Y_OVER_DELTA_Z * Z_mass;
	double X_H = 1.0 - Y_He - Z_mass;

    double f_He = Y_He / (HE_TO_H_MASS_RATIO * X_H);

	/* to prevent the overshoot */
	if (X_H < 0.0) { 
		X_H = 0.0; 
		Y_He = 1.0 - Z_mass; 
	}
    
    if (T < TEMP_MIN) {
        /* Fully neutral */
        return (HE_TO_H_MASS_RATIO * f_He + 1.0) / (1.0 + f_He);
    }
    if (T > TEMP_MAX) {
        /* Fully ionized */
        return (HE_TO_H_MASS_RATIO * f_He + 1.0) / (2.0 + 3.0 * f_He);
    }
    
    /* Find grid indices */
    double log_rho = log10(rho);
    double log_T = log10(T);
    
    /* Clamp to grid bounds */
    if (log_rho < LOG_RHO_MIN) log_rho = LOG_RHO_MIN;
    if (log_rho > LOG_RHO_MAX) log_rho = LOG_RHO_MAX;
    if (Z_metal < Z_MIN) Z_metal = Z_MIN;
    if (Z_metal > Z_MAX) Z_metal = Z_MAX;
    
    int j = find_index(g_grid.log_rho, g_grid.n_rho, log_rho);
    int i = find_index(g_grid.log_T, g_grid.n_T, log_T);
    int k = find_index(g_grid.Z, g_grid.n_Z, Z_metal);
    
    /* Simple nearest-neighbor lookup (can be extended to trilinear) */
    return g_grid.mu_table[i][j][k];
}

/**
 * Get temperature from T/mu and density
 * 
 * This is useful when internal energy is known but temperature is not,
 * since u = (3/2) * k_B * T / (mu * m_H) --> T/mu = 2*u*m_H / (3*k_B)
 * 
 * @param Tmu     Temperature divided by mean molecular weight [K]
 * @param rho     Physical density [g/cm^3]
 * @param Z_metal Metallicity in solar units
 * @return        Temperature [K]
 */
double get_temperature_from_Tmu(double Tmu, double rho, double Z_metal) {
    if (!g_grid.initialized) {
        fprintf(stderr, "Error: Grid not initialized. Call init_mu_grid() first.\n");
        return -1.0;
    }
    
    double Tmu_min = pow(10.0, g_grid.log_Tmu_min);
    double Tmu_max = pow(10.0, g_grid.log_Tmu_max);
    
    /* Handle out-of-range T/mu */
    if (Tmu < Tmu_min) {
        return Tmu * MU_MAX;
    }
    if (Tmu > Tmu_max) {
        return Tmu * MU_MIN;
    }
    
    /* Find grid indices */
    double log_rho = log10(rho);
    double log_Tmu = log10(Tmu);
    
    if (log_rho < LOG_RHO_MIN) log_rho = LOG_RHO_MIN;
    if (log_rho > LOG_RHO_MAX) log_rho = LOG_RHO_MAX;
    if (Z_metal < Z_MIN) Z_metal = Z_MIN;
    if (Z_metal > Z_MAX) Z_metal = Z_MAX;
    
    int j = find_index(g_grid.log_rho, g_grid.n_rho, log_rho);
    int i = find_index(g_grid.log_Tmu, g_grid.n_Tmu, log_Tmu);
    int k = find_index(g_grid.Z, g_grid.n_Z, Z_metal);
    
    return g_grid.T_from_Tmu[i][j][k];
}

/* ============================================================================
 *                            UTILITY FUNCTIONS
 * ============================================================================ */

/**
 * Print grid information
 */
void print_grid_info(void) {
    if (!g_grid.initialized) {
        printf("Grid not initialized.\n");
        return;
    }
    
    printf("\n");
    printf("####################################################################\n");
    printf("#       Mean Molecular Weight Grid Information                     #\n");
    printf("#       (Optimized for Galaxy Formation Simulations)               #\n");
    printf("####################################################################\n");
    printf("#  Density points:     %-6d                                       #\n", g_grid.n_rho);
    printf("#  Temperature points: %-6d                                       #\n", g_grid.n_T);
    printf("#  Metallicity points: %-6d                                       #\n", g_grid.n_Z);
    printf("#  Total grid points:  %-6d                                       #\n",
           g_grid.n_rho * g_grid.n_T * g_grid.n_Z);
    printf("####################################################################\n");
    printf("#  Density range:      10^%.0f to 10^%.0f g/cm^3                     #\n",
           LOG_RHO_MIN, LOG_RHO_MAX);
    printf("#    - IGM:            10^-31 to 10^-28 g/cm^3                      #\n");
    printf("#    - CGM:            10^-28 to 10^-26 g/cm^3                      #\n");
    printf("#    - ISM:            10^-26 to 10^-24 g/cm^3                      #\n");
    printf("#    - GMC:            10^-24 to 10^-21 g/cm^3                      #\n");
    printf("#    - Star-forming:   10^-21 to 10^-18 g/cm^3                      #\n");
    printf("#  Temperature range:  %.0f to %.0f K                              #\n",
           TEMP_MIN, TEMP_MAX);
    printf("#  Metallicity range:  %.1f to %.1f Z_solar                         #\n",
           Z_MIN, Z_MAX);
    printf("####################################################################\n\n");
}

/**
 * Save grid to binary file
 */
int save_grid_to_file(const char *filename) {
    if (!g_grid.initialized) {
        fprintf(stderr, "Error: Grid not initialized.\n");
        return -1;
    }
    
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Error: Cannot open %s for writing.\n", filename);
        return -1;
    }
    
    /* Write dimensions */
    fwrite(&g_grid.n_rho, sizeof(int), 1, fp);
    fwrite(&g_grid.n_T, sizeof(int), 1, fp);
    fwrite(&g_grid.n_Tmu, sizeof(int), 1, fp);
    fwrite(&g_grid.n_Z, sizeof(int), 1, fp);
    
    /* Write grid parameters */
    fwrite(&g_grid.log_T_step, sizeof(double), 1, fp);
    fwrite(&g_grid.log_Tmu_step, sizeof(double), 1, fp);
    fwrite(&g_grid.log_Tmu_min, sizeof(double), 1, fp);
    fwrite(&g_grid.log_Tmu_max, sizeof(double), 1, fp);
    
    /* Write coordinate arrays */
    fwrite(g_grid.log_rho, sizeof(double), g_grid.n_rho, fp);
    fwrite(g_grid.log_T, sizeof(double), g_grid.n_T, fp);
    fwrite(g_grid.log_Tmu, sizeof(double), g_grid.n_Tmu, fp);
    fwrite(g_grid.Z, sizeof(double), g_grid.n_Z, fp);
    
    /* Write tables */
    for (int i = 0; i < g_grid.n_T; i++) {
        for (int j = 0; j < g_grid.n_rho; j++) {
            fwrite(g_grid.mu_table[i][j], sizeof(double), g_grid.n_Z, fp);
        }
    }
    for (int i = 0; i < g_grid.n_T; i++) {
        for (int j = 0; j < g_grid.n_rho; j++) {
            fwrite(g_grid.Tmu_table[i][j], sizeof(double), g_grid.n_Z, fp);
        }
    }
    for (int i = 0; i < g_grid.n_Tmu; i++) {
        for (int j = 0; j < g_grid.n_rho; j++) {
            fwrite(g_grid.T_from_Tmu[i][j], sizeof(double), g_grid.n_Z, fp);
        }
    }
    
    fclose(fp);
    printf("Grid saved to %s\n", filename);
    return 0;
}

/**
 * Load grid from binary file
 */
int load_grid_from_file(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error: Cannot open %s for reading.\n", filename);
        return -1;
    }
    
    free_mu_grid();
    
    /* Read dimensions */
    if (fread(&g_grid.n_rho, sizeof(int), 1, fp) != 1) goto read_error;
    if (fread(&g_grid.n_T, sizeof(int), 1, fp) != 1) goto read_error;
    if (fread(&g_grid.n_Tmu, sizeof(int), 1, fp) != 1) goto read_error;
    if (fread(&g_grid.n_Z, sizeof(int), 1, fp) != 1) goto read_error;
    
    /* Read grid parameters */
    if (fread(&g_grid.log_T_step, sizeof(double), 1, fp) != 1) goto read_error;
    if (fread(&g_grid.log_Tmu_step, sizeof(double), 1, fp) != 1) goto read_error;
    if (fread(&g_grid.log_Tmu_min, sizeof(double), 1, fp) != 1) goto read_error;
    if (fread(&g_grid.log_Tmu_max, sizeof(double), 1, fp) != 1) goto read_error;
    
    /* Allocate arrays */
    g_grid.log_rho = (double *)malloc(g_grid.n_rho * sizeof(double));
    g_grid.log_T = (double *)malloc(g_grid.n_T * sizeof(double));
    g_grid.log_Tmu = (double *)malloc(g_grid.n_Tmu * sizeof(double));
    g_grid.Z = (double *)malloc(g_grid.n_Z * sizeof(double));
    g_grid.mu_table = allocate_3d_array(g_grid.n_T, g_grid.n_rho, g_grid.n_Z);
    g_grid.Tmu_table = allocate_3d_array(g_grid.n_T, g_grid.n_rho, g_grid.n_Z);
    g_grid.T_from_Tmu = allocate_3d_array(g_grid.n_Tmu, g_grid.n_rho, g_grid.n_Z);
    
    /* Read coordinate arrays */
    if (fread(g_grid.log_rho, sizeof(double), g_grid.n_rho, fp) != (size_t)g_grid.n_rho) goto read_error;
    if (fread(g_grid.log_T, sizeof(double), g_grid.n_T, fp) != (size_t)g_grid.n_T) goto read_error;
    if (fread(g_grid.log_Tmu, sizeof(double), g_grid.n_Tmu, fp) != (size_t)g_grid.n_Tmu) goto read_error;
    if (fread(g_grid.Z, sizeof(double), g_grid.n_Z, fp) != (size_t)g_grid.n_Z) goto read_error;
    
    /* Read tables */
    for (int i = 0; i < g_grid.n_T; i++) {
        for (int j = 0; j < g_grid.n_rho; j++) {
            if (fread(g_grid.mu_table[i][j], sizeof(double), g_grid.n_Z, fp) != (size_t)g_grid.n_Z) goto read_error;
        }
    }
    for (int i = 0; i < g_grid.n_T; i++) {
        for (int j = 0; j < g_grid.n_rho; j++) {
            if (fread(g_grid.Tmu_table[i][j], sizeof(double), g_grid.n_Z, fp) != (size_t)g_grid.n_Z) goto read_error;
        }
    }
    for (int i = 0; i < g_grid.n_Tmu; i++) {
        for (int j = 0; j < g_grid.n_rho; j++) {
            if (fread(g_grid.T_from_Tmu[i][j], sizeof(double), g_grid.n_Z, fp) != (size_t)g_grid.n_Z) goto read_error;
        }
    }
    
    fclose(fp);
    g_grid.initialized = 1;
    printf("Grid loaded from %s\n", filename);
    return 0;

read_error:
    fprintf(stderr, "Error: Failed to read grid file.\n");
    fclose(fp);
    free_mu_grid();
    return -1;
}

/* ============================================================================
 *                              TEST PROGRAM
 * ============================================================================ */

#ifdef TEST_GETMU

int main(int argc, char *argv[]) {
    printf("####################################################################\n");
    printf("#       Mean Molecular Weight Calculator - Test Program            #\n");
    printf("#       (Coupled Saha Equation with 4th-Order Polynomial)          #\n");
    printf("#       Optimized for Galaxy Formation Simulations                 #\n");
    printf("####################################################################\n\n");
    
    /* Print physical constants used */
    printf("Physical Constants (CGS):\n");
    printf("  k_B       = %.6e erg/K\n", K_BOLTZMANN);
    printf("  h         = %.6e erg*s\n", H_PLANCK);
    printf("  m_e       = %.6e g\n", M_ELECTRON);
    printf("  m_H       = %.6e g\n", M_HYDROGEN);
    printf("  chi_H     = %.3f eV\n", CHI_H_I / EV_TO_ERG);
    printf("  chi_HeI   = %.3f eV\n", CHI_HE_I / EV_TO_ERG);
    printf("  chi_HeII  = %.3f eV\n", CHI_HE_II / EV_TO_ERG);
    printf("  X_primord = %.2f\n", X_PRIMORDIAL);
    printf("  Y_primord = %.2f\n", Y_PRIMORDIAL);
    printf("  Z_solar   = %.2f\n\n", Z_SOLAR);
    
    /* Initialize grid */
    if (init_mu_grid() != 0) {
        fprintf(stderr, "Failed to initialize grid.\n");
        return 1;
    }
    
    print_grid_info();
    
    /* Test cases */
    printf("###################################################################\n");
    printf("                          TEST CASES\n");
    printf("###################################################################\n\n");
    
    /* Test 1: IGM - low density, various metallicities */
    printf("Test 1: Intergalactic Medium (IGM)\n");
    printf("  rho = 1e-29 g/cm^3(cosmic web)\n");
    {
        double rho = 1.0e-29;
        double T_vals[] = {500, 5000, 15000};
        double Z_vals[] = {0.0, 0.1, 1.0};
        
        printf("      T [K]    Z [Z#]    mu\n");
        printf("    --------  --------  ------\n");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double mu = calculate_mu(rho, T_vals[i], Z_vals[j]);
                printf("    %7.0f    %4.1f     %.4f\n", T_vals[i], Z_vals[j], mu);
            }
        }
    }
    printf("\n");
    
    /* Test 2: ISM - typical interstellar medium */
    printf("Test 2: Interstellar Medium (ISM)\n");
    printf("  rho = 1e-24 g/cm^3(diffuse ISM)\n");
    {
        double rho = 1.0e-24;
        double T_vals[] = {500, 5000, 15000};
        double Z_vals[] = {0.0, 1.0, 2.0};
        
        printf("      T [K]    Z [Z#]    mu\n");
        printf("    --------  --------  ------\n");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double mu = calculate_mu(rho, T_vals[i], Z_vals[j]);
                printf("    %7.0f    %4.1f     %.4f\n", T_vals[i], Z_vals[j], mu);
            }
        }
    }
    printf("\n");
    
    /* Test 3: Star-forming region - high density */
    printf("Test 3: Star-forming Region\n");
    printf("  rho = 1e-20 g/cm^3(dense core)\n");
    {
        double rho = 1.0e-20;
        double T_vals[] = {500, 5000, 15000};
        double Z_vals[] = {0.5, 1.0, 2.0};
        
        printf("      T [K]    Z [Z#]    mu\n");
        printf("    --------  --------  ------\n");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double mu = calculate_mu(rho, T_vals[i], Z_vals[j]);
                printf("    %7.0f    %4.1f     %.4f\n", T_vals[i], Z_vals[j], mu);
            }
        }
    }
    printf("\n");
    
    /* Temperature dependence at fixed density */
    printf("###################################################################\n");
    printf("       Temperature Dependence (rho = 1e-24 g/cm^3 Z = 1.0 Z#)\n");
    printf("###################################################################\n\n");
    printf("      T [K]        mu       Comment\n");
    printf("    ----------   ------   -------------------------\n");
    
    double rho_test = 1.0e-24;
    double Z_test = 1.0;
    double T_tests[] = {100, 500, 1000, 2000, 4000, 6000, 8000, 10000, 15000, 50000};
    const char *comments[] = {
        "Molecular cloud (neutral)",
        "Cold neutral medium",
        "Start of ionization zone",
        "Partial H ionization",
        "H ~50% ionized",
        "H mostly ionized",
        "H fully ionized",
        "End of ionization zone",
        "Fully ionized",
        "Hot plasma"
    };
    int n_tests = sizeof(T_tests) / sizeof(T_tests[0]);
    
    for (int i = 0; i < n_tests; i++) {
        double mu = calculate_mu(rho_test, T_tests[i], Z_test);
        printf("    %8.0f     %.4f   %s\n", T_tests[i], mu, comments[i]);
    }
    
    /* Density dependence at fixed temperature */
    printf("\n###################################################################\n");
    printf("         Density Dependence (T = 6000 K, Z = 1.0 Z#)\n");
    printf("###################################################################\n\n");
    printf("    rho [g/cm^3]        mu      Gas Phase\n");
    printf("    -------------    ------   -----------------\n");
    
    double T_dep = 6000.0;
    struct { double rho; const char *phase; } rho_tests[] = {
        {1e-30, "IGM (void)"},
        {1e-28, "CGM (halo)"},
        {1e-26, "Warm ISM"},
        {1e-24, "Diffuse ISM"},
        {1e-22, "Dense ISM"},
        {1e-20, "Molecular cloud"},
        {1e-18, "Star-forming core"}
    };
    int n_rho_tests = sizeof(rho_tests) / sizeof(rho_tests[0]);
    
    for (int i = 0; i < n_rho_tests; i++) {
        double mu = calculate_mu(rho_tests[i].rho, T_dep, Z_test);
        printf("    %.2e       %.4f   %s\n", rho_tests[i].rho, mu, rho_tests[i].phase);
    }
    
    /* Metallicity dependence */
    printf("\n###################################################################\n");
    printf("       Metallicity Dependence (rho = 1e-24 g/cm^3 T = 6000 K)\n");
    printf("###################################################################\n\n");
    printf("    Z [Z#]      mu      Comment\n");
    printf("    -------   ------   ----------------------\n");
    
    double Z_tests[] = {0.0, 0.01, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0};
    const char *Z_comments[] = {
        "Primordial (Pop III)",
        "Extremely metal-poor",
        "Metal-poor",
        "Sub-solar",
        "LMC-like",
        "Solar",
        "Super-solar",
        "High metallicity"
    };
    int n_Z_tests = sizeof(Z_tests) / sizeof(Z_tests[0]);
    
    for (int i = 0; i < n_Z_tests; i++) {
        double mu = calculate_mu(rho_test, T_dep, Z_tests[i]);
        printf("    %5.2f     %.4f   %s\n", Z_tests[i], mu, Z_comments[i]);
    }
    
    /* Save grid */
    printf("\n");
    save_grid_to_file("mu_grid.dat");
    
    /* Clean up */
    free_mu_grid();
    
    printf("\nTest complete.\n");
    return 0;
}

#endif /* TEST_GETMU */

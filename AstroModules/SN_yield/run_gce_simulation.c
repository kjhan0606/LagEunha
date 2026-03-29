/*
 * run_gce_simulation.c
 * * One-Zone Galaxy Chemical Evolution Simulation
 * utilizing the Cluster Supernova Library.
 * * Generates [Fe/H] vs [O/Fe] evolution track.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cluster_supernova.h"

/* Simulation Parameters */
#define MAX_STEPS 1000
#define TIME_STEP 2.0e7       // 20 Myr per step
#define MAX_TIME 10.0e9       // Run for 10 Gyr
#define SFR 1.0               // Constant Star Formation Rate (M_sun/yr)

/* Solar Abundance Reference (Asplund et al. 2009 / Grevesse & Sauval 1998 approximate)
 * Mass fractions (X, Y, Z_O, Z_Fe)
 */
#define SOLAR_Z 0.02
#define SOLAR_O 0.009618      // Mass fraction of Oxygen
#define SOLAR_FE 0.001296     // Mass fraction of Iron
#define SOLAR_H 0.70          // Mass fraction of Hydrogen

/* Structure to store history of star formation (Single Stellar Populations) */
typedef struct {
    double birth_time;
    double mass;
    double metallicity;
} SSP;

SSP history[MAX_STEPS];
int history_count = 0;

/* Global Gas Reservoir */
double gas_mass = 1.0e9;      // Initial Gas Reservoir (M_sun)
double mass_h = 0.75e9;       // Primordial H (75%)
double mass_he = 0.25e9;      // Primordial He (25%)
double mass_o = 0.0;          // Primordial O (0)
double mass_fe = 0.0;         // Primordial Fe (0)
double mass_z_total = 0.0;    // Total metals

/* Calculate [A/B] notation */
double calc_bracket(double mass_A, double mass_B, double solar_A, double solar_B) {
    if (mass_A <= 0 || mass_B <= 0) return -9.99; // Minimum floor
    
    // Convert mass ratio to number ratio? 
    // Standard notation [A/B] usually refers to number density log ratio relative to solar
    // [A/B] = log10( (N_A/N_B) / (N_A/N_B)_sun )
    //       = log10( (M_A*mu_B / M_B*mu_A) / (M_A_sun*mu_B / M_B_sun*mu_A) )
    //       = log10( (M_A/M_B) / (M_A_sun/M_B_sun) )
    
    double ratio = mass_A / mass_B;
    double solar_ratio = solar_A / solar_B;
    
    return log10(ratio / solar_ratio);
}

int main() {
    FILE *fp = fopen("gce_data.csv", "w");
    fprintf(fp, "Time_Myr,Z,Fe_H,O_Fe\n");
    
    printf("Running Galaxy Chemical Evolution Simulation...\n");
    printf("Total Time: %.1f Gyr, Step: %.1f Myr\n", MAX_TIME/1e9, TIME_STEP/1e6);
    
    double current_time = 0.0;
    
    for (int step = 0; step < MAX_STEPS; step++) {
        if (current_time >= MAX_TIME) break;
        
        /* 1. Calculate current Gas Metallicity (Z) */
        double total_gas = mass_h + mass_he + mass_z_total;
        double current_z = mass_z_total / total_gas;
        
        // Prevent Z from being exactly 0 for log calculations, but library handles 0 fine.
        if (current_z < 1.0e-8) current_z = 1.0e-8; 
        
        /* 2. Form Stars (New SSP) */
        double mass_formed = SFR * TIME_STEP;
        
        // Remove mass from gas (simplified: proportional removal)
        double h_frac = mass_h / total_gas;
        double he_frac = mass_he / total_gas;
        double z_frac = mass_z_total / total_gas;
        
        mass_h -= mass_formed * h_frac;
        mass_he -= mass_formed * he_frac;
        mass_z_total -= mass_formed * z_frac;
        // (Detailed O/Fe removal omitted for simplicity, assumes Z is well mixed)
        
        // Record SSP
        history[history_count].birth_time = current_time;
        history[history_count].mass = mass_formed;
        history[history_count].metallicity = current_z;
        history_count++;
        
        /* 3. Feedback Loop: Calculate yields from ALL past SSPs */
        // This simulates the delayed return from old stars
        double yield_h = 0, yield_he = 0, yield_o = 0, yield_fe = 0, yield_metals = 0;
        
        for (int i = 0; i < history_count; i++) {
            double age = current_time - history[i].birth_time;
            
            // Setup cluster properties for this specific past population
            ClusterProperties cluster = {
                .total_mass = history[i].mass,
                .age = age,
                .metallicity = history[i].metallicity,
                .imf_type = IMF_KROUPA,
                .sf_mode = SF_INSTANTANEOUS, // Treat as a burst
                .star_formation_rate = 0.0,
                // Enable Pop III for first generation
                .popiii_mode = (history[i].metallicity < 1e-4) ? POPIII_ENABLED : POPIII_DISABLED,
                .popiii_imf_max = 260.0,
                .popiii_mass_fraction = 0.1
            };
            
            ClusterOutput output;
            if (calculate_cluster_supernovae(&cluster, &output) == 0) {
                // The library returns rates/yields per year. Multiply by TIME_STEP.
                // NOTE: This assumes rates are constant over TIME_STEP.
                // For very young ages, TIME_STEP should be small.
                
                yield_h += output.yields.element_yields[H_ELEM] * TIME_STEP;
                yield_he += output.yields.element_yields[HE_ELEM] * TIME_STEP;
                yield_o += output.yields.element_yields[O_ELEM] * TIME_STEP;
                yield_fe += output.yields.element_yields[FE_ELEM] * TIME_STEP;
                yield_metals += output.yields.total_metals * TIME_STEP;
            }
        }
        
        /* 4. Update Gas Reservoir */
        mass_h += yield_h;
        mass_he += yield_he;
        mass_o += yield_o;
        mass_fe += yield_fe;
        mass_z_total += yield_metals;
        
        /* 5. Log Data */
        // [Fe/H] = log10( (M_Fe/M_H) / (M_Fe_sun/M_H_sun) )
        double fe_h = calc_bracket(mass_fe, mass_h, SOLAR_FE, SOLAR_H);
        
        // [O/Fe] = log10( (M_O/M_Fe) / (M_O_sun/M_Fe_sun) )
        double o_fe = calc_bracket(mass_o, mass_fe, SOLAR_O, SOLAR_FE);
        
        // Only log if we have enough metals to make sense
        if (mass_fe > 1.0e-10) {
            fprintf(fp, "%.2f,%.6e,%.3f,%.3f\n", 
                    current_time / 1.0e6, current_z, fe_h, o_fe);
        }
        
        current_time += TIME_STEP;
        
        if (step % 50 == 0) printf("Step %d: Age %.0f Myr, Z=%.4f, [Fe/H]=%.2f\n", 
                                   step, current_time/1e6, current_z, fe_h);
    }
    
    fclose(fp);
    printf("Simulation Complete. Data saved to 'gce_data.csv'\n");
    return 0;
}

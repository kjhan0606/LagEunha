/**
 * star_formation_advanced_example.c
 * 
 * Example program demonstrating advanced star formation features:
 * - Magnetic field effects
 * - Metallicity-dependent IMF
 * - Binary star populations
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "star_formation_advanced.h"

/**
 * Print detailed star particle information including binaries
 */
void print_advanced_star_info(const StarParticle *sp) {
    printf("\n=== Star Particle Information ===\n");
    printf("Total stellar mass: %.4e M_sun\n", sp->mass);
    printf("Number of stars: %d\n", sp->num_stars);
    printf("Number of binaries: %d (%.1f%%)\n", 
           sp->num_binaries, 100.0 * sp->num_binaries / sp->num_stars);
    printf("Metallicity: %.4f Z_sun\n", sp->metallicity / METALLICITY_SOLAR);
    
    if (sp->num_stars > 0) {
        /* Mass distribution statistics */
        double min_mass = sp->stellar_masses[0];
        double max_mass = sp->stellar_masses[0];
        double avg_mass = 0.0;
        
        for (int i = 0; i < sp->num_stars; i++) {
            if (sp->stellar_masses[i] < min_mass) min_mass = sp->stellar_masses[i];
            if (sp->stellar_masses[i] > max_mass) max_mass = sp->stellar_masses[i];
            avg_mass += sp->stellar_masses[i];
        }
        avg_mass /= sp->num_stars;
        
        printf("\nMass distribution:\n");
        printf("  Min mass: %.4f M_sun\n", min_mass);
        printf("  Max mass: %.4f M_sun\n", max_mass);
        printf("  Avg mass: %.4f M_sun\n", avg_mass);
        
        /* Binary statistics */
        if (sp->num_binaries > 0) {
            double avg_q = 0.0, avg_period = 0.0, avg_ecc = 0.0;
            int count = 0;
            
            for (int i = 0; i < sp->num_stars; i++) {
                if (sp->is_binary[i]) {
                    avg_q += sp->binary_params[i].mass_ratio;
                    avg_period += sp->binary_params[i].period_days;
                    avg_ecc += sp->binary_params[i].eccentricity;
                    count++;
                }
            }
            
            avg_q /= count;
            avg_period /= count;
            avg_ecc /= count;
            
            printf("\nBinary statistics:\n");
            printf("  Avg mass ratio (q): %.3f\n", avg_q);
            printf("  Avg period: %.2e days (%.2f yr)\n", avg_period, avg_period / 365.25);
            printf("  Avg eccentricity: %.3f\n", avg_ecc);
        }
    }
}

/**
 * Example 1: Magnetic field effects on star formation
 */
void example_magnetic_field_effects(void) {
    printf("\n========================================\n");
    printf("Example 1: Magnetic Field Effects\n");
    printf("========================================\n");
    
    GasElement gas;
    gas.density = 1.0e-20;
    gas.temperature = 30.0;
    gas.metallicity = METALLICITY_SOLAR;
    gas.velocity_dispersion = 1.0e5;
    gas.volume = pow(10.0 * PARSEC, 3);
    
    double timestep = 1.0e6;
    
    printf("\nComparing SF with different magnetic field strengths:\n\n");
    
    /* Case 1: No magnetic field */
    gas.magnetic_field_strength = 0.0;
    gas.magnetic_pressure = 0.0;
    gas.mass_to_flux_ratio = 10.0; /* Highly supercritical */
    
    printf("--- Case 1: No Magnetic Field ---\n");
    double sfr_nomag = calculate_sfr_magnetized(&gas, timestep);
    printf("SFR: %.4e M_sun/yr\n", sfr_nomag);
    
    /* Case 2: Weak magnetic field (β >> 1, thermal pressure dominated) */
    gas.magnetic_field_strength = estimate_magnetic_field(gas.density, gas.velocity_dispersion * 0.1);
    gas.magnetic_pressure = gas.magnetic_field_strength * gas.magnetic_field_strength / (8.0 * M_PI);
    gas.mass_to_flux_ratio = 5.0; /* Supercritical */
    
    printf("\n--- Case 2: Weak B Field (β ~ 50) ---\n");
    printf("B field: %.2e Gauss\n", gas.magnetic_field_strength);
    double sfr_weak = calculate_sfr_magnetized(&gas, timestep);
    printf("SFR: %.4e M_sun/yr (%.1f%% of no-field case)\n", 
           sfr_weak, 100.0 * sfr_weak / sfr_nomag);
    
    /* Case 3: Strong magnetic field (β ~ 1, equipartition) */
    gas.magnetic_field_strength = estimate_magnetic_field(gas.density, gas.velocity_dispersion);
    gas.magnetic_pressure = gas.magnetic_field_strength * gas.magnetic_field_strength / (8.0 * M_PI);
    gas.mass_to_flux_ratio = 2.0; /* Marginally supercritical */
    
    printf("\n--- Case 3: Strong B Field (β ~ 1, equipartition) ---\n");
    printf("B field: %.2e Gauss\n", gas.magnetic_field_strength);
    double sfr_strong = calculate_sfr_magnetized(&gas, timestep);
    printf("SFR: %.4e M_sun/yr (%.1f%% of no-field case)\n", 
           sfr_strong, 100.0 * sfr_strong / sfr_nomag);
    
    /* Case 4: Very strong magnetic field (β < 1, magnetically dominated) */
    gas.magnetic_field_strength = estimate_magnetic_field(gas.density, gas.velocity_dispersion * 2.0);
    gas.magnetic_pressure = gas.magnetic_field_strength * gas.magnetic_field_strength / (8.0 * M_PI);
    gas.mass_to_flux_ratio = 1.5; /* Marginally supercritical */
    
    printf("\n--- Case 4: Very Strong B Field (β ~ 0.25) ---\n");
    printf("B field: %.2e Gauss\n", gas.magnetic_field_strength);
    double sfr_vstrong = calculate_sfr_magnetized(&gas, timestep);
    printf("SFR: %.4e M_sun/yr (%.1f%% of no-field case)\n", 
           sfr_vstrong, 100.0 * sfr_vstrong / sfr_nomag);
    
    /* Case 5: Magnetically subcritical (cannot collapse) */
    gas.mass_to_flux_ratio = 0.5; /* Subcritical */
    
    printf("\n--- Case 5: Subcritical (λ < 1) ---\n");
    printf("Mass-to-flux ratio: %.2f (subcritical)\n", gas.mass_to_flux_ratio);
    double sfr_subcrit = calculate_sfr_magnetized(&gas, timestep);
    printf("SFR: %.4e M_sun/yr (%.1f%% of no-field case)\n", 
           sfr_subcrit, 100.0 * sfr_subcrit / sfr_nomag);
}

/**
 * Example 2: Metallicity-dependent IMF
 */
void example_metallicity_dependent_imf(void) {
    printf("\n\n========================================\n");
    printf("Example 2: Metallicity-Dependent IMF\n");
    printf("========================================\n");
    
    GasElement gas;
    gas.density = 1.0e-20;
    gas.temperature = 25.0;
    gas.velocity_dispersion = 1.2e5;
    gas.volume = pow(8.0 * PARSEC, 3);
    gas.magnetic_field_strength = estimate_magnetic_field(gas.density, gas.velocity_dispersion);
    gas.magnetic_pressure = gas.magnetic_field_strength * gas.magnetic_field_strength / (8.0 * M_PI);
    gas.mass_to_flux_ratio = 3.0;
    
    double timestep = 1.0e6;
    
    double metallicities[] = {0.001, 0.01, 0.1, 1.0, 3.0};
    const char* z_labels[] = {"0.001 Z_sun (primordial)", "0.01 Z_sun (very metal-poor)", 
                             "0.1 Z_sun (metal-poor)", "1.0 Z_sun (solar)", 
                             "3.0 Z_sun (metal-rich)"};
    
    printf("\nIMF variations with metallicity:\n");
    
    for (int i = 0; i < 5; i++) {
        gas.metallicity = metallicities[i] * METALLICITY_SOLAR;
        
        printf("\n--- %s ---\n", z_labels[i]);
        
        StarParticle sp;
        if (form_stars_advanced(&gas, timestep, IMF_METALLICITY_DEPENDENT, false, &sp)) {
            /* Calculate characteristic mass */
            double length_scale_pc = pow(gas.volume, 1.0/3.0) / PARSEC;
            double gas_mass_msun = (gas.density * gas.volume) / SOLAR_MASS;
            double surface_density = gas_mass_msun / (length_scale_pc * length_scale_pc);
            double m_char = calculate_characteristic_mass(metallicities[i], surface_density);
            
            printf("Characteristic mass: %.3f M_sun\n", m_char);
            
            /* Count massive vs low-mass stars */
            int low_mass = 0, intermediate = 0, high_mass = 0;
            for (int j = 0; j < sp.num_stars; j++) {
                if (sp.stellar_masses[j] < 0.5) low_mass++;
                else if (sp.stellar_masses[j] < 8.0) intermediate++;
                else high_mass++;
            }
            
            printf("Distribution: < 0.5 M_sun: %d (%.1f%%), 0.5-8 M_sun: %d (%.1f%%), > 8 M_sun: %d (%.1f%%)\n",
                   low_mass, 100.0 * low_mass / sp.num_stars,
                   intermediate, 100.0 * intermediate / sp.num_stars,
                   high_mass, 100.0 * high_mass / sp.num_stars);
            
            free_star_particle_advanced(&sp);
        }
    }
}

/**
 * Example 3: Binary star populations
 */
void example_binary_populations(void) {
    printf("\n\n========================================\n");
    printf("Example 3: Binary Star Populations\n");
    printf("========================================\n");
    
    GasElement gas;
    gas.density = 1.0e-20;
    gas.temperature = 30.0;
    gas.metallicity = METALLICITY_SOLAR;
    gas.velocity_dispersion = 1.0e5;
    gas.volume = pow(10.0 * PARSEC, 3);
    gas.magnetic_field_strength = estimate_magnetic_field(gas.density, gas.velocity_dispersion);
    gas.magnetic_pressure = gas.magnetic_field_strength * gas.magnetic_field_strength / (8.0 * M_PI);
    gas.mass_to_flux_ratio = 3.0;
    
    double timestep = 1.0e6;
    
    printf("\n--- Star Formation WITH Binaries ---\n");
    StarParticle sp_binary;
    if (form_stars_advanced(&gas, timestep, IMF_CHABRIER_2003, true, &sp_binary)) {
        print_advanced_star_info(&sp_binary);
        
        /* Analyze binary fraction by mass */
        printf("\nBinary fraction by mass range:\n");
        
        int count[4] = {0}, binary_count[4] = {0};
        for (int i = 0; i < sp_binary.num_stars; i++) {
            double m = sp_binary.stellar_masses[i];
            int bin = 0;
            if (m >= 0.5 && m < 1.5) bin = 1;
            else if (m >= 1.5 && m < 5.0) bin = 2;
            else if (m >= 5.0) bin = 3;
            
            count[bin]++;
            if (sp_binary.is_binary[i]) binary_count[bin]++;
        }
        
        const char* mass_ranges[] = {"< 0.5 M_sun", "0.5-1.5 M_sun", 
                                    "1.5-5 M_sun", "> 5 M_sun"};
        for (int i = 0; i < 4; i++) {
            if (count[i] > 0) {
                printf("  %s: %d/%d (%.1f%%)\n", mass_ranges[i], 
                       binary_count[i], count[i], 
                       100.0 * binary_count[i] / count[i]);
            }
        }
        
        /* Show some example binary systems */
        printf("\nExample binary systems:\n");
        int shown = 0;
        for (int i = 0; i < sp_binary.num_stars && shown < 5; i++) {
            if (sp_binary.is_binary[i]) {
                double m1 = sp_binary.stellar_masses[i];
                double m2 = m1 * sp_binary.binary_params[i].mass_ratio;
                printf("  System %d: M1=%.3f M_sun, M2=%.3f M_sun, P=%.2e d, e=%.3f, a=%.2f AU\n",
                       shown + 1, m1, m2, 
                       sp_binary.binary_params[i].period_days,
                       sp_binary.binary_params[i].eccentricity,
                       sp_binary.binary_params[i].separation_au);
                shown++;
            }
        }
        
        free_star_particle_advanced(&sp_binary);
    }
    
    printf("\n\n--- Star Formation WITHOUT Binaries (for comparison) ---\n");
    StarParticle sp_single;
    if (form_stars_advanced(&gas, timestep, IMF_CHABRIER_2003, false, &sp_single)) {
        printf("Total stellar mass: %.4e M_sun\n", sp_single.mass);
        printf("Number of stars: %d\n", sp_single.num_stars);
        printf("Number of binaries: %d\n", sp_single.num_binaries);
        free_star_particle_advanced(&sp_single);
    }
}

/**
 * Example 4: Combined effects - low metallicity with strong B field
 */
void example_combined_effects(void) {
    printf("\n\n========================================\n");
    printf("Example 4: Combined Effects\n");
    printf("========================================\n");
    printf("Low metallicity + Strong magnetic field + Binary formation\n\n");
    
    GasElement gas;
    gas.density = 5.0e-21;
    gas.temperature = 50.0;
    gas.metallicity = 0.01 * METALLICITY_SOLAR; /* Metal-poor */
    gas.velocity_dispersion = 1.5e5;
    gas.volume = pow(12.0 * PARSEC, 3);
    
    /* Strong magnetic field */
    gas.magnetic_field_strength = estimate_magnetic_field(gas.density, gas.velocity_dispersion * 1.5);
    gas.magnetic_pressure = gas.magnetic_field_strength * gas.magnetic_field_strength / (8.0 * M_PI);
    gas.mass_to_flux_ratio = 2.5;
    
    double timestep = 1.0e6;
    
    printf("Gas properties:\n");
    printf("  Metallicity: %.4f Z_sun\n", gas.metallicity / METALLICITY_SOLAR);
    printf("  B field: %.2e G\n", gas.magnetic_field_strength);
    printf("  Mass-to-flux ratio: %.2f\n", gas.mass_to_flux_ratio);
    
    double thermal_pressure = (gas.density / PROTON_MASS) * BOLTZMANN_CONSTANT * gas.temperature;
    double beta = thermal_pressure / gas.magnetic_pressure;
    printf("  Plasma beta: %.2f\n", beta);
    
    StarParticle sp;
    if (form_stars_advanced(&gas, timestep, IMF_METALLICITY_DEPENDENT, true, &sp)) {
        print_advanced_star_info(&sp);
        
        /* Expected: 
         * - Top-heavy IMF due to low metallicity
         * - Reduced SFR due to magnetic field
         * - Lower binary fraction for low-mass stars
         */
    }
}

int main(void) {
    /* Seed random number generator */
    srand(time(NULL));
    
    printf("===========================================\n");
    printf("Advanced Star Formation Module - Examples\n");
    printf("===========================================\n");
    printf("\nDemonstrating:\n");
    printf("  - Magnetic field effects (Pattle+ 2022, Kretschmer+ 2020)\n");
    printf("  - Metallicity-dependent IMF (Tanvir+ 2024, Yan+ 2018)\n");
    printf("  - Binary star populations (Moe+ 2017)\n");
    
    /* Run examples */
    example_magnetic_field_effects();
    example_metallicity_dependent_imf();
    example_binary_populations();
    example_combined_effects();
    
    printf("\n\n===========================================\n");
    printf("All examples completed successfully!\n");
    printf("===========================================\n");
    
    return 0;
}

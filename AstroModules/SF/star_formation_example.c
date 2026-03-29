/**
 * star_formation_example.c
 * 
 * Example program demonstrating the use of the star formation module
 * for galaxy formation simulations
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "star_formation.h"

/**
 * Print star particle information
 */
void print_star_particle_info(const StarParticle *sp) {
    printf("\n=== Star Particle Information ===\n");
    printf("Total stellar mass: %.4e M_sun\n", sp->mass);
    printf("Number of stars: %d\n", sp->num_stars);
    printf("Metallicity: %.4f Z_sun\n", sp->metallicity / METALLICITY_SOLAR);
    printf("Age: %.2f Myr\n", sp->age / 1.0e6);
    printf("Velocity: (%.2e, %.2e, %.2e) cm/s\n", 
           sp->velocity[0], sp->velocity[1], sp->velocity[2]);
    
    /* Print mass distribution statistics */
    if (sp->num_stars > 0) {
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
        
        /* Count stars in different mass bins */
        int count_low = 0, count_mid = 0, count_high = 0;
        for (int i = 0; i < sp->num_stars; i++) {
            if (sp->stellar_masses[i] < 0.5) count_low++;
            else if (sp->stellar_masses[i] < 8.0) count_mid++;
            else count_high++;
        }
        
        printf("\nMass bins:\n");
        printf("  < 0.5 M_sun (low-mass): %d (%.1f%%)\n", 
               count_low, 100.0 * count_low / sp->num_stars);
        printf("  0.5-8 M_sun (intermediate): %d (%.1f%%)\n", 
               count_mid, 100.0 * count_mid / sp->num_stars);
        printf("  > 8 M_sun (massive): %d (%.1f%%)\n", 
               count_high, 100.0 * count_high / sp->num_stars);
    }
}

/**
 * Example 1: Dense molecular cloud
 */
void example_dense_cloud(void) {
    printf("\n========================================\n");
    printf("Example 1: Dense Molecular Cloud\n");
    printf("========================================\n");
    
    /* Set up a dense, cold molecular cloud */
    GasElement gas;
    gas.density = 1.0e-20;              /* ~1000 H atoms/cm^3 */
    gas.temperature = 20.0;             /* 20 K - cold molecular gas */
    gas.metallicity = METALLICITY_SOLAR; /* Solar metallicity */
    gas.velocity_dispersion = 1.0e5;    /* ~1 km/s turbulent velocity */
    gas.volume = pow(10.0 * PARSEC, 3); /* 10 pc cube */
    
    printf("\nGas properties:\n");
    printf("  Density: %.2e g/cm^3 (%.1f H/cm^3)\n", 
           gas.density, gas.density / PROTON_MASS);
    printf("  Temperature: %.1f K\n", gas.temperature);
    printf("  Metallicity: %.4f Z_sun\n", gas.metallicity / METALLICITY_SOLAR);
    printf("  Velocity dispersion: %.2e cm/s (%.2f km/s)\n", 
           gas.velocity_dispersion, gas.velocity_dispersion / 1.0e5);
    printf("  Volume: %.2e cm^3 (%.1f pc^3)\n", 
           gas.volume, gas.volume / pow(PARSEC, 3));
    
    double gas_mass_msun = (gas.density * gas.volume) / SOLAR_MASS;
    printf("  Gas mass: %.2e M_sun\n", gas_mass_msun);
    
    /* Calculate derived quantities */
    double t_ff = calculate_freefall_time(gas.density);
    double m_jeans = calculate_jeans_mass(gas.density, gas.temperature, 1.3);
    
    printf("\nDerived quantities:\n");
    printf("  Free-fall time: %.2e yr (%.2f Myr)\n", t_ff, t_ff / 1.0e6);
    printf("  Jeans mass: %.2e M_sun\n", m_jeans);
    
    /* Form stars over 1 Myr timestep */
    double timestep = 1.0e6; /* 1 Myr */
    StarParticle sp;
    
    printf("\n--- Testing Chabrier IMF ---");
    if (form_stars(&gas, timestep, IMF_CHABRIER_2003, &sp)) {
        print_star_particle_info(&sp);
        free_star_particle(&sp);
    }
}

/**
 * Example 2: Comparison of different IMFs
 */
void example_imf_comparison(void) {
    printf("\n\n========================================\n");
    printf("Example 2: IMF Comparison\n");
    printf("========================================\n");
    
    /* Same gas properties for all tests */
    GasElement gas;
    gas.density = 5.0e-21;
    gas.temperature = 50.0;
    gas.metallicity = METALLICITY_SOLAR;
    gas.velocity_dispersion = 2.0e5;
    gas.volume = pow(5.0 * PARSEC, 3);
    
    double timestep = 1.0e6;
    
    /* Test each IMF type */
    IMF_Type imf_types[] = {IMF_CHABRIER_2003, IMF_KROUPA_2001, IMF_SALPETER_1955};
    const char* imf_names[] = {"Chabrier (2003)", "Kroupa (2001)", "Salpeter (1955)"};
    
    for (int i = 0; i < 3; i++) {
        printf("\n--- %s ---", imf_names[i]);
        
        StarParticle sp;
        if (form_stars(&gas, timestep, imf_types[i], &sp)) {
            print_star_particle_info(&sp);
            free_star_particle(&sp);
        }
    }
    
    /* Calculate mean masses */
    printf("\n\nMean stellar masses:\n");
    for (int i = 0; i < 3; i++) {
        double mean_mass = calculate_mean_stellar_mass(imf_types[i]);
        printf("  %s: %.4f M_sun\n", imf_names[i], mean_mass);
    }
}

/**
 * Example 3: Effect of metallicity and temperature
 */
void example_parameter_effects(void) {
    printf("\n\n========================================\n");
    printf("Example 3: Parameter Effects on SF\n");
    printf("========================================\n");
    
    /* Test different metallicities */
    double metallicities[] = {0.01 * METALLICITY_SOLAR, 0.1 * METALLICITY_SOLAR, 
                             METALLICITY_SOLAR, 3.0 * METALLICITY_SOLAR};
    const char* z_labels[] = {"0.01 Z_sun", "0.1 Z_sun", "1.0 Z_sun", "3.0 Z_sun"};
    
    GasElement gas;
    gas.density = 1.0e-20;
    gas.temperature = 30.0;
    gas.velocity_dispersion = 1.5e5;
    gas.volume = pow(8.0 * PARSEC, 3);
    
    double timestep = 1.0e6;
    
    printf("\nSFR vs Metallicity:\n");
    for (int i = 0; i < 4; i++) {
        gas.metallicity = metallicities[i];
        double sfr = calculate_sfr_volumetric(&gas, timestep);
        printf("  %s: SFR = %.4e M_sun/yr\n", z_labels[i], sfr);
    }
    
    /* Test different temperatures */
    double temperatures[] = {10.0, 30.0, 100.0, 300.0, 1000.0};
    gas.metallicity = METALLICITY_SOLAR;
    
    printf("\nSFR vs Temperature:\n");
    for (int i = 0; i < 5; i++) {
        gas.temperature = temperatures[i];
        double sfr = calculate_sfr_volumetric(&gas, timestep);
        printf("  %.1f K: SFR = %.4e M_sun/yr\n", temperatures[i], sfr);
    }
}

/**
 * Example 4: Time evolution of star formation
 */
void example_time_evolution(void) {
    printf("\n\n========================================\n");
    printf("Example 4: Time Evolution\n");
    printf("========================================\n");
    
    GasElement gas;
    gas.density = 1.0e-20;
    gas.temperature = 25.0;
    gas.metallicity = METALLICITY_SOLAR;
    gas.velocity_dispersion = 1.2e5;
    gas.volume = pow(12.0 * PARSEC, 3);
    
    double initial_gas_mass = (gas.density * gas.volume) / SOLAR_MASS;
    double total_stellar_mass = 0.0;
    
    printf("\nInitial gas mass: %.2e M_sun\n", initial_gas_mass);
    printf("\nTime evolution (timestep = 0.1 Myr):\n");
    printf("Time(Myr)  SFR(M_sun/yr)  Stellar Mass  Gas Mass  SF Efficiency\n");
    printf("------------------------------------------------------------------\n");
    
    double timestep = 1.0e5; /* 0.1 Myr */
    double total_time = 0.0;
    
    for (int step = 0; step < 50; step++) {
        double sfr = calculate_sfr_volumetric(&gas, timestep);
        double stellar_mass_formed = sfr * timestep;
        
        total_stellar_mass += stellar_mass_formed;
        total_time += timestep;
        
        /* Update gas mass (simple depletion, no feedback) */
        double current_gas_mass = initial_gas_mass - total_stellar_mass;
        gas.density = (current_gas_mass * SOLAR_MASS) / gas.volume;
        
        double sf_efficiency = total_stellar_mass / initial_gas_mass;
        
        if (step % 5 == 0) {
            printf("%.2f       %.4e       %.2e      %.2e      %.3f\n",
                   total_time / 1.0e6, sfr, total_stellar_mass, 
                   current_gas_mass, sf_efficiency);
        }
        
        if (current_gas_mass < initial_gas_mass * 0.1) {
            printf("\n(Gas mostly depleted, stopping evolution)\n");
            break;
        }
    }
}

int main(void) {
    /* Seed random number generator */
    srand(time(NULL));
    
    printf("===========================================\n");
    printf("Star Formation Module - Example Programs\n");
    printf("===========================================\n");
    printf("\nThis demonstrates a comprehensive star formation\n");
    printf("implementation for galaxy simulations.\n");
    printf("\nKey features:\n");
    printf("  - Kennicutt-Schmidt law\n");
    printf("  - Volumetric SF with turbulence\n");
    printf("  - Multiple IMF options (Chabrier, Kroupa, Salpeter)\n");
    printf("  - Individual stellar mass sampling\n");
    printf("  - Metallicity and temperature dependence\n");
    
    /* Run examples */
    example_dense_cloud();
    example_imf_comparison();
    example_parameter_effects();
    example_time_evolution();
    
    printf("\n\n===========================================\n");
    printf("Examples completed successfully!\n");
    printf("===========================================\n");
    
    return 0;
}

/*
 * Example Usage of Stellar Cluster Supernova Library
 * 
 * This program demonstrates typical usage patterns for
 * galaxy formation simulations
 */

#include "cluster_supernova.h"

/* ===== Example 1: Basic Usage ===== */
void example_basic_usage() {
    printf("=======================================================\n");
    printf("EXAMPLE 1: Basic Usage\n");
    printf("=======================================================\n\n");
    
    /* Define a young star cluster */
    ClusterProperties cluster = {
        .total_mass = 1.0e5,      // 100,000 M_sun
        .age = 50e6,              // 50 million years
        .metallicity = 0.02,      // Solar metallicity
        .imf_type = IMF_KROUPA,    // Kroupa IMF
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    /* Calculate supernova properties */
    ClusterOutput output;
    int result = calculate_cluster_supernovae(&cluster, &output);
    
    if (result == 0) {
        printf("SUCCESS: Calculations completed\n\n");
        print_cluster_output(&output);
        
        /* Access specific values */
        printf("Quick Access Examples:\n");
        printf("  SN II rate: %.3e per year\n", output.rates.sn_ii_rate);
        printf("  Total energy: %.3e erg/yr\n", output.energy.total_energy);
        printf("  Oxygen yield: %.3e M_sun/yr\n", 
               output.yields.element_yields[O_ELEM]);
        printf("  [O/Fe] ratio: %.3f dex\n\n",
               log10(output.yields.element_yields[O_ELEM] / 
                     output.yields.element_yields[FE_ELEM]));
    }
}

/* ===== Example 2: Time Evolution ===== */
void example_time_evolution() {
    printf("=======================================================\n");
    printf("EXAMPLE 2: Time Evolution of Cluster\n");
    printf("=======================================================\n\n");
    
    printf("--- Part A: Instantaneous Burst ---\n");
    printf("(Shows when stars die at specific ages)\n\n");
    
    ClusterProperties cluster_burst = {
        .total_mass = 2.0e5,
        .metallicity = 0.02,
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0,
        .age = 0  // Will be updated
    };
    
    /* Time steps: 1, 10, 50, 100, 500, 1000, 5000 Myr */
    double ages[] = {1e6, 10e6, 50e6, 100e6, 500e6, 1000e6, 5000e6};
    int n_steps = 7;
    
    printf("%-12s %12s %12s %15s %12s\n", 
           "Age (Myr)", "SN Ia (/yr)", "SN II (/yr)", 
           "Energy (erg/yr)", "Metals (M☉/yr)");
    printf("------------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_steps; i++) {
        cluster_burst.age = ages[i];
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster_burst, &output) == 0) {
            printf("%-12.1f %12.3e %12.3e %15.3e %12.3e\n",
                   ages[i] / 1e6,
                   output.rates.sn_ia_rate,
                   output.rates.sn_ii_rate,
                   output.energy.total_energy,
                   output.yields.total_metals);
        }
    }
    
    printf("\nNote: In instantaneous burst, SNe only occur at specific ages\n");
    printf("      when stars reach end of life. Zero = no stars dying.\n\n");
    
    printf("--- Part B: Continuous Star Formation ---\n");
    printf("(Steady SN rate from ongoing star formation)\n\n");
    
    ClusterProperties cluster_cont = {
        .total_mass = 2.0e5,
        .metallicity = 0.02,
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_CONTINUOUS,
        .star_formation_rate = 0.002,  // 0.002 M_sun/yr → total mass in 100 Myr
        .age = 0  // Will be updated
    };
    
    printf("%-12s %12s %12s %15s %12s\n", 
           "Age (Myr)", "SN Ia (/yr)", "SN II (/yr)", 
           "Energy (erg/yr)", "Metals (M☉/yr)");
    printf("------------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_steps; i++) {
        cluster_cont.age = ages[i];
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster_cont, &output) == 0) {
            printf("%-12.1f %12.3e %12.3e %15.3e %12.3e\n",
                   ages[i] / 1e6,
                   output.rates.sn_ia_rate,
                   output.rates.sn_ii_rate,
                   output.energy.total_energy,
                   output.yields.total_metals);
        }
    }
    
    printf("\nNote: Continuous SF shows steady enrichment once SN start\n");
    printf("      (after ~10 Myr for Type II, ~40 Myr for Type Ia)\n\n");
}

/* ===== Example 3: Metallicity Dependence ===== */
void example_metallicity_dependence() {
    printf("=======================================================\n");
    printf("EXAMPLE 3: Metallicity Dependence\n");
    printf("=======================================================\n\n");
    
    ClusterProperties cluster = {
        .total_mass = 1.0e5,
        .age = 20e6,  // Type II dominated age
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0,
        .metallicity = 0  // Will be updated
    };
    
    double metallicities[] = {0.0001, 0.001, 0.01, 0.02, 0.05};
    const char* labels[] = {"Very low", "Low", "0.5 solar", "Solar", "2.5 solar"};
    int n_Z = 5;
    
    printf("%-12s %10s %12s %12s %10s %10s\n",
           "Z", "Label", "O (M☉/yr)", "Fe (M☉/yr)", "[O/Fe]", "[N/O]");
    printf("-----------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_Z; i++) {
        cluster.metallicity = metallicities[i];
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster, &output) == 0) {
            double O_y = output.yields.element_yields[O_ELEM];
            double Fe_y = output.yields.element_yields[FE_ELEM];
            double N_y = output.yields.element_yields[N_ELEM];
            
            printf("%-12.4f %-10s %12.3e %12.3e %10.3f %10.3f\n",
                   metallicities[i], labels[i],
                   O_y, Fe_y,
                   log10(O_y / Fe_y),
                   log10(N_y / O_y));
        }
    }
    printf("\n");
}

/* ===== Example 4: IMF Effects ===== */
void example_imf_effects() {
    printf("=======================================================\n");
    printf("EXAMPLE 4: IMF Effects\n");
    printf("=======================================================\n\n");
    
    ClusterProperties cluster = {
        .total_mass = 1.0e5,
        .age = 30e6,
        .metallicity = 0.02,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0,
        .imf_type = IMF_SALPETER  // Will be changed
    };
    
    IMFType imfs[] = {IMF_SALPETER, IMF_KROUPA, IMF_CHABRIER};
    const char* imf_names[] = {"Salpeter", "Kroupa", "Chabrier"};
    
    printf("Different IMFs give different SN rates and yields:\n\n");
    printf("%-12s %12s %12s %12s\n", 
           "IMF", "SN II (/yr)", "O yield", "Fe yield");
    printf("-------------------------------------------------------\n");
    
    for (int i = 0; i < 3; i++) {
        cluster.imf_type = imfs[i];
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster, &output) == 0) {
            printf("%-12s %12.3e %12.3e %12.3e\n",
                   imf_names[i],
                   output.rates.sn_ii_rate,
                   output.yields.element_yields[O_ELEM],
                   output.yields.element_yields[FE_ELEM]);
        }
    }
    
    printf("\nNote: Chabrier produces more SNe (more low-mass stars, "
           "top-heavy at SN masses)\n\n");
}

/* ===== Example 5: Integration in Simulation ===== */
void example_simulation_integration() {
    printf("=======================================================\n");
    printf("EXAMPLE 5: Integration in Galaxy Simulation\n");
    printf("=======================================================\n\n");
    
    printf("Pseudo-code for typical simulation usage:\n\n");
    
    printf("// For each simulation timestep\n");
    printf("for (timestep = 0; timestep < n_timesteps; timestep++) {\n");
    printf("    current_time = timestep * dt;\n");
    printf("    \n");
    printf("    // Loop over all star clusters\n");
    printf("    for (i = 0; i < n_clusters; i++) {\n");
    printf("        cluster.total_mass = cluster_mass[i];\n");
    printf("        cluster.age = current_time - formation_time[i];\n");
    printf("        cluster.metallicity = get_local_Z(cluster_pos[i]);\n");
    printf("        cluster.imf_type = IMF_KROUPA;\n");
    printf("        \n");
    printf("        // Calculate supernova properties\n");
    printf("        calculate_cluster_supernovae(&cluster, &output);\n");
    printf("        \n");
    printf("        // Add thermal energy to gas\n");
    printf("        inject_energy(cluster_pos[i], output.energy.total_energy * dt);\n");
    printf("        \n");
    printf("        // Enrich gas with metals\n");
    printf("        for (elem = 0; elem < N_ELEMENTS; elem++) {\n");
    printf("            mass = output.yields.element_yields[elem] * dt;\n");
    printf("            enrich_gas(cluster_pos[i], elem, mass);\n");
    printf("        }\n");
    printf("        \n");
    printf("        // Update local metallicity\n");
    printf("        update_metallicity(cluster_pos[i]);\n");
    printf("    }\n");
    printf("}\n\n");
    
    /* Demonstrate with actual calculation */
    printf("Example calculation for one timestep:\n\n");
    
    ClusterProperties cluster = {
        .total_mass = 5.0e4,
        .age = 25e6,
        .metallicity = 0.015,
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    ClusterOutput output;
    calculate_cluster_supernovae(&cluster, &output);
    
    double dt = 1.0e6;  // 1 Myr timestep
    
    printf("Cluster: %.1e M_sun, age %.1f Myr, Z = %.4f\n",
           cluster.total_mass, cluster.age/1e6, cluster.metallicity);
    printf("Timestep: %.1f Myr\n\n", dt/1e6);
    
    printf("Energy injected: %.3e erg\n", output.energy.total_energy * dt);
    printf("Total mass ejected: %.3e M_sun\n", output.mass.total_mass * dt);
    printf("Oxygen added: %.3e M_sun\n", 
           output.yields.element_yields[O_ELEM] * dt);
    printf("Iron added: %.3e M_sun\n\n",
           output.yields.element_yields[FE_ELEM] * dt);
}

/* ===== Example 6: Chemical Evolution Tracking ===== */
void example_chemical_evolution() {
    printf("=======================================================\n");
    printf("EXAMPLE 6: Chemical Evolution Tracking\n");
    printf("=======================================================\n\n");
    
    printf("Tracking [O/Fe] and [N/O] evolution:\n\n");
    
    /* Simulate chemical evolution of a gas cloud */
    double gas_mass = 1.0e7;  // Initial gas mass
    double Z_gas = 0.001;     // Initial metallicity
    double O_mass = 0.0;      // Oxygen mass
    double Fe_mass = 0.0;     // Iron mass
    double N_mass = 0.0;      // Nitrogen mass
    
    /* Star formation events at different times */
    double sf_times[] = {10e6, 50e6, 100e6, 500e6};
    double sf_masses[] = {1e5, 5e4, 3e4, 2e4};
    int n_sf = 4;
    
    printf("%-12s %12s %12s %12s %12s\n",
           "Time (Myr)", "Z_gas", "O/gas (%)", "[O/Fe]", "[N/O]");
    printf("------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_sf; i++) {
        /* Form cluster */
        ClusterProperties cluster = {
            .total_mass = sf_masses[i],
            .age = sf_times[i],
            .metallicity = Z_gas,
            .imf_type = IMF_KROUPA,
            .sf_mode = SF_INSTANTANEOUS,
            .star_formation_rate = 0.0
        };
        
        /* Calculate yields */
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster, &output) == 0) {
            /* Add metals to gas (integrate over cluster lifetime) */
            double lifetime = 50e6;  // Assume 50 Myr enrichment
            O_mass += output.yields.element_yields[O_ELEM] * lifetime;
            Fe_mass += output.yields.element_yields[FE_ELEM] * lifetime;
            N_mass += output.yields.element_yields[N_ELEM] * lifetime;
            
            /* Update gas metallicity */
            double total_metals = O_mass + Fe_mass + N_mass;  // Simplified
            Z_gas = total_metals / gas_mass;
            
            /* Calculate abundance ratios */
            double O_to_Fe = (Fe_mass > 0) ? log10(O_mass / Fe_mass) : 0.0;
            double N_to_O = (O_mass > 0) ? log10(N_mass / O_mass) : 0.0;
            
            printf("%-12.1f %12.6f %12.3f %12.3f %12.3f\n",
                   sf_times[i] / 1e6,
                   Z_gas,
                   100.0 * O_mass / gas_mass,
                   O_to_Fe,
                   N_to_O);
        }
    }
    printf("\n");
}

/* ===== Main ===== */
int main(int argc, char *argv[]) {
    printf("===========================================================\n");
    printf(" Stellar Cluster Supernova Library - Usage Examples\n");
    printf("===========================================================\n\n");
    
    example_basic_usage();
    example_time_evolution();
    example_metallicity_dependence();
    example_imf_effects();
    example_simulation_integration();
    example_chemical_evolution();
    
    printf("===========================================================\n");
    printf(" All Examples Completed\n");
    printf("===========================================================\n");
    
    return 0;
}

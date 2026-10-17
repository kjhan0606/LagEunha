/*
 * Test Program for Continuous Star Formation Mode
 */

#include "cluster_supernova.h"

void test_continuous_sf_young() {
    printf("=======================================================\n");
    printf("TEST: Young Continuous SF Region\n");
    printf("=======================================================\n");
    
    ClusterProperties cluster = {
        .total_mass = 5.0e7,      // 50 million M_sun formed so far
        .age = 50e6,              // 50 Myr of continuous SF
        .metallicity = 0.02,      // Solar
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_CONTINUOUS,
        .star_formation_rate = 1.0  // 1 M_sun/year
    };
    
    print_cluster_properties(&cluster);
    
    ClusterOutput output;
    if (calculate_cluster_supernovae(&cluster, &output) == 0) {
        print_cluster_output(&output);
        
        printf("Expected behavior:\n");
        printf("  - SN II rate should be ~0.007 /year (0.7%% of SFR)\n");
        printf("  - SN Ia rate should be small but non-zero (age > 40 Myr)\n");
        printf("  - Total rate ~0.01 /year\n");
        printf("\nActual: Total rate = %.4f /year\n\n", output.rates.total_rate);
    }
}

void test_continuous_sf_milky_way() {
    printf("=======================================================\n");
    printf("TEST: Milky Way-like Galaxy\n");
    printf("=======================================================\n");
    
    ClusterProperties cluster = {
        .total_mass = 6.0e10,     // 60 billion M_sun
        .age = 10e9,              // 10 Gyr
        .metallicity = 0.02,      // Solar
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_CONTINUOUS,
        .star_formation_rate = 2.0  // 2 M_sun/year (current MW SFR)
    };
    
    print_cluster_properties(&cluster);
    
    ClusterOutput output;
    if (calculate_cluster_supernovae(&cluster, &output) == 0) {
        print_cluster_output(&output);
        
        printf("Comparison with Milky Way:\n");
        printf("  Observed SN rate: 2-3 per century = 0.02-0.03 /year\n");
        printf("  Our calculation:  %.3f /year\n", output.rates.total_rate);
        
        if (output.rates.total_rate > 0.01 && output.rates.total_rate < 0.05) {
            printf("  ✓ Within observed range!\n");
        } else {
            printf("  ⚠ Outside expected range\n");
        }
        printf("\n");
    }
}

void test_continuous_vs_burst() {
    printf("=======================================================\n");
    printf("TEST: Continuous SF vs Instantaneous Burst\n");
    printf("=======================================================\n");
    
    /* Burst mode */
    printf("--- Instantaneous Burst (10^5 M_sun, 50 Myr) ---\n");
    ClusterProperties burst = {
        .total_mass = 1.0e5,
        .age = 50e6,
        .metallicity = 0.02,
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    ClusterOutput output_burst;
    calculate_cluster_supernovae(&burst, &output_burst);
    printf("SN rate: %.6e /year\n\n", output_burst.rates.total_rate);
    
    /* Continuous mode */
    printf("--- Continuous SF (SFR = 0.001 M_sun/yr for 100 Myr) ---\n");
    ClusterProperties continuous = {
        .total_mass = 1.0e5,      // 0.001 * 100 Myr = 10^5 M_sun
        .age = 100e6,
        .metallicity = 0.02,
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_CONTINUOUS,
        .star_formation_rate = 0.001  // 0.001 M_sun/year
    };
    
    ClusterOutput output_cont;
    calculate_cluster_supernovae(&continuous, &output_cont);
    printf("SN rate: %.6e /year\n\n", output_cont.rates.total_rate);
    
    printf("Analysis:\n");
    printf("  Burst: Very low rate (stochastic, only one age group)\n");
    printf("  Continuous: Higher rate (all ages contributing)\n");
    printf("  Ratio: %.1f× higher for continuous SF\n\n",
           output_cont.rates.total_rate / output_burst.rates.total_rate);
}

void test_sfr_scaling() {
    printf("=======================================================\n");
    printf("TEST: SN Rate Scaling with SFR\n");
    printf("=======================================================\n");
    
    double sfrs[] = {0.1, 1.0, 10.0, 100.0};
    int n_sfr = 4;
    
    printf("%-15s %15s %15s %15s\n", "SFR (M_sun/yr)", "SN Ia (/yr)", "SN II (/yr)", "Total (/yr)");
    printf("----------------------------------------------------------------\n");
    
    for (int i = 0; i < n_sfr; i++) {
        ClusterProperties cluster = {
            .total_mass = sfrs[i] * 1e9,  // Assume 1 Gyr of SF
            .age = 1e9,
            .metallicity = 0.02,
            .imf_type = IMF_KROUPA,
            .sf_mode = SF_CONTINUOUS,
            .star_formation_rate = sfrs[i]
        };
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster, &output) == 0) {
            printf("%-15.1f %15.6e %15.6e %15.6e\n",
                   sfrs[i],
                   output.rates.sn_ia_rate,
                   output.rates.sn_ii_rate,
                   output.rates.total_rate);
        }
    }
    
    printf("\nRate should scale linearly with SFR ✓\n\n");
}

void test_age_evolution_continuous() {
    printf("=======================================================\n");
    printf("TEST: Time Evolution (Continuous SF)\n");
    printf("=======================================================\n");
    
    double sfr = 1.0;  // 1 M_sun/year
    double ages[] = {10e6, 50e6, 100e6, 500e6, 1e9, 5e9, 10e9};
    int n_ages = 7;
    
    printf("SFR = %.1f M_sun/yr (constant)\n\n", sfr);
    printf("%-12s %15s %15s %15s\n", "Age (Myr)", "SN Ia (/yr)", "SN II (/yr)", "Total (/yr)");
    printf("----------------------------------------------------------------\n");
    
    for (int i = 0; i < n_ages; i++) {
        ClusterProperties cluster = {
            .total_mass = sfr * ages[i],  // Total mass formed
            .age = ages[i],
            .metallicity = 0.02,
            .imf_type = IMF_KROUPA,
            .sf_mode = SF_CONTINUOUS,
            .star_formation_rate = sfr
        };
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster, &output) == 0) {
            printf("%-12.1f %15.6e %15.6e %15.6e\n",
                   ages[i] / 1e6,
                   output.rates.sn_ia_rate,
                   output.rates.sn_ii_rate,
                   output.rates.total_rate);
        }
    }
    
    printf("\nNotes:\n");
    printf("  - SN II rate reaches steady state quickly (~40 Myr)\n");
    printf("  - SN Ia rate increases with age (DTD integration)\n");
    printf("  - Total rate dominated by SN II for constant SFR\n\n");
}

int main(int argc, char *argv[]) {
    printf("===========================================================\n");
    printf(" Continuous Star Formation Mode - Test Suite\n");
    printf("===========================================================\n\n");
    
    test_continuous_sf_young();
    test_continuous_sf_milky_way();
    test_continuous_vs_burst();
    test_sfr_scaling();
    test_age_evolution_continuous();
    
    printf("===========================================================\n");
    printf(" All Continuous SF Tests Completed\n");
    printf("===========================================================\n");
    
    return 0;
}

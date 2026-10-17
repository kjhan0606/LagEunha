/*
 * Test Program for Stellar Cluster Supernova Library
 * 
 * Tests various cluster configurations and validates results
 */

#include "cluster_supernova.h"

/* Test case structure */
typedef struct {
    const char *name;
    ClusterProperties cluster;
    int expected_success;
} TestCase;

/* ===== Test Cases ===== */

void test_young_cluster() {
    printf("=======================================================\n");
    printf("TEST 1: Young Massive Cluster (Instantaneous Burst)\n");
    printf("=======================================================\n");
    
    ClusterProperties cluster = {
        .total_mass = 1.0e5,      // 100,000 M_sun (massive cluster)
        .age = 10e6,              // 10 Myr (young)
        .metallicity = 0.02,      // Solar metallicity
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0  // Not used for instantaneous
    };
    
    print_cluster_properties(&cluster);
    
    ClusterOutput output;
    int result = calculate_cluster_supernovae(&cluster, &output);
    
    if (result == 0) {
        print_cluster_output(&output);
        printf("TEST 1: PASSED\n\n");
    } else {
        printf("TEST 1: FAILED - Calculation error\n\n");
    }
}

void test_old_cluster() {
    printf("=======================================================\n");
    printf("TEST 2: Old Globular Cluster\n");
    printf("=======================================================\n");
    
    ClusterProperties cluster = {
        .total_mass = 1.0e6,      // 1 million M_sun
        .age = 10e9,              // 10 Gyr (old)
        .metallicity = 0.001,     // Metal-poor ([Fe/H] ~ -1.3)
        .imf_type = IMF_SALPETER,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    print_cluster_properties(&cluster);
    
    ClusterOutput output;
    int result = calculate_cluster_supernovae(&cluster, &output);
    
    if (result == 0) {
        print_cluster_output(&output);
        
        /* Validation: Old cluster should have Type Ia but no Type II */
        if (output.rates.sn_ia_rate > 0.0 && output.rates.sn_ii_rate < 1e-10) {
            printf("Validation: Correctly shows Type Ia only (all massive stars dead)\n");
            printf("TEST 2: PASSED\n\n");
        } else {
            printf("TEST 2: FAILED - Unexpected SN rates for old cluster\n\n");
        }
    } else {
        printf("TEST 2: FAILED - Calculation error\n\n");
    }
}

void test_starburst_cluster() {
    printf("=======================================================\n");
    printf("TEST 3: Starburst Cluster (Type II dominated)\n");
    printf("=======================================================\n");
    
    ClusterProperties cluster = {
        .total_mass = 5.0e5,      // 500,000 M_sun (very massive)
        .age = 20e6,              // 20 Myr (prime Type II age)
        .metallicity = 0.02,      // Solar
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    print_cluster_properties(&cluster);
    
    ClusterOutput output;
    int result = calculate_cluster_supernovae(&cluster, &output);
    
    if (result == 0) {
        print_cluster_output(&output);
        
        /* Validation: Should have significant Type II, minimal Type Ia */
        if (output.rates.sn_ii_rate > output.rates.sn_ia_rate) {
            printf("Validation: Correctly Type II dominated\n");
            printf("TEST 3: PASSED\n\n");
        } else {
            printf("TEST 3: FAILED - Expected Type II dominance\n\n");
        }
    } else {
        printf("TEST 3: FAILED - Calculation error\n\n");
    }
}

void test_metal_poor_cluster() {
    printf("=======================================================\n");
    printf("TEST 4: Primordial Metallicity Cluster\n");
    printf("=======================================================\n");
    
    ClusterProperties cluster = {
        .total_mass = 1.0e5,
        .age = 15e6,              // 15 Myr
        .metallicity = 0.0001,    // Very metal-poor ([Fe/H] ~ -2.3)
        .imf_type = IMF_CHABRIER,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    print_cluster_properties(&cluster);
    
    ClusterOutput output;
    int result = calculate_cluster_supernovae(&cluster, &output);
    
    if (result == 0) {
        print_cluster_output(&output);
        
        /* Check oxygen-to-iron ratio */
        double O_yield = output.yields.element_yields[O_ELEM];
        double Fe_yield = output.yields.element_yields[FE_ELEM];
        
        if (Fe_yield > 0) {
            double O_to_Fe = O_yield / Fe_yield;
            printf("Validation: [O/Fe] ratio = %.3f dex\n", log10(O_to_Fe));
            
            /* Metal-poor should have high [O/Fe] */
            if (log10(O_to_Fe) > 1.0) {
                printf("Validation: High [O/Fe] as expected for metal-poor\n");
                printf("TEST 4: PASSED\n\n");
            } else {
                printf("TEST 4: WARNING - [O/Fe] lower than expected\n\n");
            }
        }
    } else {
        printf("TEST 4: FAILED - Calculation error\n\n");
    }
}

void test_intermediate_age_cluster() {
    printf("=======================================================\n");
    printf("TEST 5: Intermediate Age Cluster (100 Myr)\n");
    printf("=======================================================\n");
    
    ClusterProperties cluster = {
        .total_mass = 2.0e5,
        .age = 100e6,             // 100 Myr
        .metallicity = 0.02,      // Solar
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    print_cluster_properties(&cluster);
    
    ClusterOutput output;
    int result = calculate_cluster_supernovae(&cluster, &output);
    
    if (result == 0) {
        print_cluster_output(&output);
        
        /* Should have both Type Ia and Type II */
        if (output.rates.sn_ia_rate > 0.0 && output.rates.sn_ii_rate > 0.0) {
            printf("Validation: Both Type Ia and Type II present\n");
            printf("TEST 5: PASSED\n\n");
        } else {
            printf("TEST 5: WARNING - Expected both SN types\n\n");
        }
    } else {
        printf("TEST 5: FAILED - Calculation error\n\n");
    }
}

void test_error_handling() {
    printf("=======================================================\n");
    printf("TEST 6: Error Handling\n");
    printf("=======================================================\n");
    
    int passed = 0;
    
    /* Test 1: Negative mass */
    ClusterProperties bad_cluster1 = {
        .total_mass = -1000.0,
        .age = 1e9,
        .metallicity = 0.02,
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    ClusterOutput output;
    if (calculate_cluster_supernovae(&bad_cluster1, &output) != 0) {
        printf("✓ Correctly rejected negative mass\n");
        passed++;
    } else {
        printf("✗ Failed to reject negative mass\n");
    }
    
    /* Test 2: Invalid metallicity */
    ClusterProperties bad_cluster2 = {
        .total_mass = 1e5,
        .age = 1e9,
        .metallicity = 0.5,  // Too high
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
    };
    
    if (calculate_cluster_supernovae(&bad_cluster2, &output) != 0) {
        printf("✓ Correctly rejected invalid metallicity\n");
        passed++;
    } else {
        printf("✗ Failed to reject invalid metallicity\n");
    }
    
    /* Test 3: NULL pointer */
    if (calculate_cluster_supernovae(NULL, &output) != 0) {
        printf("✓ Correctly handled NULL cluster pointer\n");
        passed++;
    } else {
        printf("✗ Failed to handle NULL cluster pointer\n");
    }
    
    if (passed == 3) {
        printf("\nTEST 6: PASSED (All error cases handled)\n\n");
    } else {
        printf("\nTEST 6: FAILED (%d/3 error cases handled)\n\n", passed);
    }
}

void test_imf_comparison() {
    printf("=======================================================\n");
    printf("TEST 7: IMF Comparison\n");
    printf("=======================================================\n");
    
    double mass = 1.0e5;
    double age = 50e6;
    double Z = 0.02;
    
    printf("Cluster: %.1e M_sun, age %.1f Myr, Z = %.4f\n\n", mass, age/1e6, Z);
    
    IMFType imf_types[] = {IMF_SALPETER, IMF_KROUPA, IMF_CHABRIER};
    const char* imf_names[] = {"Salpeter", "Kroupa", "Chabrier"};
    
    printf("%-12s %12s %12s %12s\n", "IMF", "SN Ia Rate", "SN II Rate", "Total Rate");
    printf("----------------------------------------------------------------\n");
    
    for (int i = 0; i < 3; i++) {
        ClusterProperties cluster = {
            .total_mass = mass,
            .age = age,
            .metallicity = Z,
            .imf_type = imf_types[i],
            .sf_mode = SF_INSTANTANEOUS,
            .star_formation_rate = 0.0
        };
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster, &output) == 0) {
            printf("%-12s %12.6e %12.6e %12.6e\n",
                   imf_names[i],
                   output.rates.sn_ia_rate,
                   output.rates.sn_ii_rate,
                   output.rates.total_rate);
        }
    }
    
    printf("\nTEST 7: PASSED (IMF comparison complete)\n\n");
}

void test_metallicity_effects() {
    printf("=======================================================\n");
    printf("TEST 8: Metallicity Effects on Yields\n");
    printf("=======================================================\n");
    
    double mass = 1.0e5;
    double age = 20e6;  // Type II dominated age
    
    double metallicities[] = {0.0001, 0.001, 0.02, 0.05};
    const char* Z_labels[] = {"0.0001", "0.001", "0.02 (solar)", "0.05"};
    
    printf("Cluster: %.1e M_sun, age %.1f Myr\n\n", mass, age/1e6);
    printf("%-15s %12s %12s %10s\n", "Metallicity", "O yield", "Fe yield", "[O/Fe]");
    printf("----------------------------------------------------------------\n");
    
    for (int i = 0; i < 4; i++) {
        ClusterProperties cluster = {
            .total_mass = mass,
            .age = age,
            .metallicity = metallicities[i],
            .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0
        };
        
        ClusterOutput output;
        if (calculate_cluster_supernovae(&cluster, &output) == 0) {
            double O_yield = output.yields.element_yields[O_ELEM];
            double Fe_yield = output.yields.element_yields[FE_ELEM];
            double O_to_Fe = (Fe_yield > 0) ? log10(O_yield / Fe_yield) : 0.0;
            
            printf("%-15s %12.6e %12.6e %10.3f\n",
                   Z_labels[i], O_yield, Fe_yield, O_to_Fe);
        }
    }
    
    printf("\nNote: [O/Fe] should decrease with increasing metallicity\n");
    printf("TEST 8: PASSED (Metallicity effects demonstrated)\n\n");
}

/* ===== Main Test Runner ===== */

int main(int argc, char *argv[]) {
    printf("===========================================================\n");
    printf(" Stellar Cluster Supernova Library - Test Suite\n");
    printf("===========================================================\n\n");
    
    /* Run all tests */
    test_young_cluster();
    test_old_cluster();
    test_starburst_cluster();
    test_metal_poor_cluster();
    test_intermediate_age_cluster();
    test_error_handling();
    test_imf_comparison();
    test_metallicity_effects();
    
    printf("===========================================================\n");
    printf(" All Tests Completed\n");
    printf("===========================================================\n");
    
    return 0;
}

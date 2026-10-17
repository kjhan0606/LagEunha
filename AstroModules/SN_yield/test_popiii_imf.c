#include "cluster_supernova.h"
#include <stdio.h>

int main() {
    printf("=======================================================\n");
    printf(" Pop III IMF Comparison Test\n");
    printf("=======================================================\n\n");
    
    /* Test 1: Normal Kroupa IMF vs Pop III top-heavy IMF */
    printf("=== TEST 1: IMF Comparison (Z=0, Age=3 Myr) ===\n\n");
    
    /* Same cluster properties, different IMFs */
    ClusterProperties kroupa_cluster = {
        .total_mass = 1e6,
        .age = 3e6,
        .metallicity = 0.0,
        .imf_type = IMF_KROUPA,
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0,
        .popiii_mode = POPIII_ENABLED,
        .popiii_imf_max = 260.0,
        .popiii_mass_fraction = 0.10  // 10%
    };
    
    ClusterProperties popiii_cluster = kroupa_cluster;
    popiii_cluster.imf_type = IMF_POPIII_TOPHEAVY;
    
    ClusterOutput kroupa_output, popiii_output;
    
    printf("--- Kroupa IMF (standard) ---\n");
    print_cluster_properties(&kroupa_cluster);
    calculate_cluster_supernovae(&kroupa_cluster, &kroupa_output);
    print_cluster_output(&kroupa_output);
    
    printf("\n--- Pop III Top-Heavy IMF ---\n");
    print_cluster_properties(&popiii_cluster);
    calculate_cluster_supernovae(&popiii_cluster, &popiii_output);
    print_cluster_output(&popiii_output);
    
    /* Comparison table */
    printf("=== IMF Comparison Summary ===\n\n");
    printf("%-25s %15s %15s %12s\n", 
           "Property", "Kroupa IMF", "Pop III IMF", "Ratio");
    printf("------------------------------------------------------------------------\n");
    printf("%-25s %15.3e %15.3e %12.2f\n",
           "Type II rate (/yr)",
           kroupa_output.rates.sn_ii_rate,
           popiii_output.rates.sn_ii_rate,
           popiii_output.rates.sn_ii_rate / (kroupa_output.rates.sn_ii_rate + 1e-30));
    printf("%-25s %15.3e %15.3e %12.2f\n",
           "PISN rate (/yr)",
           kroupa_output.rates.pisn_rate,
           popiii_output.rates.pisn_rate,
           popiii_output.rates.pisn_rate / (kroupa_output.rates.pisn_rate + 1e-30));
    printf("%-25s %15.3e %15.3e %12.2f\n",
           "Total rate (/yr)",
           kroupa_output.rates.total_rate,
           popiii_output.rates.total_rate,
           popiii_output.rates.total_rate / (kroupa_output.rates.total_rate + 1e-30));
    printf("%-25s %15.3e %15.3e %12.2f\n",
           "Total energy (erg/yr)",
           kroupa_output.energy.total_energy,
           popiii_output.energy.total_energy,
           popiii_output.energy.total_energy / (kroupa_output.energy.total_energy + 1e-30));
    printf("%-25s %15.3e %15.3e %12.2f\n",
           "Fe yield (M#/yr)",
           kroupa_output.yields.element_yields[FE_ELEM],
           popiii_output.yields.element_yields[FE_ELEM],
           popiii_output.yields.element_yields[FE_ELEM] / 
           (kroupa_output.yields.element_yields[FE_ELEM] + 1e-30));
    printf("\n");
    
    /* Test 2: Different mass fractions with Pop III IMF */
    printf("=== TEST 2: Pop III Mass Fraction Sensitivity ===\n\n");
    
    double fractions[] = {0.01, 0.05, 0.10, 0.20, 0.30, 0.50};
    int n_fractions = 6;
    
    printf("%-15s %15s %15s %15s\n",
           "Mass Fraction", "PISN Rate", "Total Energy", "Fe Yield");
    printf("---------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_fractions; i++) {
        ClusterProperties test_cluster = popiii_cluster;
        test_cluster.popiii_mass_fraction = fractions[i];
        
        ClusterOutput test_output;
        calculate_cluster_supernovae(&test_cluster, &test_output);
        
        printf("%-15.0f%% %15.3e %15.3e %15.3e\n",
               fractions[i] * 100.0,
               test_output.rates.pisn_rate,
               test_output.energy.total_energy,
               test_output.yields.element_yields[FE_ELEM]);
    }
    printf("\n");
    
    /* Test 3: Different metallicities */
    printf("=== TEST 3: Metallicity Dependence ===\n\n");
    
    double metallicities[] = {0.0, 0.0001, 0.001, 0.02};
    const char* z_labels[] = {"Z=0 (primordial)", "Z=10^-4", "Z=10^-3", "Z=0.02 (solar)"};
    int n_metals = 4;
    
    printf("%-20s %15s %15s %12s\n",
           "Metallicity", "PISN Rate", "Type II Rate", "Ratio");
    printf("---------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_metals; i++) {
        ClusterProperties test_cluster = popiii_cluster;
        test_cluster.metallicity = metallicities[i];
        test_cluster.popiii_mass_fraction = 0.10;
        
        ClusterOutput test_output;
        calculate_cluster_supernovae(&test_cluster, &test_output);
        
        double ratio = 0.0;
        if (test_output.rates.pisn_rate > 0) {
            ratio = test_output.rates.pisn_rate / 
                   (test_output.rates.sn_ii_rate + 1e-30);
        }
        
        printf("%-20s %15.3e %15.3e %12.2f\n",
               z_labels[i],
               test_output.rates.pisn_rate,
               test_output.rates.sn_ii_rate,
               ratio);
    }
    printf("\n");
    printf("Note: PISN only occur at Z < 0.001\n\n");
    
    /* Test 4: Time evolution with Pop III IMF */
    printf("=== TEST 4: Time Evolution (Pop III IMF) ===\n\n");
    
    double ages[] = {1e6, 2e6, 3e6, 4e6, 5e6, 10e6};
    int n_ages = 6;
    
    printf("%-10s %15s %15s %15s\n",
           "Age (Myr)", "PISN Rate", "Total Rate", "PISN Fraction");
    printf("-------------------------------------------------------------\n");
    
    for (int i = 0; i < n_ages; i++) {
        ClusterProperties test_cluster = popiii_cluster;
        test_cluster.age = ages[i];
        test_cluster.popiii_mass_fraction = 0.20;  // 20% for more visible effect
        
        ClusterOutput test_output;
        calculate_cluster_supernovae(&test_cluster, &test_output);
        
        double pisn_frac = test_output.rates.pisn_rate / 
                          (test_output.rates.total_rate + 1e-30);
        
        printf("%-10.1f %15.3e %15.3e %15.3f\n",
               ages[i] / 1e6,
               test_output.rates.pisn_rate,
               test_output.rates.total_rate,
               pisn_frac);
    }
    printf("\n");
    
    /* Summary */
    printf("=======================================================\n");
    printf(" Summary:\n");
    printf("=======================================================\n\n");
    printf("Pop III Top-Heavy IMF produces:\n");
    printf("  - Higher PISN rates (more massive stars)\n");
    printf("  - Peak around 100 M# instead of ~0.3 M#\n");
    printf("  - Slope M^-1.0 (flatter than M^-2.3)\n");
    printf("  - Only for Z < 0.001 (primordial/metal-poor)\n\n");
    
    printf("Mass fraction recommendations:\n");
    printf("  - Conservative: 10%% (similar to normal SF)\n");
    printf("  - Moderate: 20%%-30%% (Pop III dominated)\n");
    printf("  - Extreme: 50%% (pure Pop III cluster)\n\n");
    
    printf("=======================================================\n");
    printf(" All Tests Completed!\n");
    printf("=======================================================\n");
    
    return 0;
}

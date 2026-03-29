/*
 * ============================================================================
 * STELLAR WIND YIELD LIBRARY - Test and Example Program
 * ============================================================================
 * 
 * This program demonstrates the usage of the stellar wind yield library
 * for calculating metal yields from stellar winds in galaxy simulations.
 * 
 * Compilation:
 *   make test_stellar_wind
 * 
 * Usage:
 *   ./test_stellar_wind [imf_type] [metallicity] [cluster_mass] [cluster_age]
 * 
 * Example:
 *   ./test_stellar_wind 2 0.02 1e6 1e7
 *   (Chabrier IMF, solar Z, 10^6 M_sun cluster, 10 Myr age)
 * 
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stellar_wind_lib.h"

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================ */

void print_separator(const char* title) {
    printf("\n");
    printf("================================================================================\n");
    printf("  %s\n", title);
    printf("================================================================================\n");
}

void print_subseparator(const char* title) {
    printf("\n--- %s ---\n", title);
}

void print_usage(const char* progname) {
    printf("Usage: %s [imf_type] [metallicity] [cluster_mass_Msun] [cluster_age_yr]\n\n", progname);
    printf("Arguments:\n");
    printf("  imf_type      : IMF model (0=Salpeter, 1=Kroupa, 2=Chabrier, 3=Baldry-Glazebrook)\n");
    printf("  metallicity   : Initial metallicity Z (e.g., 0.02 for solar, 0.001 for low-Z)\n");
    printf("  cluster_mass  : Initial cluster mass in M_sun (e.g., 1e6)\n");
    printf("  cluster_age   : Cluster age in years (e.g., 1e7 for 10 Myr)\n");
    printf("\nExample:\n");
    printf("  %s 2 0.02 1e6 1e7\n", progname);
    printf("  (Chabrier IMF, solar metallicity, 10^6 M_sun cluster at 10 Myr)\n");
}

/* ============================================================================
 * TEST FUNCTIONS
 * ============================================================================ */

/*
 * Test 1: Basic library initialization and utility functions
 */
void test_basic_functions(void) {
    print_separator("TEST 1: Basic Library Functions");
    
    /* Initialize library */
    SWL_ErrorCode err = swl_init();
    printf("Library initialization: %s\n", swl_get_error_message(err));
    
    /* Print element names */
    printf("\nTracked elements:\n");
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        printf("  [%2d] %-8s  Solar abundance: %.4e\n", 
               i, swl_get_element_name((SWL_Element)i),
               swl_get_solar_abundance((SWL_Element)i));
    }
    
    /* Print IMF names */
    printf("\nAvailable IMF models:\n");
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        printf("  [%d] %s\n", i, swl_get_imf_name((SWL_IMFType)i));
    }
}

/*
 * Test 2: Stellar lifetime and turnoff mass
 */
void test_stellar_lifetimes(double Z) {
    print_separator("TEST 2: Stellar Lifetimes and Turnoff Masses");
    
    printf("Metallicity: Z = %.4f (Z/Z_sun = %.2f)\n\n", Z, Z / SWL_Z_SUN);
    
    /* Stellar lifetimes */
    printf("Stellar Lifetimes:\n");
    printf("%-12s %-15s %-12s\n", "Mass (M_sun)", "Lifetime (yr)", "Lifetime");
    printf("%-12s %-15s %-12s\n", "-----------", "-------------", "--------");
    
    double masses[] = {0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 15.0, 25.0, 40.0, 60.0, 100.0, 150.0};
    int n_masses = sizeof(masses) / sizeof(masses[0]);
    
    for (int i = 0; i < n_masses; i++) {
        double tau = swl_stellar_lifetime(masses[i], Z);
        char time_str[32];
        if (tau >= 1e9) {
            sprintf(time_str, "%.2f Gyr", tau / 1e9);
        } else if (tau >= 1e6) {
            sprintf(time_str, "%.2f Myr", tau / 1e6);
        } else {
            sprintf(time_str, "%.2e yr", tau);
        }
        printf("%-12.1f %-15.4e %-12s\n", masses[i], tau, time_str);
    }
    
    /* Turnoff masses */
    printf("\nMain Sequence Turnoff Masses:\n");
    printf("%-15s %-15s\n", "Age", "M_TO (M_sun)");
    printf("%-15s %-15s\n", "---", "-----------");
    
    double ages[] = {1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 1e10};
    int n_ages = sizeof(ages) / sizeof(ages[0]);
    
    for (int i = 0; i < n_ages; i++) {
        double m_to = swl_turnoff_mass(ages[i], Z);
        char age_str[32];
        if (ages[i] >= 1e9) {
            sprintf(age_str, "%.1f Gyr", ages[i] / 1e9);
        } else if (ages[i] >= 1e6) {
            sprintf(age_str, "%.1f Myr", ages[i] / 1e6);
        } else {
            sprintf(age_str, "%.0e yr", ages[i]);
        }
        printf("%-15s %-15.3f\n", age_str, m_to);
    }
}

/*
 * Test 3: Single star wind yields
 */
void test_single_star_yields(double Z) {
    print_separator("TEST 3: Single Star Wind Yields");
    
    printf("Metallicity: Z = %.4f\n\n", Z);
    
    double masses[] = {1.5, 3.0, 5.0, 7.0, 15.0, 25.0, 40.0, 60.0, 100.0, 150.0};
    int n_masses = sizeof(masses) / sizeof(masses[0]);
    
    /* Header */
    printf("%-8s %-12s %-10s %-10s %-10s %-10s %-10s\n",
           "Mass", "Wind Mass", "C", "N", "O", "Mg", "Fe");
    printf("%-8s %-12s %-10s %-10s %-10s %-10s %-10s\n",
           "(M_sun)", "(M_sun)", "(M_sun)", "(M_sun)", "(M_sun)", "(M_sun)", "(M_sun)");
    printf("--------------------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_masses; i++) {
        SWL_StarWindYield result;
        swl_calculate_star_wind_yield(masses[i], Z, -1, &result);
        
        printf("%-8.1f %-12.4e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e\n",
               masses[i],
               result.total_wind_mass,
               result.element_yields.yields[SWL_EL_C],
               result.element_yields.yields[SWL_EL_N],
               result.element_yields.yields[SWL_EL_O],
               result.element_yields.yields[SWL_EL_MG],
               result.element_yields.yields[SWL_EL_FE]);
    }
}

/*
 * Test 4: Time-dependent cluster wind yields
 */
void test_cluster_time_evolution(SWL_IMFType imf, double Z, double M_cluster) {
    print_separator("TEST 4: Time Evolution of Cluster Wind Yields");
    
    printf("IMF: %s\n", swl_get_imf_name(imf));
    printf("Metallicity: Z = %.4f\n", Z);
    printf("Cluster Mass: %.2e M_sun\n\n", M_cluster);
    
    double ages[] = {1e6, 3e6, 5e6, 1e7, 2e7, 5e7, 1e8, 3e8, 1e9, 3e9, 1e10};
    int n_ages = sizeof(ages) / sizeof(ages[0]);
    
    /* Header */
    printf("%-10s %-10s %-12s %-12s %-12s %-12s %-12s\n",
           "Age", "M_TO", "Wind Mass", "C yield", "N yield", "O yield", "Z_total");
    printf("%-10s %-10s %-12s %-12s %-12s %-12s %-12s\n",
           "(yr)", "(M_sun)", "(M_sun)", "(M_sun)", "(M_sun)", "(M_sun)", "(M_sun)");
    printf("-------------------------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_ages; i++) {
        SWL_ClusterWindYield result;
        swl_calculate_cluster_wind_yield(imf, Z, M_cluster, ages[i], &result);
        
        char age_str[16];
        if (ages[i] >= 1e9) {
            sprintf(age_str, "%.1fGyr", ages[i] / 1e9);
        } else if (ages[i] >= 1e6) {
            sprintf(age_str, "%.0fMyr", ages[i] / 1e6);
        } else {
            sprintf(age_str, "%.0eyr", ages[i]);
        }
        
        printf("%-10s %-10.2f %-12.4e %-12.4e %-12.4e %-12.4e %-12.4e\n",
               age_str,
               result.turnoff_mass,
               result.total_wind_mass,
               result.element_yields.yields[SWL_EL_C],
               result.element_yields.yields[SWL_EL_N],
               result.element_yields.yields[SWL_EL_O],
               result.element_yields.yields[SWL_EL_Z_TOTAL]);
    }
}

/*
 * Test 5: Detailed cluster analysis at specific age
 */
void test_detailed_cluster(SWL_IMFType imf, double Z, double M_cluster, double age) {
    print_separator("TEST 5: Detailed Cluster Analysis");
    
    SWL_ClusterWindYield result;
    SWL_ErrorCode err = swl_calculate_cluster_wind_yield(imf, Z, M_cluster, age, &result);
    
    if (err != SWL_SUCCESS) {
        printf("Error: %s\n", swl_get_error_message(err));
        return;
    }
    
    /* Print input parameters */
    print_subseparator("Input Parameters");
    printf("IMF Model:           %s\n", swl_get_imf_name(imf));
    printf("Metallicity:         Z = %.4e (Z/Z_sun = %.3f)\n", Z, Z / SWL_Z_SUN);
    printf("Cluster Mass:        %.4e M_sun\n", M_cluster);
    printf("Cluster Age:         %.4e years", age);
    if (age >= 1e9) printf(" (%.2f Gyr)", age / 1e9);
    else if (age >= 1e6) printf(" (%.2f Myr)", age / 1e6);
    printf("\n");
    
    /* Print mass budget */
    print_subseparator("Mass Budget");
    printf("Initial cluster mass:    %12.4e M_sun\n", result.cluster_mass);
    printf("Mass in living stars:    %12.4e M_sun (%.1f%%)\n", 
           result.mass_in_living_stars, 100.0 * result.mass_in_living_stars / M_cluster);
    printf("Mass in remnants:        %12.4e M_sun (%.1f%%)\n",
           result.mass_in_remnants, 100.0 * result.mass_in_remnants / M_cluster);
    printf("Mass in winds (to ISM):  %12.4e M_sun (%.1f%%)\n",
           result.total_wind_mass, 100.0 * result.total_wind_mass / M_cluster);
    
    double accounted = result.mass_in_living_stars + result.mass_in_remnants + result.total_wind_mass;
    printf("Total accounted:         %12.4e M_sun (%.1f%%)\n", accounted, 100.0 * accounted / M_cluster);
    
    /* Print stellar population info */
    print_subseparator("Stellar Population");
    printf("Total stars formed:      %12.4e\n", result.n_stars_total);
    printf("Stars that have died:    %12.4e (%.1f%%)\n", 
           result.n_stars_dead, 100.0 * result.n_stars_dead / result.n_stars_total);
    printf("Main sequence turnoff:   %12.3f M_sun\n", result.turnoff_mass);
    
    /* Print current rates */
    print_subseparator("Current Wind Rates (at cluster age)");
    printf("Wind mass-loss rate:     %12.4e M_sun/yr\n", result.wind_mass_rate);
    
    /* Print element yields */
    print_subseparator("Element Yields from Stellar Winds");
    printf("%-10s %15s %15s %15s\n", "Element", "Yield (M_sun)", "Yield/M_cluster", "Rate (M_sun/yr)");
    printf("%-10s %15s %15s %15s\n", "-------", "-------------", "---------------", "---------------");
    
    for (int i = 0; i < SWL_N_ELEMENTS; i++) {
        double y = result.element_yields.yields[i];
        double y_norm = y / M_cluster;
        double rate = result.yield_rate.yields[i];
        
        printf("%-10s %15.6e %15.6e %15.6e\n",
               swl_get_element_name((SWL_Element)i), y, y_norm, rate);
    }
    
    /* Print comparison to solar abundance */
    print_subseparator("Yield Enhancement Relative to Solar");
    printf("(Net yield / (M_cluster * solar_abundance))\n\n");
    printf("%-10s %15s\n", "Element", "Enhancement");
    printf("%-10s %15s\n", "-------", "-----------");
    
    for (int i = SWL_EL_C; i < SWL_N_ELEMENTS; i++) {
        double y = result.element_yields.yields[i];
        double solar = swl_get_solar_abundance((SWL_Element)i);
        double enhancement = y / (M_cluster * solar);
        printf("%-10s %15.4f\n", swl_get_element_name((SWL_Element)i), enhancement);
    }
}

/*
 * Test 6: Compare all IMFs
 */
void test_imf_comparison(double Z, double M_cluster, double age) {
    print_separator("TEST 6: IMF Comparison");
    
    printf("Metallicity: Z = %.4f\n", Z);
    printf("Cluster Mass: %.2e M_sun\n", M_cluster);
    printf("Cluster Age: %.2e yr\n\n", age);
    
    /* Header */
    printf("%-12s", "Quantity");
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        printf(" %18s", swl_get_imf_name((SWL_IMFType)i));
    }
    printf("\n");
    
    printf("%-12s", "--------");
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        printf(" %18s", "------------------");
    }
    printf("\n");
    
    /* Calculate for all IMFs */
    SWL_ClusterWindYield results[SWL_N_IMF_TYPES];
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        swl_calculate_cluster_wind_yield((SWL_IMFType)i, Z, M_cluster, age, &results[i]);
    }
    
    /* Print turnoff mass */
    printf("%-12s", "M_TO");
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        printf(" %18.3f", results[i].turnoff_mass);
    }
    printf("\n");
    
    /* Print wind mass */
    printf("%-12s", "Wind mass");
    for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
        printf(" %18.4e", results[i].total_wind_mass);
    }
    printf("\n");
    
    /* Print selected element yields */
    int elements[] = {SWL_EL_C, SWL_EL_N, SWL_EL_O, SWL_EL_MG, SWL_EL_SI, SWL_EL_FE, SWL_EL_Z_TOTAL};
    int n_el = sizeof(elements) / sizeof(elements[0]);
    
    for (int e = 0; e < n_el; e++) {
        int el = elements[e];
        printf("%-12s", swl_get_element_name((SWL_Element)el));
        for (int i = 0; i < SWL_N_IMF_TYPES; i++) {
            printf(" %18.4e", results[i].element_yields.yields[el]);
        }
        printf("\n");
    }
}

/*
 * Test 7: Metallicity dependence
 */
void test_metallicity_dependence(SWL_IMFType imf, double M_cluster, double age) {
    print_separator("TEST 7: Metallicity Dependence");
    
    printf("IMF: %s\n", swl_get_imf_name(imf));
    printf("Cluster Mass: %.2e M_sun\n", M_cluster);
    printf("Cluster Age: %.2e yr\n\n", age);
    
    double metallicities[] = {0.0001, 0.0003, 0.001, 0.003, 0.008, 0.014, 0.02, 0.03, 0.05};
    int n_Z = sizeof(metallicities) / sizeof(metallicities[0]);
    
    /* Header */
    printf("%-10s %-10s %-12s %-12s %-12s %-12s %-12s\n",
           "Z", "Z/Z_sun", "Wind Mass", "C", "N", "O", "Fe");
    printf("--------------------------------------------------------------------------------\n");
    
    for (int i = 0; i < n_Z; i++) {
        SWL_ClusterWindYield result;
        swl_calculate_cluster_wind_yield(imf, metallicities[i], M_cluster, age, &result);
        
        printf("%-10.4f %-10.3f %-12.4e %-12.4e %-12.4e %-12.4e %-12.4e\n",
               metallicities[i],
               metallicities[i] / SWL_Z_SUN,
               result.total_wind_mass,
               result.element_yields.yields[SWL_EL_C],
               result.element_yields.yields[SWL_EL_N],
               result.element_yields.yields[SWL_EL_O],
               result.element_yields.yields[SWL_EL_FE]);
    }
}

/*
 * Test 8: Export results to CSV
 */
void test_csv_export(SWL_IMFType imf, double Z, double M_cluster, double age) {
    print_separator("TEST 8: Export Results to CSV");
    
    SWL_ClusterWindYield result;
    swl_calculate_cluster_wind_yield(imf, Z, M_cluster, age, &result);
    
    const char* filename = "cluster_wind_yields.csv";
    SWL_ErrorCode err = swl_export_to_csv(&result, filename);
    
    if (err == SWL_SUCCESS) {
        printf("Results exported to: %s\n", filename);
    } else {
        printf("Export failed: %s\n", swl_get_error_message(err));
    }
}

/* ============================================================================
 * MAIN PROGRAM
 * ============================================================================ */

int main(int argc, char* argv[]) {
    /* Default parameters */
    SWL_IMFType imf_type = SWL_IMF_CHABRIER;
    double metallicity = 0.02;    /* Solar */
    double cluster_mass = 1.0e6;  /* 10^6 M_sun */
    double cluster_age = 1.0e7;   /* 10 Myr */
    
    /* Parse command line arguments */
    if (argc >= 2) {
        if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        int imf_val = atoi(argv[1]);
        if (imf_val >= 0 && imf_val < SWL_N_IMF_TYPES) {
            imf_type = (SWL_IMFType)imf_val;
        }
    }
    if (argc >= 3) {
        metallicity = atof(argv[2]);
        if (metallicity <= 0 || metallicity > 0.1) metallicity = 0.02;
    }
    if (argc >= 4) {
        cluster_mass = atof(argv[3]);
        if (cluster_mass <= 0) cluster_mass = 1.0e6;
    }
    if (argc >= 5) {
        cluster_age = atof(argv[4]);
        if (cluster_age < 0) cluster_age = 1.0e7;
    }
    
    /* Print header */
    printf("\n");
    printf("================================================================================\n");
    printf("         STELLAR WIND YIELD LIBRARY - Test and Demonstration Program\n");
    printf("================================================================================\n");
    printf("\n");
    printf("This program demonstrates the stellar wind yield library functions.\n");
    printf("Yield tables based on Kobayashi et al. (2020), Limongi & Chieffi (2018),\n");
    printf("Higgins et al. (2023), and Karakas & Lugaro (2016).\n");
    printf("\n");
    printf("Selected parameters:\n");
    printf("  IMF Model:     %s\n", swl_get_imf_name(imf_type));
    printf("  Metallicity:   Z = %.4f (Z/Z_sun = %.2f)\n", metallicity, metallicity / SWL_Z_SUN);
    printf("  Cluster Mass:  %.2e M_sun\n", cluster_mass);
    printf("  Cluster Age:   %.2e years\n", cluster_age);
    
    /* Initialize library */
    swl_init();
    
    /* Run tests */
    test_basic_functions();
    test_stellar_lifetimes(metallicity);
    test_single_star_yields(metallicity);
    test_cluster_time_evolution(imf_type, metallicity, cluster_mass);
    test_detailed_cluster(imf_type, metallicity, cluster_mass, cluster_age);
    test_imf_comparison(metallicity, cluster_mass, cluster_age);
    test_metallicity_dependence(imf_type, cluster_mass, cluster_age);
    test_csv_export(imf_type, metallicity, cluster_mass, cluster_age);
    
    /* Cleanup */
    swl_cleanup();
    
    printf("\n");
    printf("================================================================================\n");
    printf("                              END OF TESTS\n");
    printf("================================================================================\n\n");
    
    return 0;
}

/*
 * ============================================================================
 * Sample Cooling Table Generator
 * ============================================================================
 * 
 * This utility generates a sample cooling/heating rate table in the binary
 * format used by cooling_rate_mpi.c, for testing the exact integrator.
 *
 * The cooling rates are based on simplified analytical fits to realistic
 * ISM cooling curves.
 *
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_Temperature_Array 128

/* Grid parameters structure - must match cooling_rate_mpi.c exactly */
typedef struct {
    int L, M, N;                        /* 3D grid dimensions */
    double density_min, density_max;
    double metallicity_min, metallicity_max;
    double temp_min, temp_max;
    double temp[N_Temperature_Array];
} GridParams;

/* Physical constants */
#define K_BOLTZMANN 1.3806504e-16
#define M_PROTON    1.6726219e-24

/* Calculate simplified cooling rate [erg/g/s] */
double calc_cooling_rate(double nH, double Z, double T) {
    double rho = 1.22 * M_PROTON * nH;
    double n2_over_rho = nH * nH / rho;
    
    /* Primordial H+He cooling (Cen 1992 style) */
    double Lambda_v = 0.0;
    
    if (T > 1e4) {
        double sqrtT5 = sqrt(T / 1e5);
        /* H excitation cooling */
        Lambda_v += 7.5e-19 * exp(-118348.0/T) / (1.0 + sqrtT5);
        /* He excitation cooling */
        Lambda_v += 5.54e-17 * pow(T, -0.397) * exp(-473638.0/T) / (1.0 + sqrtT5);
        /* Bremsstrahlung */
        Lambda_v += 1.42e-27 * sqrt(T);
    }
    
    /* Metal cooling (simplified Sutherland & Dopita style) */
    if (Z > 0 && T > 100) {
        double Lm = 0.0;
        if (T > 1e4 && T < 1e5) {
            Lm = 1.0e-22 * pow(T/1e4, 0.5);
        } else if (T >= 1e5 && T < 1e7) {
            /* Peak near 2e5 K */
            Lm = 2.5e-22 * exp(-pow(log10(T/2e5), 2) / 0.5);
        } else if (T >= 1e7) {
            Lm = 3.0e-23 * pow(T/1e7, -0.7);
        }
        /* Fine-structure cooling at low T */
        if (T > 100 && T < 1e4) {
            Lm += 1e-26 * (T/100.0) * exp(-90.0/T);
        }
        Lambda_v += Z * Lm;
    }
    
    /* Convert to [erg/g/s] */
    return Lambda_v * n2_over_rho;
}

/* Calculate simplified heating rate [erg/g/s] */
double calc_heating_rate(double nH, double Z, double T) {
    double rho = 1.22 * M_PROTON * nH;
    double n_over_rho = nH / rho;
    
    /* UV background photoheating (HM12 z=0 style) */
    double Gamma_H = 0.0;
    
    double n_crit = 0.01;  /* self-shielding threshold */
    double shield = exp(-nH / n_crit);
    
    if (T < 5e4) {
        /* Photoheating rate scales with recombination rate */
        double alpha_B = 2.6e-13 * pow(T/1e4, -0.7);
        double eps = 5e-12;  /* ~3 eV per ionization */
        Gamma_H = eps * alpha_B * nH * (1.0 - shield);
    }
    
    /* Convert to [erg/g/s] */
    return Gamma_H * n_over_rho;
}

int main(int argc, char *argv[]) {
    const char *outfile = (argc > 1) ? argv[1] : "sample_cooling_table.bin";
    
    printf("================================================================\n");
    printf("  Sample Cooling Table Generator\n");
    printf("================================================================\n\n");
    
    /* Set grid parameters */
    GridParams params;
    params.L = 10;   /* Density points */
    params.M = 10;   /* Metallicity points */
    params.N = N_Temperature_Array;
    
    params.density_min = 1e-8;
    params.density_max = 1e3;
    params.metallicity_min = 1e-6;
    params.metallicity_max = 1.0;
    params.temp_min = 1e1;
    params.temp_max = 1e10;
    
    /* Build temperature grid (log-spaced) */
    double log_T_min = log10(params.temp_min * 1.02);
    double log_T_max = log10(params.temp_max * 0.98);
    double log_T_step = (log_T_max - log_T_min) / (params.N - 1);
    
    for (int k = 0; k < params.N; k++) {
        params.temp[k] = pow(10.0, log_T_min + k * log_T_step);
    }
    
    printf("Grid dimensions: %d x %d x %d\n", params.L, params.M, params.N);
    printf("Density range: [%.2e, %.2e] cm^-3\n", params.density_min, params.density_max);
    printf("Metallicity range: [%.2e, %.2e] Z_solar\n", params.metallicity_min, params.metallicity_max);
    printf("Temperature range: [%.2e, %.2e] K\n\n", params.temp[0], params.temp[params.N-1]);
    
    /* Build grids */
    int total_size = params.L * params.M * params.N;
    double *cooling_grid = (double*)malloc(total_size * sizeof(double));
    double *heating_grid = (double*)malloc(total_size * sizeof(double));
    
    double log_dens_min = log10(params.density_min);
    double log_dens_max = log10(params.density_max);
    double dens_step = (log_dens_max - log_dens_min) / (params.L - 1);
    
    double log_metal_min = log10(params.metallicity_min);
    double log_metal_max = log10(params.metallicity_max);
    double metal_step = (log_metal_max - log_metal_min) / (params.M - 1);
    
    printf("Calculating cooling/heating rates...\n");
    
    for (int i = 0; i < params.L; i++) {
        double nH = pow(10.0, log_dens_min + i * dens_step);
        
        for (int j = 0; j < params.M; j++) {
            double Z = pow(10.0, log_metal_min + j * metal_step);
            
            for (int k = 0; k < params.N; k++) {
                double T = params.temp[k];
                int idx = (i * params.M + j) * params.N + k;
                
                cooling_grid[idx] = calc_cooling_rate(nH, Z, T);
                heating_grid[idx] = calc_heating_rate(nH, Z, T);
            }
        }
    }
    
    /* Write binary file */
    FILE *fp = fopen(outfile, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create file '%s'\n", outfile);
        free(cooling_grid);
        free(heating_grid);
        return 1;
    }
    
    fwrite(&params, sizeof(GridParams), 1, fp);
    fwrite(cooling_grid, sizeof(double), total_size, fp);
    fwrite(heating_grid, sizeof(double), total_size, fp);
    
    fclose(fp);
    
    printf("Written to: %s\n", outfile);
    printf("File size: %.2f KB\n\n", (sizeof(GridParams) + 2 * total_size * sizeof(double)) / 1024.0);
    
    /* Print sample values */
    printf("Sample values at nH=0.1 cm^-3, Z=0.1 Z_solar:\n");
    printf("%12s %16s %16s %16s\n", "T [K]", "Cooling", "Heating", "Net");
    printf("----------------------------------------------------------------\n");
    
    double nH = 0.1, Z = 0.1;
    double temps[] = {1e2, 1e3, 1e4, 3e4, 1e5, 1e6, 1e7};
    for (int i = 0; i < 7; i++) {
        double T = temps[i];
        double cool = calc_cooling_rate(nH, Z, T);
        double heat = calc_heating_rate(nH, Z, T);
        printf("%12.2e %16.4e %16.4e %16.4e\n", T, cool, heat, cool - heat);
    }
    
    free(cooling_grid);
    free(heating_grid);
    
    printf("\n================================================================\n");
    printf("  Done. Use with: ./test_exact_cooling %s\n", outfile);
    printf("================================================================\n");
    
    return 0;
}

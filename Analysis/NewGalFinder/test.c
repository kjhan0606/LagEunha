#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sfftw.h>
#include <srfftw.h>

// define 3D indexing macro
#define INDEX(i, j, k) ((i) + mx * ((j) + ny * (k)))

void gaussian_smoothing(fftw_real *data, int nx, int ny, int nz, double cellsize, double RG) {
    int i, j, k;
    int mx = 2 * (nx/2 + 1); // padded x-dim for real->complex in-place FFT
    long total_size = (long)mx * ny * nz;

    // Create forward and inverse FFT plans
    rfftwnd_plan plan_forward = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
    rfftwnd_plan plan_backward = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

    // Cast to complex pointer (in-place FFT)
    fftw_complex *data_c = (fftw_complex *)data;

    // Forward FFT (real to complex)
    rfftwnd_one_real_to_complex(plan_forward, data, NULL);

    // Prepare Gaussian filter factor
    int nxh = nx/2 + 1;
    int nyh = ny/2;
    int nzh = nz/2;
    double scale = 1.0 / total_size;
    double fac = - (M_PI * RG / cellsize) * (M_PI * RG / cellsize);

    for (k = 0; k < nz; ++k) {
        int kz = (k <= nzh) ? k : k - nz;
        for (j = 0; j < ny; ++j) {
            int ky = (j <= nyh) ? j : j - ny;
            for (i = 0; i < nxh; ++i) {
                int kx = i;
                int idx = i + nxh * (j + ny * k);

                double k2 = kx*kx + ky*ky + kz*kz;
                double gauss = exp(fac * k2);
                data_c[idx].re *= gauss * scale;
                data_c[idx].im *= gauss * scale;
            }
        }
    }

    // Inverse FFT (complex to real)
    rfftwnd_one_complex_to_real(plan_backward, data_c, NULL);

    // Destroy plans
    rfftwnd_destroy_plan(plan_forward);
    rfftwnd_destroy_plan(plan_backward);
}

// Test main function
int main() {
    int nx = 32, ny = 32, nz = 32;
    double cellsize = 1.0;
    double RG = 2.0;

    int mx = 2 * (nx/2 + 1);
    long total_size = (long)mx * ny * nz;

    // Allocate memory
    fftw_real *density = (fftw_real *) malloc(sizeof(fftw_real) * total_size);
    if (!density) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Initialize test data: delta function at center
    for (int i = 0; i < total_size; ++i)
        density[i] = 0.0;
    int cx = nx/2, cy = ny/2, cz = nz/2;
    density[INDEX(cx, cy, cz)] = 1.0;

    // Apply Gaussian smoothing
    gaussian_smoothing(density, nx, ny, nz, cellsize, RG);

    // Print central slice
    printf("Smoothed density at z = %d:\n", nz/2);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            printf("%8.5f ", density[INDEX(i, j, nz/2)]);
        }
        printf("\n");
    }
	FILE *fp = fopen("test.out","w");
	fwrite(&mx, sizeof(int), 1,fp);
	fwrite(&ny, sizeof(int), 1,fp);
	fwrite(&nz, sizeof(int), 1,fp);
	fwrite(density, sizeof(float), mx*ny*nz,fp);
	fclose(fp);

    free(density);
    return 0;
}


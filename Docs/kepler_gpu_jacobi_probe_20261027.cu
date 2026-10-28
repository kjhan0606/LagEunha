/* Isolated CUDA probe for J Wdot = -Vdot_fixed on an existing power FaceCSR.
 * This deliberately does not alter particle weights or launch the hydro force.
 * Boundary and MPI padding faces need an explicit treatment before production.
 */
#include <cuda_runtime.h>
#include <math.h>
#include "exam_gpu.h"

__global__ static void volume_rate_rows(FaceCSR faces, ParticleSoA parts,
    const unsigned char *active, int np, int nt, double *conductance,
    double *rhs, double *diagonal, unsigned long long *unsupported)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= np) return;
    if (!active[i]) {
        rhs[i] = 0.; diagonal[i] = 1.;
        return;
    }
    const double xi = parts.x[i], yi = parts.y[i];
    const double vxi = parts.vx[i], vyi = parts.vy[i];
    double area_rate = 0., diag = 0.;
    for (int f = faces.face_offset[i]; f < faces.face_offset[i+1]; ++f) {
        const int j = faces.neighbor_idx[f];
        if (j < 0 || j >= nt || faces.is_ghost[f]) {
            atomicAdd(unsupported, 1ull);
            conductance[f] = 0.;
            continue;
        }
        const double dx = parts.x[j]-xi, dy = parts.y[j]-yi;
        const double d = sqrt(dx*dx+dy*dy);
        const double ex = faces.c2x[f]-faces.c1x[f];
        const double ey = faces.c2y[f]-faces.c1y[f];
        const double length = sqrt(ex*ex+ey*ey);
        if (!(d > 0.) || !(length > 0.)) {
            atomicAdd(unsupported, 1ull);
            conductance[f] = 0.;
            continue;
        }
        const double cx = xi + .5*(faces.c1x[f]+faces.c2x[f]);
        const double cy = yi + .5*(faces.c1y[f]+faces.c2y[f]);
        const double dvx = parts.vx[j]-vxi;
        const double dvy = parts.vy[j]-vyi;
        const double un = ((parts.x[j]*parts.vx[j] + parts.y[j]*parts.vy[j])
                         -(xi*vxi+yi*vyi) - dvx*cx-dvy*cy)/d;
        const double g = length/(2.*d);
        conductance[f] = g;
        diag += g;
        area_rate += length*un;
    }
    if (!(diag > 0.)) atomicAdd(unsupported, 1ull);
    rhs[i] = -area_rate;
    diagonal[i] = diag > 0. ? diag : 1.;
}

__global__ static void jacobi_rate_rows(FaceCSR faces,
    const unsigned char *active, int np, int nt,
    const double *conductance, const double *rhs, const double *diagonal,
    const double *previous, double *next, double omega)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= np) return;
    if (!active[i]) { next[i] = 0.; return; }
    double jx = 0.;
    const double xi = previous[i];
    for (int f = faces.face_offset[i]; f < faces.face_offset[i+1]; ++f) {
        const int j = faces.neighbor_idx[f];
        const double xj = (j >= 0 && j < nt && active[j] && j < np)
                        ? previous[j] : 0.;
        jx += conductance[f]*(xi-xj);
    }
    next[i] = xi + omega*(rhs[i]-jx)/diagonal[i];
}

extern "C" int gpu_probe_volume_rate_jacobi(GPUContext *ctx,
    const unsigned char *host_active, int n_particles, int n_total,
    int iterations, double omega, double *host_rate,
    unsigned long long *host_unsupported)
{
    if (!ctx || !ctx->initialized || ctx->d_faces.n_faces_total <= 0 ||
        !host_active || !host_rate || !host_unsupported ||
        n_particles <= 0 || n_total < n_particles ||
        n_particles != ctx->d_faces.n_particles ||
        iterations <= 0 || !(omega > 0. && omega <= 1.)) return -1;

    cudaStream_t stream = (cudaStream_t)ctx->stream;
    unsigned char *d_active = 0;
    double *d_g = 0, *d_rhs = 0, *d_diag = 0, *d_a = 0, *d_b = 0;
    unsigned long long *d_unsupported = 0;
    int status = 1;
    cudaError_t err;
    const int block = 256;
    const int grid = (n_particles+block-1)/block;
#define CUDA_TRY(call) do { err = (call); if (err != cudaSuccess) goto done; } while (0)
    CUDA_TRY(cudaMalloc((void**)&d_active, (size_t)n_total*sizeof(*d_active)));
    CUDA_TRY(cudaMalloc((void**)&d_g, (size_t)ctx->d_faces.n_faces_total*sizeof(*d_g)));
    CUDA_TRY(cudaMalloc((void**)&d_rhs, (size_t)n_particles*sizeof(*d_rhs)));
    CUDA_TRY(cudaMalloc((void**)&d_diag, (size_t)n_particles*sizeof(*d_diag)));
    CUDA_TRY(cudaMalloc((void**)&d_a, (size_t)n_particles*sizeof(*d_a)));
    CUDA_TRY(cudaMalloc((void**)&d_b, (size_t)n_particles*sizeof(*d_b)));
    CUDA_TRY(cudaMalloc((void**)&d_unsupported, sizeof(*d_unsupported)));
    CUDA_TRY(cudaMemcpyAsync(d_active, host_active, (size_t)n_total,
                             cudaMemcpyHostToDevice, stream));
    CUDA_TRY(cudaMemsetAsync(d_a, 0, (size_t)n_particles*sizeof(*d_a), stream));
    CUDA_TRY(cudaMemsetAsync(d_unsupported, 0, sizeof(*d_unsupported), stream));
    volume_rate_rows<<<grid,block,0,stream>>>(ctx->d_faces, ctx->d_parts,
        d_active, n_particles, n_total, d_g, d_rhs, d_diag, d_unsupported);
    CUDA_TRY(cudaGetLastError());
    for (int k = 0; k < iterations; ++k) {
        jacobi_rate_rows<<<grid,block,0,stream>>>(ctx->d_faces, d_active,
            n_particles, n_total, d_g, d_rhs, d_diag, d_a, d_b, omega);
        CUDA_TRY(cudaGetLastError());
        double *tmp = d_a; d_a = d_b; d_b = tmp;
    }
    CUDA_TRY(cudaMemcpyAsync(host_rate, d_a, (size_t)n_particles*sizeof(*host_rate),
                             cudaMemcpyDeviceToHost, stream));
    CUDA_TRY(cudaMemcpyAsync(host_unsupported, d_unsupported,
                             sizeof(*host_unsupported), cudaMemcpyDeviceToHost, stream));
    CUDA_TRY(cudaStreamSynchronize(stream));
    status = *host_unsupported ? -2 : 0;
done:
    cudaFree(d_unsupported);
    cudaFree(d_b); cudaFree(d_a); cudaFree(d_diag); cudaFree(d_rhs);
    cudaFree(d_g); cudaFree(d_active);
#undef CUDA_TRY
    return status;
}

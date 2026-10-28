/* Experimental GPU FaceCSR pressure-work gradient and host scalar projection.
 * Reuses the independently validated Jacobi16 kernel in this translation unit.
 * Does not change the force, RK4 state, or particle weights.
 */
#include "kepler_gpu_jacobi_probe_20261027.cu"
#include <float.h>

__global__ static void pressure_weight_gradient_rows(
    FaceCSR faces, ParticleSoA parts, const unsigned char *active,
    int np, int nt, double *gradient)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=np) return;
    if(!active[i]) { gradient[i]=0.; return; }
    const double xi=parts.x[i], yi=parts.y[i];
    const double pi=(double)parts.pressure[i];
    double h=0.;
    for(int f=faces.face_offset[i]; f<faces.face_offset[i+1]; ++f) {
        int j=faces.neighbor_idx[f];
        if(j<0 || j>=nt || faces.is_ghost[f]) { gradient[i]=NAN; return; }
        double dx=parts.x[j]-xi, dy=parts.y[j]-yi;
        double d=sqrt(dx*dx+dy*dy);
        double ex=faces.c2x[f]-faces.c1x[f];
        double ey=faces.c2y[f]-faces.c1y[f];
        double length=sqrt(ex*ex+ey*ey);
        if(!(d>0.) || !(length>0.)) { gradient[i]=NAN; return; }
        h += (length/(2*d))*(pi-(double)parts.pressure[j]);
    }
    gradient[i]=h;
}

extern "C" int gpu_probe_energy_neutral_rate(GPUContext *ctx,
    const unsigned char *host_active, int n_particles, int n_total,
    int iterations, double omega, double *host_rate,
    double *host_gradient, double *host_scalar,
    unsigned long long *host_unsupported)
{
    if(!ctx || !host_gradient || !host_scalar || !ctx->d_parts.pressure)
        return -1;
    int status=gpu_probe_volume_rate_jacobi(ctx, host_active,
        n_particles, n_total, iterations, omega, host_rate,
        host_unsupported);
    if(status) return status;
    cudaStream_t stream=(cudaStream_t)ctx->stream;
    unsigned char *d_active=0;
    double *d_h=0;
    cudaError_t error;
    status=1;
    error=cudaMalloc((void**)&d_active,(size_t)n_total);
    if(error!=cudaSuccess) goto done;
    error=cudaMalloc((void**)&d_h,(size_t)n_particles*sizeof(double));
    if(error!=cudaSuccess) goto done;
    error=cudaMemcpyAsync(d_active,host_active,(size_t)n_total,
                          cudaMemcpyHostToDevice,stream);
    if(error!=cudaSuccess) goto done;
    pressure_weight_gradient_rows<<<(n_particles+255)/256,256,0,stream>>>(
        ctx->d_faces,ctx->d_parts,d_active,n_particles,n_total,d_h);
    error=cudaGetLastError();
    if(error!=cudaSuccess) goto done;
    error=cudaMemcpyAsync(host_gradient,d_h,
        (size_t)n_particles*sizeof(double),cudaMemcpyDeviceToHost,stream);
    if(error!=cudaSuccess) goto done;
    error=cudaStreamSynchronize(stream);
    if(error!=cudaSuccess) goto done;
    {
        long double dot=0.,norm=0.;
        for(int i=0;i<n_particles;++i) {
            if(!host_active[i]) continue;
            if(!isfinite(host_gradient[i]) || !isfinite(host_rate[i])) {
                status=-2; goto done;
            }
            dot+=(long double)host_gradient[i]*host_rate[i];
            norm+=(long double)host_gradient[i]*host_gradient[i];
        }
        if(!(norm>0.)) { status=-3; goto done; }
        double alpha=(double)(dot/norm);
        for(int i=0;i<n_particles;++i) {
            if(host_active[i]) host_rate[i]-=alpha*host_gradient[i];
        }
        *host_scalar=alpha;
    }
    status=0;
done:
    cudaFree(d_h);
    cudaFree(d_active);
    return status;
}

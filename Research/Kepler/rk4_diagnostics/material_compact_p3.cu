#include "material_compact_p3_cuda.h"
#include <cfloat>
namespace {
__global__ void sample(const compact_p3::Triangle*triangles,size_t nt,
    const STCoefficient*coefficient,const compact_p3::ThermalDensity*fields,
    const compact_p3::Sample*points,size_t n,compact_p3::SampleValue*out){
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        auto&value=out[i];value={};const auto&q=points[i];
        if(q.triangle>=nt)continue;
        double sum=0;bool inside=true;
        for(int k=0;k<3;k++){inside=inside&&isfinite(q.l[k])&&q.l[k]>=0;sum+=q.l[k];}
        if(!inside||fabs(sum-1)>32*DBL_EPSILON)continue;
        if(!compact_p3::point(triangles[q.triangle],coefficient,q.l,&value.map))continue;
        double B[10];compact_p3::bernstein(q.l,B);
        value.rho_reference=compact_p3::field(fields[q.triangle].mass,B);
        value.energy_reference=compact_p3::field(fields[q.triangle].internal_energy,B);
        value.valid=isfinite(value.rho_reference)&&value.rho_reference>0&&isfinite(value.energy_reference)&&value.energy_reference>0;
        for(int k=0;k<2;k++)value.valid=value.valid&&isfinite(value.map.position[k])&&isfinite(value.map.velocity[k]);
        for(int k=0;k<4;k++)value.valid=value.valid&&isfinite(value.map.F[k])&&isfinite(value.map.grad_v[k]);
    }
}
}
extern "C" int compact_p3_sample_device(const compact_p3::Triangle*t,size_t nt,
    const STCoefficient*c,const compact_p3::ThermalDensity*f,
    const compact_p3::Sample*q,size_t n,compact_p3::SampleValue*out,cudaStream_t stream){
    if(!n)return int(cudaSuccess);
    if(!t||!c||!f||!q||!out||!nt)return int(cudaErrorInvalidValue);
    size_t wanted=(n-1)/128+1;unsigned blocks=unsigned(wanted>65535?65535:wanted);
    sample<<<blocks,128,0,stream>>>(t,nt,c,f,q,n,out);
    return int(cudaPeekAtLastError());
}

#include "material_compact_geometry_cuda.h"
namespace {
__device__ bool evaluate(const compact_geometry::Patch*p,size_t np,const STCoefficient*c,size_t nc,
    const compact_geometry::Sample&q,compact_geometry::Basis*b,GfsMaterialPoint12*x){
    if(q.patch>=np)return false;double total=0;
    for(double l:q.l){if(!(l>=0)||!isfinite(l))return false;total+=l;}
    if(fabs(total-1)>32*DBL_EPSILON||!compact_geometry::basis(p[q.patch],q.l,b))return false;
    for(int a=0;a<b->n;a++)if(b->id[a]>=nc)return false;
    return compact_geometry::point(*b,c,x);
}
__global__ void sample(const compact_geometry::Patch*p,size_t np,const STCoefficient*c,size_t nc,
    const compact_p3::ThermalDensity*f,const compact_geometry::Sample*q,size_t n,compact_geometry::SampleValue*out){
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        auto&v=out[i];v={};compact_geometry::Basis b;
        if(!evaluate(p,np,c,nc,q[i],&b,&v.map))continue;
        double B[10];compact_p3::bernstein(q[i].l,B);
        v.rho_reference=compact_p3::field(f[q[i].patch].mass,B);
        v.energy_reference=compact_p3::field(f[q[i].patch].internal_energy,B);
        v.valid=v.rho_reference>0&&isfinite(v.rho_reference)&&v.energy_reference>0&&isfinite(v.energy_reference);
    }
}
__global__ void stress(const compact_geometry::Patch*p,size_t np,const STCoefficient*c,size_t nc,
    const compact_geometry::StressInput*q,size_t n,compact_geometry::StressValue*out){
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        out[i]={};compact_geometry::Basis b;GfsMaterialPoint12 x;
        if(evaluate(p,np,c,nc,q[i].point,&b,&x))out[i]=compact_geometry::stress(b,c,x,q[i].rho_reference,q[i].weight,q[i].sigma,q[i].gravity[0],q[i].gravity[1]);
    }
}
__global__ void exterior(const compact_geometry::Patch*p,size_t np,const STCoefficient*c,size_t nc,
    const compact_geometry::BoundaryInput*q,size_t n,compact_geometry::BoundaryValue*out){
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        out[i]={};compact_geometry::Basis b;GfsMaterialPoint12 x;
        if(evaluate(p,np,c,nc,q[i].point,&b,&x))out[i]=compact_geometry::exterior(b,c,x,q[i].normal[0],q[i].normal[1],q[i].weight,q[i].pressure);
    }
}
unsigned blocks(size_t n){size_t wanted=(n-1)/128+1;return unsigned(wanted>65535?65535:wanted);}
}
extern "C" int compact_geometry_sample_device(const compact_geometry::Patch*p,size_t np,const STCoefficient*c,size_t nc,
    const compact_p3::ThermalDensity*f,const compact_geometry::Sample*q,size_t n,compact_geometry::SampleValue*out,cudaStream_t stream){
    if(!n)return int(cudaSuccess);if(!p||!np||!c||!nc||!f||!q||!out)return int(cudaErrorInvalidValue);
    sample<<<blocks(n),128,0,stream>>>(p,np,c,nc,f,q,n,out);return int(cudaPeekAtLastError());
}
extern "C" int compact_geometry_stress_device(const compact_geometry::Patch*p,size_t np,const STCoefficient*c,size_t nc,
    const compact_geometry::StressInput*q,size_t n,compact_geometry::StressValue*out,cudaStream_t stream){
    if(!n)return int(cudaSuccess);if(!p||!np||!c||!nc||!q||!out)return int(cudaErrorInvalidValue);
    stress<<<blocks(n),128,0,stream>>>(p,np,c,nc,q,n,out);return int(cudaPeekAtLastError());
}
extern "C" int compact_geometry_exterior_device(const compact_geometry::Patch*p,size_t np,const STCoefficient*c,size_t nc,
    const compact_geometry::BoundaryInput*q,size_t n,compact_geometry::BoundaryValue*out,cudaStream_t stream){
    if(!n)return int(cudaSuccess);if(!p||!np||!c||!nc||!q||!out)return int(cudaErrorInvalidValue);
    exterior<<<blocks(n),128,0,stream>>>(p,np,c,nc,q,n,out);return int(cudaPeekAtLastError());
}

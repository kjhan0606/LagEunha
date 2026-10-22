#include "check_compact_relative_probe.h"
__global__ void probe(const compact_geometry::Patch*p,const RelativeInput*input,size_t n,RelativeOutput*out){
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n)return;const auto&q=input[i];auto&v=out[i];v={};
    compact_geometry::Basis basis;
    v.valid=compact_geometry::relative_trace(p->trace,STJet(q.t,1,0),q.lo,q.hi,q.B,q.C,v.relative)&&compact_geometry::basis(*p,q.l,&basis);
    if(v.valid)for(int k=0;k<2;k++)v.mode[k]=basis.N[10+k];
}
extern "C" int compact_relative_probe(const compact_geometry::Patch*p,const RelativeInput*q,size_t n,RelativeOutput*out){
    if(!n||!p||!q||!out)return int(cudaErrorInvalidValue);probe<<<unsigned((n+127)/128),128>>>(p,q,n,out);return int(cudaGetLastError());
}

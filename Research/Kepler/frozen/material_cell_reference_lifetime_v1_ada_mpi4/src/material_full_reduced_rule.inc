#include "material_full_reduced_rule.h"
#include <climits>
#include <cmath>
extern "C" int mf_evolved_rule_upload(void*ptr,const uint64_t*offset,size_t nt,const uint32_t*index,const double*weight,size_t nodes){
    auto*s=static_cast<Context*>(ptr);if(!s||s->reduced_reference_only||s->reduced_offset||nt!=s->nt||offset[0]!=0||offset[nt]!=nodes)return cudaErrorInvalidValue;
    for(size_t i=0;i<nt;i++)if(offset[i+1]<offset[i]||offset[i+1]-offset[i]>INT_MAX)return cudaErrorInvalidValue;
    for(size_t i=0;i<nodes;i++)if(!(weight[i]>0)||!std::isfinite(weight[i]))return cudaErrorInvalidValue;
    cudaError_t e=cudaMalloc((void**)&s->reduced_offset,(nt+1)*sizeof(uint64_t));
    if(e==cudaSuccess)e=cudaMemcpy(s->reduced_offset,offset,(nt+1)*sizeof(uint64_t),cudaMemcpyHostToDevice);
    if(e==cudaSuccess&&nodes)e=cudaMalloc((void**)&s->reduced_index,nodes*sizeof(uint32_t));
    if(e==cudaSuccess&&nodes)e=cudaMalloc((void**)&s->reduced_weight,nodes*sizeof(double));
    if(e==cudaSuccess&&nodes)e=cudaMemcpy(s->reduced_index,index,nodes*sizeof(uint32_t),cudaMemcpyHostToDevice);
    if(e==cudaSuccess&&nodes)e=cudaMemcpy(s->reduced_weight,weight,nodes*sizeof(double),cudaMemcpyHostToDevice);
    return e;
}
extern "C" int mf_evolved_rule_use(void*ptr,int use){auto*s=static_cast<Context*>(ptr);if(!s||use<0||use>1||(s->reduced_reference_only&&!use)||(use&&!s->reduced_offset))return cudaErrorInvalidValue;s->reduced_active=bool(use);s->cache_start=s->cache_count=0;return 0;}

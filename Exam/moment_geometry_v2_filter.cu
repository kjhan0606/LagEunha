/* Same scalar interval predicates on CPU/OpenMP/CUDA. This file intentionally
 * builds without fast-math. It evaluates queries, not a GPU topology builder. */
#include "moment_geometry_interval.h"
#include <cuda_runtime.h>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <cfenv>
#include <omp.h>
template<class T>struct MG2Device{
 T*p=nullptr;bool set(size_t n,const T*a=nullptr){if(!n)return true;if(cudaMalloc(&p,n*sizeof(T))!=cudaSuccess)return false;return !a||cudaMemcpy(p,a,n*sizeof(T),cudaMemcpyHostToDevice)==cudaSuccess;}
 ~MG2Device(){if(p)cudaFree(p);}
};
MG2_HD int8_t support_query(const int64_t*q,const double*xy,const double*W,const double*b){
 return mg2_support_sign(mg2_plane(q[0],q[1],xy,W,b),mg2_plane(q[0],q[2],xy,W,b),mg2_plane(q[0],q[3],xy,W,b));
}
MG2_HD int8_t feasibility_query(int64_t owner,const double*enclosure,int j,const double*xy,const double*W,const double*b){
 auto p=mg2_plane(owner,j,xy,W,b);MG2Interval x{enclosure[0],enclosure[1]},y{enclosure[2],enclosure[3]};
 return mg2_sign(mg2_sub(mg2_add(mg2_mul(p.x,x),mg2_mul(p.y,y)),p.b));
}
__global__ void support_queries(int n,const int64_t*q,const double*xy,const double*W,const double*b,int8_t*out){
 int k=blockIdx.x*blockDim.x+threadIdx.x;if(k<n)out[k]=support_query(q+4*k,xy,W,b);
}
__global__ void feasibility_queries(int64_t n,int nr,const int64_t*owners,const double*enclosures,const double*xy,const double*W,const double*b,int8_t*out){
 int64_t k=int64_t(blockIdx.x)*blockDim.x+threadIdx.x;if(k<n)out[k]=feasibility_query(owners[k/nr],enclosures+4*(k/nr),k%nr,xy,W,b);
}
/* mode0: n support queries[4], mode1: n vertices x registry Cartesian product.
 * In mode1 idx[n] is owner and enc[n,4] is outward exact local vertex enclosure.
 * Output sign0 is unresolved; never treat it as an inside/outside decision. */
extern "C" int moment_geometry_v2_filter(int backend,int threads,int mode,int nr,
 const double*xy,const double*W,const double*b,int n,const int64_t*idx,const double*enc,int8_t*output){
 if((backend!=0&&backend!=1)||threads<1||(mode!=0&&mode!=1)||nr<1||n<1||!xy||!W||!b||!idx||!output||(mode==1&&!enc)||std::fegetround()!=FE_TONEAREST)return 1;
 for(int j=0;j<nr;j++)if(!isfinite(xy[2*j])||!isfinite(xy[2*j+1])||!isfinite(W[j]))return 1;
 for(int k=0;k<4;k++)if(!isfinite(b[k]))return 1;
 if(mode==0){for(int k=0;k<n;k++){if(idx[4*k]<0||idx[4*k]>=nr)return 1;for(int t=1;t<4;t++)if(idx[4*k+t]<-4||idx[4*k+t]>=nr||idx[4*k+t]==idx[4*k])return 1;}}
 else for(int k=0;k<n;k++)if(idx[k]<0||idx[k]>=nr||!isfinite(enc[4*k])||!isfinite(enc[4*k+1])||!isfinite(enc[4*k+2])||!isfinite(enc[4*k+3])||enc[4*k]>enc[4*k+1]||enc[4*k+2]>enc[4*k+3])return 1;
 try{
  int64_t count=mode==0?n:int64_t(n)*nr;std::vector<int8_t>out(count);
  if(backend==0){
   #pragma omp parallel for num_threads(threads)
   for(int64_t k=0;k<count;k++)out[k]=mode==0?support_query(idx+4*k,xy,W,b):feasibility_query(idx[k/nr],enc+4*(k/nr),k%nr,xy,W,b);
  }else{
   MG2Device<double>dxy,dw,db,de;MG2Device<int64_t>di;MG2Device<int8_t>dout;
   if(!dxy.set(2*size_t(nr),xy)||!dw.set(nr,W)||!db.set(4,b)||!di.set(size_t(n)*(mode==0?4:1),idx)||!dout.set(count)||(mode==1&&!de.set(4*size_t(n),enc)))return 8;
   if(mode==0)support_queries<<<(count+255)/256,256>>>(n,di.p,dxy.p,dw.p,db.p,dout.p);
   else feasibility_queries<<<(count+255)/256,256>>>(count,nr,di.p,de.p,dxy.p,dw.p,db.p,dout.p);
   if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 7;
   if(cudaMemcpy(out.data(),dout.p,count*sizeof(int8_t),cudaMemcpyDeviceToHost)!=cudaSuccess)return 7;
  }
  std::copy(out.begin(),out.end(),output);return 0;
 }catch(...){return 8;}
}

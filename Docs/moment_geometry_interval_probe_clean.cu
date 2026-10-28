/* Intermediate intervals and nextafter bits for one stored predicate replay. */
#include "moment_geometry_interval.h"
#include <cuda_runtime.h>
MG2_HD void interval_probe_clean(double one,const int64_t*q,const double*xy,const double*W,const double*b,double*out){
 out[0]=::nextafter(one,HUGE_VAL);out[1]=::nextafter(one,-HUGE_VAL);
 out[2]=::nextafter(one,HUGE_VAL);out[3]=::nextafter(one,-HUGE_VAL);
 auto add=mg2_add(mg2_point(one),mg2_point(one)),mul=mg2_mul(mg2_point(one),mg2_point(one));
 out[4]=add.lo;out[5]=add.hi;out[6]=mul.lo;out[7]=mul.hi;
 auto p=mg2_plane(q[0],q[1],xy,W,b),r=mg2_plane(q[0],q[2],xy,W,b),s=mg2_plane(q[0],q[3],xy,W,b);
 MG2Interval fields[9]={p.x,p.y,p.b,r.x,r.y,r.b,s.x,s.y,s.b};
 for(int i=0;i<9;i++){out[8+2*i]=fields[i].lo;out[9+2*i]=fields[i].hi;}
 auto D=mg2_sub(mg2_mul(p.x,r.y),mg2_mul(p.y,r.x));
 auto X=mg2_sub(mg2_mul(p.b,r.y),mg2_mul(p.y,r.b));
 auto Y=mg2_sub(mg2_mul(p.x,r.b),mg2_mul(p.b,r.x));
 auto N=mg2_sub(mg2_add(mg2_mul(s.x,X),mg2_mul(s.y,Y)),mg2_mul(s.b,D));
 fields[0]=D;fields[1]=X;fields[2]=Y;fields[3]=N;
 for(int i=0;i<4;i++){out[26+2*i]=fields[i].lo;out[27+2*i]=fields[i].hi;}
 out[34]=mg2_sign(D);out[35]=mg2_sign(N);out[36]=mg2_support_sign(p,r,s);
}
__global__ void probe_kernel_clean(double one,const int64_t*q,const double*xy,const double*W,const double*b,int n,double*out){int i=threadIdx.x;if(i<n)interval_probe_clean(one,q+4*i,xy,W,b,out+37*i);}
extern "C" int moment_geometry_interval_probe(int backend,int nr,const double*xy,const double*W,const double*b,int n,const int64_t*q,double*out){
 if(backend==0){for(int i=0;i<n;i++)interval_probe_clean(1.,q+4*i,xy,W,b,out+37*i);return 0;}
 double *dx=nullptr,*dw=nullptr,*db=nullptr,*dout=nullptr;int64_t*dq=nullptr;int code=0;
 if(cudaMalloc(&dx,2*nr*sizeof(double))!=cudaSuccess||cudaMalloc(&dw,nr*sizeof(double))!=cudaSuccess||cudaMalloc(&db,4*sizeof(double))!=cudaSuccess||cudaMalloc(&dq,4*n*sizeof(int64_t))!=cudaSuccess||cudaMalloc(&dout,37*n*sizeof(double))!=cudaSuccess)code=8;
 if(!code&&(cudaMemcpy(dx,xy,2*nr*sizeof(double),cudaMemcpyHostToDevice)!=cudaSuccess||cudaMemcpy(dw,W,nr*sizeof(double),cudaMemcpyHostToDevice)!=cudaSuccess||cudaMemcpy(db,b,4*sizeof(double),cudaMemcpyHostToDevice)!=cudaSuccess||cudaMemcpy(dq,q,4*n*sizeof(int64_t),cudaMemcpyHostToDevice)!=cudaSuccess))code=7;
 if(!code){probe_kernel_clean<<<1,128>>>(1.,dq,dx,dw,db,n,dout);if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess||cudaMemcpy(out,dout,37*n*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)code=7;}
 if(dx)cudaFree(dx);if(dw)cudaFree(dw);if(db)cudaFree(db);if(dq)cudaFree(dq);if(dout)cudaFree(dout);return code;
}

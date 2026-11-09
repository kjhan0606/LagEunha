#include <cooperative_groups.h>
namespace material_cg = cooperative_groups;
/* All blocks are cooperatively resident. Every grid barrier is unconditional
 * across participating threads; scalar decisions are broadcast before branching. */
static __device__ void weight_grid_sum(double local,double*part,double*s,int slot,double*work){
 double v=block_sum(local,work);if(threadIdx.x==0)part[blockIdx.x]=v;
 material_cg::this_grid().sync();
 if(blockIdx.x==0){v=0;for(int k=threadIdx.x;k<gridDim.x;k+=blockDim.x)v+=part[k];v=block_sum(v,work);if(threadIdx.x==0)s[slot]=v;}
 material_cg::this_grid().sync();
}
template<bool anchored>
static __global__ void weight_persistent(int n,const int*off,const int*col,const double*a,const double*diag,double*q,double*r,double*p,double*ap,double*s,double*part,double norm,int limit){
 __shared__ double work[256];int tid=blockIdx.x*blockDim.x+threadIdx.x,stride=gridDim.x*blockDim.x;
 if(tid==0){s[5]=1;s[6]=0;s[7]=0;}
 double zsum=0;
 if(!anchored){for(int i=tid;i<n;i+=stride)zsum+=r[i]/diag[i];weight_grid_sum(zsum,part,s,3,work);}
 if(tid==0)s[3]=anchored?0:s[3]/n;material_cg::this_grid().sync();
 double rz=0;for(int i=tid;i<n;i+=stride){p[i]=r[i]/diag[i]-s[3];rz+=r[i]*p[i];}
 weight_grid_sum(rz,part,s,0,work);
 if(tid==0&&(!(s[0]>0)||!isfinite(s[0])))s[5]=0;
 material_cg::this_grid().sync();if(!s[5])return;
 for(int it=0;it<limit;it++){
  double pap=0;
  for(int i=tid;i<n;i+=stride){double v=0;for(int k=off[i];k<off[i+1];k++)v+=a[k]*(p[i]-(anchored&&col[k]<0?0:p[col[k]]));ap[i]=v;pap+=p[i]*v;}
  weight_grid_sum(pap,part,s,1,work);
  if(tid==0){double curvature=s[1];s[1]=s[0]/curvature;if(!(curvature>0)||!isfinite(curvature)||!isfinite(s[1]))s[5]=0;}
  material_cg::this_grid().sync();if(!s[5])return;
  double rr=0;zsum=0;
  for(int i=tid;i<n;i+=stride){q[i]+=s[1]*p[i];r[i]-=s[1]*ap[i];rr+=r[i]*r[i];if(!anchored)zsum+=r[i]/diag[i];}
  weight_grid_sum(rr,part,s,2,work);if(!anchored)weight_grid_sum(zsum,part,s,3,work);
  if(tid==0){s[3]=anchored?0:s[3]/n;s[7]=it+1;if(!isfinite(s[2]))s[5]=0;else if(s[2]<=1e-20*norm)s[6]=1;}
  material_cg::this_grid().sync();if(!s[5]||s[6])return;
  double next=0;for(int i=tid;i<n;i+=stride)next+=r[i]*(r[i]/diag[i]-s[3]);
  weight_grid_sum(next,part,s,4,work);
  if(tid==0){double nrz=s[4];s[4]=nrz/s[0];s[0]=nrz;if(!(nrz>0)||!isfinite(nrz)||!isfinite(s[4]))s[5]=0;}
  material_cg::this_grid().sync();if(!s[5])return;
  for(int i=tid;i<n;i+=stride)p[i]=r[i]/diag[i]-s[3]+s[4]*p[i];
  material_cg::this_grid().sync();
 }
}

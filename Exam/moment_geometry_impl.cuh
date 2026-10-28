#ifndef EXAM_MOMENT_GEOMETRY_IMPL_CUH
#define EXAM_MOMENT_GEOMETRY_IMPL_CUH
#include "moment_geometry.h"
#include "power_polygon_clip.h"
#include <cuda_runtime.h>
#include <vector>
#include <algorithm>
#include <new>
#include <omp.h>
#include <float.h>

PPC_HD int moment_geometry_cell(int row,const double *xy,const double *w,
 const double *v,const double *z,const int64_t *ids,const int *kind,
 const int64_t *selected,const int64_t *off,const int64_t *cand,
 const double *bounds,uint64_t epoch,MomentGeometryPolygon *out){
 const int cap=MOMENT_GEOMETRY_CORNERS;int64_t i=selected[row];
 double xi=xy[2*i],yi=xy[2*i+1];
 double x[cap],y[cap],tx[cap],ty[cap];int64_t tag[cap],tt[cap];
 x[0]=bounds[0]-xi;y[0]=bounds[2]-yi;tag[0]=-1;
 x[1]=bounds[1]-xi;y[1]=bounds[2]-yi;tag[1]=-2;
 x[2]=bounds[1]-xi;y[2]=bounds[3]-yi;tag[2]=-3;
 x[3]=bounds[0]-xi;y[3]=bounds[3]-yi;tag[3]=-4;
 int nc=4;
 for(int64_t k=off[row];k<off[row+1];k++){
  int64_t j=cand[k];double dx=xy[2*j]-xi,dy=xy[2*j+1]-yi;
  double d2=dx*dx+dy*dy,threshold=.5*(d2+w[i]-w[j]);
  if(!isfinite(d2)||d2<=0.||!isfinite(threshold))return 6;
  int next=0,err=power_polygon_clip(nc,x,y,tag,dx,dy,threshold,j,cap,tx,ty,tt,1,0.,&next);
  if(err)return err==1?2:3;
  nc=next;if(nc<3)return 4;
  for(int e=0;e<nc;e++){x[e]=tx[e];y[e]=ty[e];tag[e]=tt[e];}
 }
 double twice=0.,cx=0.,cy=0.,xx=0.,yy=0.,xyint=0.;
 for(int e=0;e<nc;e++){
  int en=(e+1)%nc;double ax=x[e],ay=y[e],bx=x[en],by=y[en],c=ax*by-bx*ay;
  twice+=c;cx+=(ax+bx)*c;cy+=(ay+by)*c;
  xx+=(ax*ax+ax*bx+bx*bx)*c;yy+=(ay*ay+ay*by+by*by)*c;
  xyint+=(2.*ax*ay+ax*by+bx*ay+2.*bx*by)*c;
 }
 if(!isfinite(twice)||twice<=0.)return 6;
 out->epoch=epoch;out->count=nc;out->area=.5*twice;
 out->centroid_relative[0]=cx/(3.*twice);out->centroid_relative[1]=cy/(3.*twice);
 out->second_relative[0]=xx/12.;out->second_relative[1]=xyint/24.;
 out->second_relative[2]=xyint/24.;out->second_relative[3]=yy/12.;
 for(int e=0;e<nc;e++){
  int en=(e+1)%nc;int64_t j=tag[e];if(j<0)return 5;
  double dx=xy[2*j]-xi,dy=xy[2*j+1]-yi,d=hypot(dx,dy);
  double ex=x[en]-x[e],ey=y[en]-y[e],len=hypot(ex,ey);
  if(!(len>0.)||!isfinite(len)||!(d>0.))return 6;
  /* CCW edge outward normal must agree with its original power plane. */
  double dot=ey*dx-ex*dy;
  if(!(dot>0.)||!isfinite(dot))return 6;
  double px=.5*(x[e]+x[en]),py=.5*(y[e]+y[en]);
  double wn=(dx*v[2*j]+dy*v[2*j+1]
    -(v[2*j]-v[2*i])*px-(v[2*j+1]-v[2*i+1])*py)/d+(z[i]-z[j])/(2.*d);
  if(!isfinite(wn))return 6;
  out->x[e]=x[e];out->y[e]=y[e];out->neighbor[e]=ids[j];out->neighbor_kind[e]=kind[j];
  out->nx[e]=dx/d;out->ny[e]=dy/d;out->length[e]=len;out->wn_midpoint[e]=wn;
 }
 for(int k=0;k<2;k++)if(!isfinite(out->centroid_relative[k]))return 6;
 for(int k=0;k<4;k++)if(!isfinite(out->second_relative[k]))return 6;
 return 0;
}
__global__ void moment_geometry_kernel(int n,const double *xy,const double *w,
 const double *v,const double *z,const int64_t *ids,const int *kind,
 const int64_t *selected,const int64_t *off,const int64_t *cand,const double *bounds,
 uint64_t epoch,MomentGeometryPolygon *out,int *status){
 int row=blockIdx.x*blockDim.x+threadIdx.x;
 if(row<n)status[row]=moment_geometry_cell(row,xy,w,v,z,ids,kind,selected,off,cand,bounds,epoch,out+row);
}
template<class T> struct MomentGeometryDevice {
 T *p=nullptr;
 bool set(size_t n,const T *host=nullptr){
  if(!n)return true;
  if(cudaMalloc(&p,n*sizeof(T))!=cudaSuccess)return false;
  return !host||cudaMemcpy(p,host,n*sizeof(T),cudaMemcpyHostToDevice)==cudaSuccess;
 }
 ~MomentGeometryDevice(){if(p)cudaFree(p);}
};
extern "C" size_t gpu_moment_geometry_polygon_size(void){return sizeof(MomentGeometryPolygon);}
extern "C" int gpu_moment_geometry_export(int backend,int threads,int n,int nr,
 const int64_t *ids,const double *xy,const double *w,const double *v,const double *z,
 const int *kind,const int64_t *selected,const int64_t *off,const int64_t *cand,
 const double *bounds,uint64_t epoch,MomentGeometryPolygon *output,int *status){
 if((backend!=0&&backend!=1)||threads<1||n<1||nr<1||!ids||!xy||!w||!v||!z||!kind||!selected||!off||!cand||!bounds||!output||!status)return 1;
 if(off[0]!=0||off[n]<0||!isfinite(bounds[0])||!isfinite(bounds[1])||!isfinite(bounds[2])||!isfinite(bounds[3])||bounds[0]>=bounds[1]||bounds[2]>=bounds[3])return 1;
 try{
  std::vector<int64_t> sorted(ids,ids+nr);std::sort(sorted.begin(),sorted.end());
  if(sorted[0]<0||std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())return 1;
  for(int j=0;j<nr;j++)if(!isfinite(xy[2*j])||!isfinite(xy[2*j+1])||!isfinite(w[j])||!isfinite(v[2*j])||!isfinite(v[2*j+1])||!isfinite(z[j])||kind[j]<0||kind[j]>2)return 1;
  for(int i=0;i<n;i++){
   if(selected[i]<0||selected[i]>=nr||off[i]>off[i+1])return 1;
   int64_t previous=-1;
   for(int64_t k=off[i];k<off[i+1];k++){
    int64_t j=cand[k];if(j<0||j>=nr||j==selected[i]||ids[j]<=previous)return 1;previous=ids[j];
   }
  }
  std::vector<MomentGeometryPolygon> tmp(n);std::vector<int> st(n,0);
  if(backend==0){
   #pragma omp parallel for num_threads(threads)
   for(int i=0;i<n;i++)st[i]=moment_geometry_cell(i,xy,w,v,z,ids,kind,selected,off,cand,bounds,epoch,tmp.data()+i);
  }else{
   MomentGeometryDevice<double> dxy,dw,dv,dz,db;MomentGeometryDevice<int64_t> did,ds,doff,dc;
   MomentGeometryDevice<int> dk,dst;MomentGeometryDevice<MomentGeometryPolygon> dout;
   if(!dxy.set(2*size_t(nr),xy)||!dw.set(nr,w)||!dv.set(2*size_t(nr),v)||!dz.set(nr,z)||!db.set(4,bounds)||!did.set(nr,ids)||!ds.set(n,selected)||!doff.set(n+1,off)||!dc.set(off[n],cand)||!dk.set(nr,kind)||!dst.set(n)||!dout.set(n))return 8;
   moment_geometry_kernel<<<(n+63)/64,64>>>(n,dxy.p,dw.p,dv.p,dz.p,did.p,dk.p,ds.p,doff.p,dc.p,db.p,epoch,dout.p,dst.p);
   if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 7;
   if(cudaMemcpy(st.data(),dst.p,n*sizeof(int),cudaMemcpyDeviceToHost)!=cudaSuccess)return 7;
   std::copy(st.begin(),st.end(),status);
   for(int s:st)if(s)return s;
   if(cudaMemcpy(tmp.data(),dout.p,n*sizeof(MomentGeometryPolygon),cudaMemcpyDeviceToHost)!=cudaSuccess)return 7;
  }
  std::copy(st.begin(),st.end(),status);for(int s:st)if(s)return s;
  std::copy(tmp.begin(),tmp.end(),output);return 0;
 }catch(...){return 8;}
}
#endif

#include "gfs_polygon_moment_20261027.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <vector>

__global__ void moment_kernel(const double *x1,const double *y1,const double *x2,const double *y2,
                              const int *offset,const double *mass,int n,gfs_polygon_moment *out) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n)out[i]=gfs_polygon_moments(x1,y1,x2,y2,offset[i],offset[i+1],mass[i]);
}

static bool ok(cudaError_t e,const char *where) {
    if(e==cudaSuccess)return true;
    std::fprintf(stderr,"%s: %s\n",where,cudaGetErrorString(e));return false;
}

template<class T> static T *upload(const std::vector<T> &v) {
    T *d=nullptr;
    if(!ok(cudaMalloc((void**)&d,v.size()*sizeof(T)),"cudaMalloc"))return nullptr;
    if(!ok(cudaMemcpy(d,v.data(),v.size()*sizeof(T),cudaMemcpyHostToDevice),"cudaMemcpy")) {
        cudaFree(d);return nullptr;
    }
    return d;
}

int main() {
    const int n=55351,faces=4*n,iterations=100;
    std::vector<double> x1(faces),y1(faces),x2(faces),y2(faces),mass(n),width(n),height(n);
    std::vector<int> offset(n+1);
    for(int i=0;i<n;++i) {
        double w=.01+(i%7)*.0003,h=.011+(i%5)*.0002;
        double cx=((i%11)-5)*.0002,cy=((i%13)-6)*.0002;
        double x[4]={cx-w/2,cx+w/2,cx+w/2,cx-w/2};
        double y[4]={cy-h/2,cy-h/2,cy+h/2,cy+h/2};
        width[i]=w;height[i]=h;mass[i]=.001+(i%3)*.0001;offset[i]=4*i;
        for(int k=0;k<4;++k) {
            int f=4*i+k,j=(k+1)%4;x1[f]=x[k];y1[f]=y[k];x2[f]=x[j];y2[f]=y[j];
        }
    }
    offset[n]=faces;
    double *dx1=upload(x1),*dy1=upload(y1),*dx2=upload(x2),*dy2=upload(y2),*dm=upload(mass);
    int *doffset=upload(offset);
    gfs_polygon_moment *dout=nullptr;
    if(!dx1||!dy1||!dx2||!dy2||!dm||!doffset ||
       !ok(cudaMalloc((void**)&dout,n*sizeof(*dout)),"output allocation"))return 1;
    dim3 block(256),grid((n+255)/256);
    moment_kernel<<<grid,block>>>(dx1,dy1,dx2,dy2,doffset,dm,n,dout);
    if(!ok(cudaGetLastError(),"warmup launch")||!ok(cudaDeviceSynchronize(),"warmup"))return 1;
    cudaEvent_t begin,end;
    if(!ok(cudaEventCreate(&begin),"event create")||!ok(cudaEventCreate(&end),"event create"))return 1;
    cudaEventRecord(begin);
    for(int k=0;k<iterations;++k)moment_kernel<<<grid,block>>>(dx1,dy1,dx2,dy2,doffset,dm,n,dout);
    cudaEventRecord(end);
    if(!ok(cudaGetLastError(),"timed launch")||!ok(cudaEventSynchronize(end),"timed sync"))return 1;
    float elapsed=0.;cudaEventElapsedTime(&elapsed,begin,end);
    std::vector<gfs_polygon_moment> out(n);
    if(!ok(cudaMemcpy(out.data(),dout,n*sizeof(*dout),cudaMemcpyDeviceToHost),"result copy"))return 1;
    double max_area_rel=0.,max_inertia_rel=0.;int bad=0;
    for(int i=0;i<n;++i) {
        double area=width[i]*height[i];
        double inertia=mass[i]*(width[i]*width[i]+height[i]*height[i])/12.;
        if(out[i].status!=0)++bad;
        max_area_rel=std::fmax(max_area_rel,std::fabs(out[i].area/area-1.));
        max_inertia_rel=std::fmax(max_inertia_rel,std::fabs(out[i].inertia/inertia-1.));
    }
    std::printf("cells=%d faces=%d iterations=%d mean_kernel_ms=%.9g max_area_rel=%.9g max_inertia_rel=%.9g invalid=%d\n",
                n,faces,iterations,elapsed/iterations,max_area_rel,max_inertia_rel,bad);
    cudaEventDestroy(end);cudaEventDestroy(begin);
    cudaFree(dout);cudaFree(doffset);cudaFree(dm);cudaFree(dy2);cudaFree(dx2);cudaFree(dy1);cudaFree(dx1);
    if(bad || max_area_rel>1e-11 || max_inertia_rel>1e-11)return 1;
    std::printf("PASS: CUDA polygon centroid and inertia synthetic batch\n");
    return 0;
}

/* Synthetic connectivity microbenchmark; not a full simulation benchmark. */
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Exam/gfs_weight_exchange_gpu.cuh"
static void ck(cudaError_t e){if(e!=cudaSuccess){fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));exit(1);}}
int main(){
    const int nx=235,ny=235,nc=nx*ny;std::vector<GfsWeightExchangeFace> faces;
    const double dx=.01;
    auto add=[&](int i,int j){
        double xi=(i%nx)*dx,yi=(i/nx)*dx,xj=(j%nx)*dx,yj=(j/nx)*dx;
        double sx=xj-xi,sy=yj-yi,d=std::hypot(sx,sy),mx=.5*(xi+xj),my=.5*(yi+yj);
        GfsWeightExchangeFace f={};f.i=i;f.j=j;f.distance=d;
        f.wdot_i=1e-5*std::sin(.013*i);f.wdot_j=1e-5*std::sin(.013*j);
        f.rho=1.;f.ux=.2;f.uy=-.3;f.specific_energy=2.;
        f.ax=mx+.5*dx*sy/d;f.ay=my-.5*dx*sx/d;f.bx=mx-.5*dx*sy/d;f.by=my+.5*dx*sx/d;
        f.ox_i=xi;f.oy_i=yi;f.ox_j=xj;f.oy_j=yj;faces.push_back(f);
    };
    for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
        int i=y*nx+x;if(x+1<nx)add(i,i+1);if(y+1<ny)add(i,i+nx);
        if(x+1<nx && y+1<ny)add(i,i+nx+1);
    }
    int nf=(int)faces.size();std::vector<GfsMomentState> expected(nc),actual(nc);
    for(const auto &f:faces){double q[4];
        if(!gfs_weight_exchange_flux(f.wdot_i,f.wdot_j,f.distance,f.rho,f.ux,f.uy,f.specific_energy,q) ||
           !gfs_weight_exchange_segment(&expected[f.i],&expected[f.j],f.ax,f.ay,f.bx,f.by,f.ox_i,f.oy_i,f.ox_j,f.oy_j,q))return 2;
    }
    GfsWeightExchangeFace *df;GfsMomentState *dr;int *invalid;
    ck(cudaMalloc(&df,faces.size()*sizeof(faces[0])));ck(cudaMalloc(&dr,nc*sizeof(*dr)));ck(cudaMalloc(&invalid,sizeof(int)));
    ck(cudaMemcpy(df,faces.data(),faces.size()*sizeof(faces[0]),cudaMemcpyHostToDevice));
    cudaEvent_t a,b;ck(cudaEventCreate(&a));ck(cudaEventCreate(&b));std::vector<float> ms;
    for(int rep=0;rep<13;rep++){
        ck(cudaMemset(dr,0,nc*sizeof(*dr)));ck(cudaMemset(invalid,0,sizeof(int)));
        ck(cudaEventRecord(a));gfs_weight_exchange_faces_kernel<<<(nf+127)/128,128>>>(df,nf,nc,dr,invalid);
        ck(cudaGetLastError());ck(cudaEventRecord(b));ck(cudaEventSynchronize(b));
        float elapsed;ck(cudaEventElapsedTime(&elapsed,a,b));if(rep>=3)ms.push_back(elapsed);
    }
    int bad;ck(cudaMemcpy(&bad,invalid,sizeof(bad),cudaMemcpyDeviceToHost));if(bad)return 3;
    ck(cudaMemcpy(actual.data(),dr,nc*sizeof(*dr),cudaMemcpyDeviceToHost));
    long double err2=0,ref2=0,mass=0,px=0,py=0,ene=0,L=0;double maxabs=0;
    for(int k=0;k<nc;k++){
        for(int q=0;q<GFM_NSTATE;q++){
            double z=actual[k].a[q]-expected[k].a[q];maxabs=std::max(maxabs,std::abs(z));
            err2+=(long double)z*z;ref2+=(long double)expected[k].a[q]*expected[k].a[q];
        }
        mass+=actual[k].a[GFM_M];px+=actual[k].a[GFM_PX];py+=actual[k].a[GFM_PY];ene+=actual[k].a[GFM_E];
        L+=gfs_moment_Lz(&actual[k],(k%nx)*dx,(k/nx)*dx);
    }
    double rel=std::sqrt((double)(err2/ref2));if(!(rel<1e-11) || maxabs>1e-12)return 4;
    std::sort(ms.begin(),ms.end());
    printf("cells=%d faces=%d median_kernel_ms=%.6f min_ms=%.6f max_ms=%.6f rel_rms=%.3g max_abs=%.3g mass=%.3Lg px=%.3Lg py=%.3Lg energy=%.3Lg L=%.3Lg face_bytes=%zu rhs_bytes=%zu\n",
           nc,nf,.5*(ms[4]+ms[5]),ms.front(),ms.back(),rel,maxabs,mass,px,py,ene,L,
           faces.size()*sizeof(faces[0]),nc*sizeof(*dr));
    ck(cudaEventDestroy(a));ck(cudaEventDestroy(b));ck(cudaFree(df));ck(cudaFree(dr));ck(cudaFree(invalid));
}

/* Real CUDA kernel versus double exhaustive search; no simulation evolution. */
#include <cstdio>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cfloat>
extern "C" {
#include "exam_gpu.h"
}
int main(){
    GPUContext ctx={};
    const int n=64, groups=8, total=n+1+groups;
    gpu_init(&ctx,256,512,0);
    gpu_alloc_tree_buffers(&ctx,256);
    for(int test=0;test<4;test++){
        double shift=test==2?1e8:0, scale=test==2?1e-5:1;
        for(int i=0;i<n;i++){
            ctx.h_parts.x[i]=shift+scale*(i/8+.01*std::sin(7.*i));
            ctx.h_parts.y[i]=shift+scale*(.1*(i%8)+.001*std::cos(11.*i));
            ctx.h_parts.indx[i]=i;
        }
        if(test==1){ ctx.h_parts.x[1]=ctx.h_parts.x[0];ctx.h_parts.y[1]=ctx.h_parts.y[0]; }
        if(test==3){ // explicit periodic image near the opposite boundary
            ctx.h_parts.x[0]=0.001;ctx.h_parts.y[0]=0;
            ctx.h_parts.x[63]=-0.001;ctx.h_parts.y[63]=0;
        }
        for(int i=0;i<n;i++){
            int group=i/8;
            ctx.h_tree_type[i]=1;ctx.h_tree_ptl_idx[i]=i;
            ctx.h_tree_cx[i]=ctx.h_parts.x[i];ctx.h_tree_cy[i]=ctx.h_parts.y[i];
            ctx.h_tree_size[i]=0;ctx.h_tree_daughter[i]=-1;
            ctx.h_tree_sibling[i]=(i%8!=7)?i+1:(group<7?n+2+group:-1);
        }
        for(int g=0;g<groups;g++){
            int k=n+1+g;
            double cx=0,cy=0;
            for(int j=8*g;j<8*g+8;j++){cx+=ctx.h_parts.x[j]/8;cy+=ctx.h_parts.y[j]/8;}
            double radius=0;
            for(int j=8*g;j<8*g+8;j++) radius=std::max(radius,std::hypot(ctx.h_parts.x[j]-cx,ctx.h_parts.y[j]-cy));
            ctx.h_tree_type[k]=0;ctx.h_tree_ptl_idx[k]=-1;
            ctx.h_tree_cx[k]=cx;ctx.h_tree_cy[k]=cy;ctx.h_tree_size[k]=std::nextafter(radius,INFINITY);
            ctx.h_tree_daughter[k]=8*g;ctx.h_tree_sibling[k]=g<7?k+1:-1;
        }
        ctx.h_tree_type[n]=0;ctx.h_tree_ptl_idx[n]=-1;
        ctx.h_tree_cx[n]=shift;ctx.h_tree_cy[n]=shift;ctx.h_tree_size[n]=100;
        ctx.h_tree_daughter[n]=n+1;ctx.h_tree_sibling[n]=-1;
        ctx.tree_root_idx=n;
        gpu_upload_particles(&ctx,n,0);gpu_upload_tree(&ctx,total);
        gpu_launch_nearest_neighbor_tree(&ctx,n,n,0);
        double maxrel=0;
        for(int i=0;i<n;i++){
            double ref=1e23;
            for(int j=0;j<n;j++)if(j!=i){
                double dx=ctx.h_parts.x[i]-ctx.h_parts.x[j],dy=ctx.h_parts.y[i]-ctx.h_parts.y[j];
                ref=std::min(ref,dx*dx+dy*dy);
            }
            double got=ctx.h_tree_nn2[i];assert(std::isfinite(got));
            if(ref==0)assert(got==0);
            else{double rel=std::abs(got-ref)/ref;maxrel=std::max(maxrel,rel);assert(rel<32*DBL_EPSILON);}
        }
        std::printf("case=%d queries=%d max_relative_error=%.17g PASS\n",test,n,maxrel);
    }
    ctx.h_tree_daughter[n]=n; // broken topology must not return a finite result
    gpu_upload_tree(&ctx,total);gpu_launch_nearest_neighbor_tree(&ctx,n,n,0);
    for(int i=0;i<n;i++)assert(std::isnan(ctx.h_tree_nn2[i]));
    std::puts("cycle_detection PASS");
    gpu_free(&ctx);
}

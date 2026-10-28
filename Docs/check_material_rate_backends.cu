#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include <cuda_runtime.h>
#include "moment_material_drift.h"

static void check(cudaError_t e) {
    if(e!=cudaSuccess){fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));exit(2);}
}
static void save(const char *path,const std::vector<double>& a) {
    FILE *f=fopen(path,"wx");if(!f){perror(path);exit(3);}
    if(fwrite(a.data(),sizeof(double),a.size(),f)!=a.size()){perror("write");exit(4);}
    if(fclose(f)!=0)exit(5);
}
__global__ void kernel(const double *a,double *b,size_t n) {
    const size_t k=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if(k<n)md_material_rate(a+17*k,b+17*k);
}
int main(int argc,char **argv) {
    if(argc!=5)return 1;
    FILE *f=fopen(argv[1],"rb");if(!f)return 2;
    if(fseek(f,0,SEEK_END)!=0)return 2;
    long size=ftell(f);if(size<=0 || size%(17*sizeof(double)))return 2;
    rewind(f);size_t n=size/(17*sizeof(double));
    std::vector<double>a(17*n),b(17*n);
    if(fread(a.data(),sizeof(double),a.size(),f)!=a.size())return 2;
    fclose(f);
    for(int mode=0;mode<2;mode++) {
        omp_set_dynamic(0);omp_set_num_threads(mode?4:1);
        int observed=0;
        #pragma omp parallel
        {
            #pragma omp single
            observed=omp_get_num_threads();
            #pragma omp for
            for(size_t k=0;k<n;k++)md_material_rate(a.data()+17*k,b.data()+17*k);
        }
        if(observed!=(mode?4:1))return 6;
        printf("CPU threads=%d faces=%zu\n",observed,n);save(argv[2+mode],b);
    }
    double *da=nullptr,*db=nullptr;
    check(cudaMalloc(&da,a.size()*sizeof(double)));check(cudaMalloc(&db,b.size()*sizeof(double)));
    check(cudaMemcpy(da,a.data(),a.size()*sizeof(double),cudaMemcpyHostToDevice));
    kernel<<<(n+255)/256,256>>>(da,db,n);check(cudaGetLastError());check(cudaDeviceSynchronize());
    check(cudaMemcpy(b.data(),db,b.size()*sizeof(double),cudaMemcpyDeviceToHost));save(argv[4],b);
    check(cudaFree(da));check(cudaFree(db));printf("CUDA faces=%zu\n",n);return 0;
}

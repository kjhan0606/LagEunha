/* Standalone backend contract probe, not wired into live hydro. */
#include <cuda_runtime.h>
#include <omp.h>
#include <cstring>
__global__ void assemble_gpu(int nr,const int *ptr,const int *col,const double *b,
        const double *v,const double *z,double *dz,double *force){
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=nr)return;
    double a=0,c=0;
    for(int k=ptr[i];k<ptr[i+1];k++){
        int j=col[k];double q=b[k];
        a+=q*v[2*j];c+=q*v[2*j+1];
        atomicAdd(force+2*j,-q*z[2*i]);atomicAdd(force+2*j+1,-q*z[2*i+1]);
    }
    dz[2*i]=a;dz[2*i+1]=c;
}
extern "C" int memory_assemble(int n,int nr,const int *ptr,const int *col,const double *b,
        const double *v,const double *z,double *dz,double *force,int mode){
    if(n<1 || nr<1 || mode<0 || mode>2)return -1;
    int nnz=ptr[nr];
    if(mode<2){
        std::memset(force,0,2*n*sizeof(double));
        #pragma omp parallel for if(mode==1) schedule(static) num_threads(4)
        for(int i=0;i<nr;i++){
            double a=0,c=0;
            for(int k=ptr[i];k<ptr[i+1];k++){
                int j=col[k];double q=b[k];a+=q*v[2*j];c+=q*v[2*j+1];
                #pragma omp atomic update
                force[2*j]-=q*z[2*i];
                #pragma omp atomic update
                force[2*j+1]-=q*z[2*i+1];
            }
            dz[2*i]=a;dz[2*i+1]=c;
        }
        return 0;
    }
    int *ip=nullptr;double *dp=nullptr;int result=0;cudaError_t err;
    size_t count=(size_t)nnz+4*(size_t)n+4*(size_t)nr;
    #define CHECK(x) do{err=(x);if(err!=cudaSuccess){result=-1000-(int)err;goto done;}}while(0)
    CHECK(cudaMalloc((void**)&ip,((size_t)nr+1+nnz)*sizeof(int)));
    CHECK(cudaMalloc((void**)&dp,count*sizeof(double)));
    CHECK(cudaMemcpy(ip,ptr,(nr+1)*sizeof(int),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(ip+nr+1,col,nnz*sizeof(int),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dp,b,nnz*sizeof(double),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dp+nnz,v,2*n*sizeof(double),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dp+nnz+2*n,z,2*nr*sizeof(double),cudaMemcpyHostToDevice));
    CHECK(cudaMemset(dp+nnz+2*n+4*nr,0,2*n*sizeof(double)));
    assemble_gpu<<<(nr+127)/128,128>>>(nr,ip,ip+nr+1,dp,dp+nnz,dp+nnz+2*n,
                                      dp+nnz+2*n+2*nr,dp+nnz+2*n+4*nr);
    CHECK(cudaGetLastError());CHECK(cudaDeviceSynchronize());
    CHECK(cudaMemcpy(dz,dp+nnz+2*n+2*nr,2*nr*sizeof(double),cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(force,dp+nnz+2*n+4*nr,2*n*sizeof(double),cudaMemcpyDeviceToHost));
done:
    if(ip)cudaFree(ip);if(dp)cudaFree(dp);return result;
    #undef CHECK
}

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <algorithm>
#ifdef MGS_MPI
#include <mpi.h>
#endif
#include <omp.h>
#include <cuda_runtime.h>
#include "moment_gravity_quadrature.h"
static void check(cudaError_t e){if(e!=cudaSuccess){fprintf(stderr,"CUDA %s\n",cudaGetErrorString(e));exit(2);}}
template<class T> std::vector<T> read(const char *p){
    FILE *f=fopen(p,"rb");if(!f)exit(3);fseek(f,0,SEEK_END);long n=ftell(f);rewind(f);
    if(n<=0 || n%sizeof(T))exit(4);std::vector<T>a(n/sizeof(T));
    if(fread(a.data(),sizeof(T),a.size(),f)!=a.size())exit(5);fclose(f);return a;
}
static void save(const char *p,const std::vector<double>&a){FILE*f=fopen(p,"wx");if(!f)exit(6);if(fwrite(a.data(),sizeof(double),a.size(),f)!=a.size())exit(7);if(fclose(f))exit(8);}
static __host__ __device__ void cell(const double*a,const uint64_t*offset,double*b,size_t i){
    for(int k=0;k<11;k++)b[11*i+k]=12345.;
    b[11*i]=mgs_integrate(a+7*offset[i],offset[i+1]-offset[i],b+11*i+1);
}
__global__ void kernel(const double*a,const uint64_t*o,double*b,size_t n,int rank,int world){size_t i=rank+((size_t)blockIdx.x*blockDim.x+threadIdx.x)*world;if(i<n)cell(a,o,b,i);}
static void collect(std::vector<double>&b){
#ifdef MGS_MPI
    if(MPI_Allreduce(MPI_IN_PLACE,b.data(),(int)b.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD)!=MPI_SUCCESS)MPI_Abort(MPI_COMM_WORLD,12);
#endif
}
int main(int argc,char**argv){
    int rank=0,world=1;
#ifdef MGS_MPI
    int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    if(provided<MPI_THREAD_FUNNELED)MPI_Abort(MPI_COMM_WORLD,13);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&world);
#endif
    if(argc!=6)return 1;auto a=read<double>(argv[1]);auto o=read<uint64_t>(argv[2]);
    if(o.size()<2 || o[0]!=0 || o.back()*7!=a.size())return 9;
    size_t n=o.size()-1;for(size_t i=0;i<n;i++)if(o[i+1]<o[i])return 10;
    std::vector<double>b(n*11);
    for(int mode=0;mode<2;mode++){
        std::fill(b.begin(),b.end(),0.);
        omp_set_dynamic(0);omp_set_num_threads(mode?4:1);int observed=0;
        #pragma omp parallel
        {
            #pragma omp single
            observed=omp_get_num_threads();
            #pragma omp for
            for(size_t i=rank;i<n;i+=world)cell(a.data(),o.data(),b.data(),i);
        }
        if(observed!=(mode?4:1))return 11;printf("CPU rank=%d/%d threads=%d total_cells=%zu\n",rank,world,observed,n);collect(b);if(rank==0)save(argv[3+mode],b);
    }
    double *da,*db;uint64_t*doo;check(cudaMalloc(&da,a.size()*sizeof(double)));check(cudaMalloc(&db,b.size()*sizeof(double)));check(cudaMalloc(&doo,o.size()*sizeof(uint64_t)));
    check(cudaMemcpy(da,a.data(),a.size()*sizeof(double),cudaMemcpyHostToDevice));check(cudaMemcpy(doo,o.data(),o.size()*sizeof(uint64_t),cudaMemcpyHostToDevice));
    check(cudaMemset(db,0,b.size()*sizeof(double)));
    kernel<<<(n+127)/128,128>>>(da,doo,db,n,rank,world);check(cudaGetLastError());check(cudaDeviceSynchronize());
    check(cudaMemcpy(b.data(),db,b.size()*sizeof(double),cudaMemcpyDeviceToHost));collect(b);if(rank==0)save(argv[5],b);
    check(cudaFree(da));check(cudaFree(db));check(cudaFree(doo));printf("CUDA rank=%d/%d total_cells=%zu\n",rank,world,n);
#ifdef MGS_MPI
    MPI_Finalize();
#endif
    return 0;
}

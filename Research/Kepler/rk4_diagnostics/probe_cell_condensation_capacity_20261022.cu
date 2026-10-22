#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
int main(int argc,char**argv){
    if(argc!=3)return 2;
    uint64_t bytes=strtoull(argv[1],nullptr,10);cudaDeviceProp p{};size_t available=0,total=0;
    cudaError_t e=cudaGetDeviceProperties(&p,0);if(e==cudaSuccess)e=cudaMemGetInfo(&available,&total);
    bool capacity=e==cudaSuccess&&strstr(p.name,"H200")&&bytes<=uint64_t(total/5*4)&&bytes<=uint64_t(available/5*4);
    void*ptr=nullptr;cudaError_t allocation=cudaErrorInvalidValue;
    if(capacity)allocation=cudaMalloc(&ptr,size_t(bytes));
    bool pass=capacity&&allocation==cudaSuccess;if(ptr)cudaFree(ptr);
    FILE*f=fopen(argv[2],"wx");if(!f)return 2;
    fprintf(f,"{\"scope\":\"single planned aggregate allocation proof; actual MPI runtime accounting remains required\",\"device\":\"%s\",\"runtime_status\":%d,\"allocation_status\":%d,\"total_bytes\":%llu,\"free_before_bytes\":%llu,\"planned_all_ranks_bytes\":%llu,\"twenty_percent_headroom\":%s,\"planned_allocation_pass\":%s,\"physical_evolution\":false}\n",p.name,int(e),int(allocation),(unsigned long long)total,(unsigned long long)available,(unsigned long long)bytes,capacity?"true":"false",pass?"true":"false");
    if(fclose(f))return 2;return pass?0:2;
}

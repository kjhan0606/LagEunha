/* Execute the same new dimensional FV primitives on the actual CUDA device.
 * Backend parity supplements, rather than replaces, the independent C audit.
 * No geometry or PDE is launched. */
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "../gfs_fv_ale.h"

#define CHECK(call) do {cudaError_t e=(call);if(e!=cudaSuccess){ \
    std::fprintf(stderr,"CUDA audit: %s\n",cudaGetErrorString(e));return 2;}}while(0)

__host__ __device__ static void evaluate(int k,double out[16])
{
    int dim=1+k%3,test=k/3;double n[3]={0};
    GfsFvPrimitive l={1,1,{0,0,0}},r={.125,.1,{0,0,0}};
    for(int a=0;a<dim;a++){
        n[a]=1/sqrt((double)dim);
        l.u[a]=sin(.3*(test+1)*(a+1));r.u[a]=cos(.4*(test+1)*(a+1));
        if(test%4==0){l.u[a]=-8*n[a];r.u[a]=8*n[a];}
        if(test%4==1){l.u[a]=8*n[a];r.u[a]=-8*n[a];}
        if(test%4==2){l.u[a]=r.u[a]=.5*n[a];r.p=1;}
    }
    GfsFvFan fan=gfs_fv_fan(dim,l,r,n,1.4);
    GfsFvSample sample=gfs_fv_sample(&fan,.2*(test-3));
    for(int c=0;c<5;c++)out[c]=sample.flux[c];
    double rate[5]={0};rate[4]=-fan.state[0][4]-1;
    for(int a=0;a<dim;a++)rate[1+a]=.2*(a+1);
    out[5]=sample.density;out[6]=gfs_fv_admissible_interval(dim,fan.state[0],rate);
    out[7]=fan.valid;out[8]=sample.valid;out[9]=fan.hll;
    out[10]=gfs_fv_internal(dim,fan.state[0]);
    double qs=0,rs=0,x[5]={0};
    for(int c=0;c<5;c++){qs=fmax(qs,fabs(fan.state[0][c]));rs=fmax(rs,fabs(rate[c]));}
    for(int c=0;c<5;c++)x[c]=fan.state[0][c]/qs;
    out[11]=qs;out[12]=rs;out[13]=gfs_fv_internal(dim,x);out[14]=x[0];out[15]=2*x[0]*out[13];
}
__global__ static void evaluate_kernel(double *out)
{
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k<96)evaluate(k,out+16*k);
}
int main(void)
{
    double *device,actual[1536];CHECK(cudaMalloc(&device,sizeof(actual)));
    evaluate_kernel<<<1,128>>>(device);CHECK(cudaGetLastError());
    CHECK(cudaMemcpy(actual,device,sizeof(actual),cudaMemcpyDeviceToHost));CHECK(cudaFree(device));
    cudaDeviceProp properties;CHECK(cudaGetDeviceProperties(&properties,0));
    double error=0;int bad=0,fallback=0,shown=0;
    FILE *trace=std::fopen("fv_device_values.bin","wbx");if(!trace)return 2;
    if(std::fwrite(actual,sizeof(actual),1,trace)!=1)return 2;
    for(int k=0;k<96;k++){
        double expected[16];evaluate(k,expected);
        if(std::fwrite(expected,sizeof(expected),1,trace)!=1)return 2;
        if(actual[16*k+7]!=1||actual[16*k+8]!=1||!(actual[16*k+6]>0))bad=1;
        fallback+=(int)actual[16*k+9];
        for(int c=0;c<16;c++){
            double a=actual[16*k+c],b=expected[c];
            if(!std::isfinite(a)||!std::isfinite(b))bad=1;
            double mismatch=fabs(a-b)/(1+fabs(a)+fabs(b));
            error=fmax(error,mismatch);
            if(mismatch>1e-10&&shown++<12)
                std::fprintf(stderr,"FV_DEVICE_DIFFERENCE case=%d dim=%d component=%d GPU=%.17g CPU=%.17g normalized=%.17g\n",k,1+k%3,c,a,b,mismatch);
        }
    }
    if(std::fclose(trace))return 2;
    int pass=!bad&&error<1e-12;
    std::printf("{\"scope\":\"common dimensional fan and cone executed on CUDA; no PDE\","
                "\"device\":\"%s\",\"cases\":96,\"HLL_cases\":%d,\"CPU_GPU_normalized_max_error\":%.17g,"
                "\"PDE_steps\":0,\"pass\":%s}\n",properties.name,fallback,error,pass?"true":"false");
    return pass?0:1;
}

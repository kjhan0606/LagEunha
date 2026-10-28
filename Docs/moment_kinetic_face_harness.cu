#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
#include <float.h>
#include "moment_kinetic_face_flux.h"
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define KTEST static __host__ __device__
#else
#define KTEST static
#endif

KTEST int synthetic_checks(){
    double a[33]={0.},f[18];int failures=0;
    a[0]=.3;a[1]=-.4;a[2]=.6;a[3]=.8;a[4]=.25;a[31]=2.;a[32]=.7;
    failures+=(lmk_face_node_flux(a,f)!=0);
    // No relative transport: pressure force, torque and work still remain.
    failures+=(fabs(f[7]-.3)>1e-15 || fabs(f[8]-.4)>1e-15 || fabs(f[17]-.35)>1e-15);
    failures+=(fabs((f[10]-f[11])-.24)>1e-15);
    for(int k=0;k<7;k++)failures+=(f[k]!=0.);
    for(int k=13;k<17;k++)failures+=(f[k]!=0.);
    for(int k=0;k<18;k++)failures+=((f[k]+(-f[k]))!=0.); // same owner-computed flux
    a[31]=0.;a[29]=a[30]=.3;a[21]=.2;a[25]=-.1;a[4]=.7;
    for(int k=0;k<4;k++){a[5+2*k]=a[13+2*k]=1.25;a[6+2*k]=a[14+2*k]=-.5;}
    failures+=(lmk_face_node_flux(a,f)!=0);
    const double energy=.5*(1.25*1.25+.5*.5)+.3;
    failures+=(fabs(f[7]-1.25*f[0])>1e-15 || fabs(f[8]+.5*f[0])>1e-15 || fabs(f[17]-energy*f[0])>1e-15);
    failures+=(f[4]!=f[5] || f[14]!=f[15]);
    for(int kind=0;kind<5;kind++){
        double bad[33],out[18];for(int k=0;k<33;k++)bad[k]=a[k];
        for(int k=0;k<18;k++)out[k]=-713.25;
        int expected=1;
        if(kind==0)bad[4]=-1.;
        if(kind==1)bad[21]=-1.;
        if(kind==2)bad[25]=1.;
        if(kind==3)bad[0]=NAN;
        if(kind==4){bad[5]=DBL_MAX;bad[21]=DBL_MAX;expected=2;}
        failures+=(lmk_face_node_flux(bad,out)!=expected);
        for(int k=0;k<18;k++)failures+=(out[k]!=-713.25);
    }
    return failures;
}

KTEST void one(const double *input,double *output,int index){
    output[18]=lmk_face_node_flux(input,output);
    if(index==0)output[18]+=synthetic_checks();
}
#ifdef __CUDACC__
__global__ void kernel(const double *in,double *out,int n){
    const int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)one(in+33*i,out+19*i,i);
}
#endif
int main(int argc,char **argv){
    if(argc!=4)return 2;
    std::ifstream f(argv[1],std::ios::binary|std::ios::ate);if(!f)return 3;
    const auto bytes=f.tellg();if(bytes<=0 || bytes%(33*sizeof(double)))return 4;
    const int n=int(bytes/(33*sizeof(double)));std::vector<double> in(33*n),out(19*n);
    f.seekg(0);f.read(reinterpret_cast<char*>(in.data()),bytes);if(!f)return 5;
    if(std::strcmp(argv[3],"cpu")==0){
        #pragma omp parallel for
        for(int i=0;i<n;i++)one(in.data()+33*i,out.data()+19*i,i);
    }else{
#ifdef __CUDACC__
        double *di=nullptr,*dout=nullptr;
        if(cudaMalloc(&di,in.size()*sizeof(double))!=cudaSuccess || cudaMalloc(&dout,out.size()*sizeof(double))!=cudaSuccess)return 6;
        if(cudaMemcpy(di,in.data(),in.size()*sizeof(double),cudaMemcpyHostToDevice)!=cudaSuccess)return 7;
        kernel<<<(n+127)/128,128>>>(di,dout,n);
        if(cudaGetLastError()!=cudaSuccess || cudaDeviceSynchronize()!=cudaSuccess)return 8;
        if(cudaMemcpy(out.data(),dout,out.size()*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return 9;
        cudaFree(di);cudaFree(dout);
#else
        return 10;
#endif
    }
    for(int i=0;i<n;i++)if(out[19*i+18]!=0.){std::fprintf(stderr,"node %d status %g\n",i,out[19*i+18]);return 11;}
    std::ofstream result(argv[2],std::ios::binary);result.write(reinterpret_cast<const char*>(out.data()),out.size()*sizeof(double));if(!result)return 12;
    std::printf("{\"passed\":true,\"nodes\":%d,\"backend\":\"%s\"}\n",n,argv[3]);return 0;
}

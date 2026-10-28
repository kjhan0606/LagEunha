#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstring>
#include "moment_state_layout.h"
#include "weight_remap_moments.h"
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define TINLINE static __host__ __device__
#else
#define TINLINE static
#endif

TINLINE int reject_raw(const double *input,int kind) {
    double raw[18],out[17];for(int k=0;k<18;k++)raw[k]=input[k];
    for(int k=0;k<17;k++)out[k]=-713.25;
    LagMomentPsi p={LMP_DIAGNOSTIC_RECONSTRUCTION,-7.25},meta={-23,-77.25};
    int expected=LMA_ASYMMETRIC;
    if(kind==0)raw[LMR_QXY]=raw[LMR_QYX]+1.;
    if(kind==1)raw[LMR_PIXY]=raw[LMR_PIYX]+1.;
    if(kind==2){p.mode=LMP_ABSENT;expected=LMA_PSI_ABSENT;}
    if(kind==3){p.mode=99;expected=LMA_INVALID_PSI_MODE;}
    if(kind==4){raw[LMR_M]=NAN;expected=LMA_NONFINITE;}
    if(kind==5){p.value=INFINITY;expected=LMA_NONFINITE;}
    if(kind==6){raw[LMR_QXY]=0.;raw[LMR_QYX]=-0.;}
    if(kind==7){raw[LMR_PIXY]=0.;raw[LMR_PIYX]=-0.;}
    const int status=lml_raw18_to_wm17(raw,&p,out,&meta);
    if(status!=expected || meta.mode!=-23 || meta.value!=-77.25)return 1;
    for(int k=0;k<17;k++)if(out[k]!=-713.25)return 1;
    return 0;
}

TINLINE void one(const double *raw,double *out) {
    for(int k=0;k<64;k++)out[k]=0.;
    double compact[17],back[18];LagMomentPsi psi={LMP_DIAGNOSTIC_RECONSTRUCTION,-7.25},meta={0,0.},backmeta={0,0.};
    out[1]=lml_raw18_to_wm17(raw,&psi,compact,&meta);
    out[2]=lml_wm17_to_raw18(compact,meta.mode,back,&backmeta);
    out[0]+=(out[1]!=LMA_OK || out[2]!=LMA_OK);
    for(int k=0;k<18;k++){out[3+k]=back[k];out[0]+=(back[k]!=raw[k]);}
    out[21]=backmeta.mode;out[22]=backmeta.value;
    out[0]+=(backmeta.mode!=psi.mode || backmeta.value!=psi.value);
    for(int kind=0;kind<8;kind++){out[23+kind]=reject_raw(raw,kind);out[0]+=out[23+kind];}
    double alias[18];for(int k=0;k<18;k++)alias[k]=raw[k];
    out[39]=lml_raw18_to_wm17(alias,&psi,alias,&meta);
    out[40]=lml_wm17_to_raw18(alias,meta.mode,alias,&backmeta);
    out[0]+=(out[39]!=LMA_OK || out[40]!=LMA_OK);
    for(int k=0;k<18;k++)out[0]+=(alias[k]!=raw[k]);
    // Signed rates may have negative mass/energy; no state projection occurs.
    double signedraw[18];for(int k=0;k<18;k++)signedraw[k]=raw[k];
    signedraw[0]=-2.;signedraw[17]=-3.;
    out[31]=lml_raw18_to_wm17(signedraw,&psi,compact,&meta);out[0]+=(out[31]!=LMA_OK);
    // Independently generated legacy fixture with asymmetric D.
    const double atom[7]={2.,.3,-.7,1.1,-.4,.2,-3.};double wm[17];wm_atom(atom,wm);
    out[32]=lml_wm17_to_raw18(wm,LMP_TRANSPORTED_MOMENT,back,&backmeta);
    out[0]+=(out[32]!=LMA_OK || backmeta.mode!=LMP_TRANSPORTED_MOMENT);
    out[33]=back[LMR_DXY];out[34]=back[LMR_DYX];
    out[0]+=(back[LMR_DXY]!=wm[9] || back[LMR_DYX]!=wm[10] || back[LMR_DXY]==back[LMR_DYX]);
    // Reverse failures must also be transactional, including metadata.
    for(int kind=0;kind<3;kind++){
        wm_atom(atom,wm);for(int k=0;k<18;k++)back[k]=-713.25;
        backmeta.mode=-23;backmeta.value=-77.25;
        int mode=LMP_DIAGNOSTIC_RECONSTRUCTION,expected=LMA_NONFINITE;
        if(kind==0)wm[16]=NAN;
        if(kind==1){mode=LMP_ABSENT;expected=LMA_PSI_ABSENT;}
        if(kind==2){mode=99;expected=LMA_INVALID_PSI_MODE;}
        int status=lml_wm17_to_raw18(wm,mode,back,&backmeta);
        out[0]+=(status!=expected || backmeta.mode!=-23 || backmeta.value!=-77.25);
        for(int k=0;k<18;k++)out[0]+=(back[k]!=-713.25);
    }
}

#ifdef __CUDACC__
__global__ void kernel(const double *in,double *out,int n){
    const int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)one(in+18*i,out+64*i);
}
#endif
int main(int argc,char **argv){
    if(argc!=4)return 2;
    std::ifstream f(argv[1],std::ios::binary|std::ios::ate);if(!f)return 3;
    const auto bytes=f.tellg();if(bytes<=0 || bytes%(18*sizeof(double)))return 4;
    int n=int(bytes/(18*sizeof(double)));std::vector<double> in(18*n),out(64*n);
    f.seekg(0);f.read(reinterpret_cast<char*>(in.data()),bytes);if(!f)return 5;
    if(std::strcmp(argv[3],"cpu")==0){
        #pragma omp parallel for
        for(int i=0;i<n;i++)one(in.data()+18*i,out.data()+64*i);
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
    for(int i=0;i<n;i++)if(out[64*i]!=0.)return 11;
    std::ofstream result(argv[2],std::ios::binary);result.write(reinterpret_cast<const char*>(out.data()),out.size()*sizeof(double));
    if(!result)return 12;
    std::printf("{\"passed\":true,\"fixtures\":%d,\"backend\":\"%s\"}\n",n,argv[3]);return 0;
}

#define main hydro_reference_main
#include "moment_hydro_cell_harness.cu"
#undef main
#include "moment_factor_update.h"
enum {FWIDTH=63};
HD int factor_fixtures(){
    double x[4]={1.,2.,0.,1.},ex[4],ed[3];
    if(lmf_exponential(x,ex,ed)||ed[0]<=0.)return 16;
    const double ee=exp(1.);
    if(fabs(ex[0]-ee)>2e-13||fabs(ex[1]-2*ee)>2e-13||ex[2]!=0.||fabs(ex[3]-ee)>2e-13)return 17;
    x[0]=x[3]=0.;x[1]=-2.;x[2]=2.;
    if(lmf_exponential(x,ex,ed))return 18;
    if(fabs(ex[0]-cos(2.))>2e-13||fabs(ex[1]+sin(2.))>2e-13||fabs(ex[2]-sin(2.))>2e-13||fabs(ex[3]-cos(2.))>2e-13)return 19;
    double base[4]={1.,.2,.3,.8},linear[4]={.9,.1,-.2,.7},A[4]={.03,.2,-.1,-.04},K[4]={.2,.03,.03,.1},q[4]={.1,-.2,.3,.1};
    double raw[18]={0.};raw[0]=1.;raw[3]=raw[6]=1.;raw[13]=raw[16]=.2;raw[17]=2.;
    double l[4],r[18],d[40];
    for(int bypass=0;bypass<2;bypass++){
        double rankzero[4]={0.};
        if(lmf_update(base,rankzero,A,K,q,2,bypass?.1:0.,raw,bypass,l,r,d))return 1;
        for(int k=0;k<4;k++)if(l[k]!=base[k])return 2;
        for(int k=0;k<18;k++)if(r[k]!=raw[k])return 3;
    }
    if(lmf_update(base,linear,A,K,q,2,.01,raw,0,l,r,d))return 4;
    if(l[1]!=0.||l[0]<=0.||l[3]<=0.||d[31]>1e-14||d[36]>1e-13)return 5;
    double R[4],J[4];
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        R[2*i+j]=linear[2*i]*linear[2*j]+linear[2*i+1]*linear[2*j+1];
        J[2*i+j]=(linear[2*i]*K[0]+linear[2*i+1]*K[2])*linear[2*j]+(linear[2*i]*K[1]+linear[2*i+1]*K[3])*linear[2*j+1];
    }
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        double v=0.;for(int k=0;k<2;k++)v+=d[2*i+k]*R[2*k+j]+R[2*i+k]*d[2*k+j];
        if(fabs(v-J[2*i+j])>1e-14)return 6;
    }
    for(int k=0;k<18;k++)if(k<13||k>16)if(r[k]!=raw[k])return 7;
    if(r[14]!=r[15])return 8;
    if(lmf_update(base,linear,A,K,nullptr,0,.01,raw,0,l,r,d))return 9;
    double fromBase[4],fromLinear[4],contrast=0.;
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        fromBase[2*i+j]=d[12+2*i]*base[j]+d[13+2*i]*base[2+j];
        fromLinear[2*i+j]=d[12+2*i]*linear[j]+d[13+2*i]*linear[2+j];
    }
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        const double rb=fromBase[2*i]*fromBase[2*j]+fromBase[2*i+1]*fromBase[2*j+1];
        const double rl=fromLinear[2*i]*fromLinear[2*j]+fromLinear[2*i+1]*fromLinear[2*j+1];
        if(fabs(d[16+2*i+j]-rb)>1e-13)return 14;
        contrast=fmax(contrast,fabs(rb-rl));
    }
    if(contrast<.01)return 15;
    for(int kind=0;kind<5;kind++){
        double bl[4],ba[4],bk[4],br[18];for(int k=0;k<4;k++){bl[k]=linear[k];ba[k]=A[k];bk[k]=K[k];l[k]=-97.;}
        for(int k=0;k<18;k++){br[k]=raw[k];r[k]=-97.;}for(int k=0;k<40;k++)d[k]=-97.;
        double h=.01;
        if(kind==0)for(int k=0;k<4;k++)bl[k]=0.;
        if(kind==1)ba[0]=NAN;if(kind==2)h=-1.;if(kind==3)br[17]=0.;if(kind==4)bk[1]=1.;
        if(lmf_update(base,bl,ba,bk,q,2,h,br,0,l,r,d)==0)return 10;
        for(int k=0;k<4;k++)if(l[k]!=-97.)return 11;
        for(int k=0;k<18;k++)if(r[k]!=-97.)return 12;
        for(int k=0;k<40;k++)if(d[k]!=-97.)return 13;
    }
    return 0;
}
HD void factor_one(int index,int n,const double *input,const double *q,const int64_t *offset,double *out){
    const int i=index%n;const double *a=input+36*index;
    out[62]=lmf_update(a,a+4,a+8,a+12,q+2*offset[i],int(offset[i+1]-offset[i]),a[34],a+16,int(a[35]),out,out+4,out+22);
    if(index==0)out[62]+=factor_fixtures();
}
#ifdef __CUDACC__
__global__ void factor_kernel(int count,int n,const double *input,const double *q,const int64_t *offset,double *out){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<count)factor_one(i,n,input,q,offset,out+FWIDTH*i);}
#endif
int main(int argc,char **argv){
    if(argc!=5)return 1;const std::string root=argv[1],deck=argv[2],prefix=argv[3];
    auto input=read<double>(root+"/"+deck+"_input.bin"),q=read<double>(root+"/"+deck+"_Q.bin");
    auto offset=read<int64_t>(root+"/offsets.bin");int n=offset.size()-1,count=input.size()/36;
    if(input.size()!=36*size_t(count)||count!=2*n||offset[0]!=0||2*offset[n]!=int64_t(q.size()))return 2;
    std::vector<double> out(count*FWIDTH);
    if(std::strcmp(argv[4],"cpu")==0){
        #pragma omp parallel for
        for(int i=0;i<count;i++)factor_one(i,n,input.data(),q.data(),offset.data(),out.data()+FWIDTH*i);
    }else{
#ifdef __CUDACC__
        double *di=copy(input),*dq=copy(q),*dout=copy(out);int64_t *doff=copy(offset);
        factor_kernel<<<(count+63)/64,64>>>(count,n,di,dq,doff,dout);
        if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 3;
        if(cudaMemcpy(out.data(),dout,out.size()*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return 4;
        cudaFree(di);cudaFree(dq);cudaFree(dout);cudaFree(doff);
#else
        return 5;
#endif
    }
    for(int i=0;i<count;i++)if(out[FWIDTH*i+62]){std::fprintf(stderr,"map%d status%g\n",i,out[FWIDTH*i+62]);return 6;}
    write(prefix+".bin",out);std::printf("{\"passed\":true,\"maps\":%d,\"deck\":\"%s\",\"backend\":\"%s\"}\n",count,argv[2],argv[4]);return 0;
}

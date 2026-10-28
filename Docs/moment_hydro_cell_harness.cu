#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <vector>
#include "moment_kinetic_face_flux.h"
#include "moment_hydro_cell_source.h"
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define HD static __host__ __device__
#else
#define HD static
#endif
enum {WIDTH=106};
HD void node(const double *a,double *f){f[18]=lmk_face_node_flux(a,f);}
HD void cell(int i,const double *a,const double *nf,const double *raw,
             const int64_t *offset,const int64_t *inc,const int64_t *frozen,double *out){
    double dr[18]={0.},F[2]={0.},H[4]={0.},work=0.,src[18],diag[8];
    const double *r=raw+19*i;const double c[2]={r[1]/r[0],r[2]/r[0]};
    for(int64_t t=offset[i];t<offset[i+1];t++){
        const int64_t code=inc[t],k=code>0?code-1:-code-1;
        const double sign=code>0?1.:-1.;const double *v=a+33*k,*flux=nf+19*k;
        if(flux[18]!=0.){out[105]=10.+flux[18];return;}
        for(int q=0;q<18;q++)dr[q]+=sign*flux[q];
        const double traction[2]={v[4]*v[31]*v[2],v[4]*v[31]*v[3]};
        for(int p=0;p<2;p++){
            F[p]+=sign*traction[p];
            for(int q=0;q<2;q++)H[2*q+p]+=sign*(v[q]-c[q])*traction[p];
        }
        work+=sign*v[4]*v[31]*v[32];
    }
    const int status=lmh_cell_source(r,F,H,r[18],src,diag);
    if(status){out[105]=status;return;}
    double pressure[18]={0.};pressure[7]=F[0];pressure[8]=F[1];pressure[17]=work;
    for(int p=0;p<2;p++)for(int q=0;q<2;q++){
        pressure[9+2*p+q]=H[2*p+q]+(p==q?r[18]:0.)+c[p]*F[q];
        pressure[13+2*p+q]=src[13+2*p+q];
    }
    for(int q=0;q<18;q++){
        dr[q]+=src[q];out[q]=dr[q];out[18+q]=frozen[i]?0.:dr[q];
        out[36+q]=frozen[i]?0.:pressure[q];
        out[54+q]=out[18+q]-out[36+q];out[78+q]=src[q];
    }
    for(int q=0;q<2;q++)out[72+q]=F[q];
    for(int q=0;q<4;q++)out[74+q]=H[q];
    for(int q=0;q<8;q++)out[96+q]=diag[q];
    out[104]=work;out[105]=0.;
}
HD void ledger(int n,int nodes,const double *nf,const int64_t *owners,
               const int64_t *frozen,const double *out,double *result){
    for(int k=0;k<8;k++)result[k]=0.;
    for(int k=0;k<nodes;k++)if(owners[2*k+1]<0){
        const double *f=nf+19*k;
        result[0]-=f[7];result[1]-=f[8];result[2]-=f[17];result[3]-=f[10]-f[11];
    }
    for(int i=0;i<n;i++)if(frozen[i]){
        const double *r=out+WIDTH*i;
        result[4]-=r[7];result[5]-=r[8];result[6]-=r[17];result[7]-=r[10]-r[11];
    }
}
HD int fixtures(){
    double r[18]={0.},F[2]={0.},H[4]={-1.,0.,0.,-1.},src[18],diag[8];
    r[0]=1.;r[3]=r[6]=1./12.;r[17]=1.;
    if(lmh_cell_source(r,F,H,1.,src,diag))return 1;
    double summed[18]={0.};
    for(int k=0;k<4;k++){
        double a[33]={0.},flux[18];int axis=k/2;double sign=k%2?1.:-1.;
        a[axis]=.5*sign;a[2+axis]=sign;a[4]=1.;a[31]=1.;
        if(lmk_face_node_flux(a,flux))return 2;
        for(int q=0;q<18;q++)summed[q]-=flux[q];
    }
    for(int q=0;q<18;q++)if(fabs(summed[q]+src[q])>1e-15)return 3;
    for(int kind=0;kind<4;kind++){
        double bad[18],s[18],d[8];for(int q=0;q<18;q++){bad[q]=r[q];s[q]=-83.;}
        for(int q=0;q<8;q++)d[q]=-83.;
        if(kind==0)bad[0]=0.;if(kind==1)bad[6]=-1.;if(kind==2)bad[0]=NAN;
        if(kind==3)bad[14]=1.;
        if(lmh_cell_source(bad,F,H,1.,s,d)==0)return 4;
        for(int q=0;q<18;q++)if(s[q]!=-83.)return 5;
        for(int q=0;q<8;q++)if(d[q]!=-83.)return 6;
    }
    return 0;
}
#ifdef __CUDACC__
__global__ void nodes_kernel(const double *a,double *f,int n){int k=blockIdx.x*blockDim.x+threadIdx.x;if(k<n)node(a+33*k,f+19*k);}
__global__ void cells_kernel(int n,const double *a,const double *nf,const double *r,const int64_t *off,const int64_t *inc,const int64_t *fr,double *out){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){cell(i,a,nf,r,off,inc,fr,out+WIDTH*i);if(i==0)out[105]+=fixtures();}}
__global__ void ledger_kernel(int n,int nn,const double *nf,const int64_t *own,const int64_t *fr,const double *out,double *result){if(threadIdx.x==0&&blockIdx.x==0)ledger(n,nn,nf,own,fr,out,result);}
#endif
template<class T> std::vector<T> read(const std::string &p){std::ifstream f(p,std::ios::binary|std::ios::ate);if(!f)std::exit(2);auto bytes=f.tellg();if(bytes%sizeof(T))std::exit(3);std::vector<T> a(bytes/sizeof(T));f.seekg(0);f.read(reinterpret_cast<char*>(a.data()),bytes);if(!f)std::exit(4);return a;}
template<class T> void write(const std::string &p,const std::vector<T>&a){std::ofstream f(p,std::ios::binary);f.write(reinterpret_cast<const char*>(a.data()),a.size()*sizeof(T));if(!f)std::exit(5);}
#ifdef __CUDACC__
template<class T> T *copy(const std::vector<T>&a){T *p=nullptr;if(cudaMalloc(&p,a.size()*sizeof(T))!=cudaSuccess || cudaMemcpy(p,a.data(),a.size()*sizeof(T),cudaMemcpyHostToDevice)!=cudaSuccess)std::exit(6);return p;}
#endif
#ifndef LAGFORCE_HYDRO_NO_MAIN
int main(int argc,char **argv){
    if(argc!=5)return 1;std::string input=argv[1],prep=argv[2],prefix=argv[3];
    auto a=read<double>(input+"/inputs.bin"),raw=read<double>(prep+"/cells.bin");
    auto own=read<int64_t>(prep+"/owners.bin"),off=read<int64_t>(prep+"/offsets.bin"),inc=read<int64_t>(prep+"/incidence.bin"),fr=read<int64_t>(prep+"/frozen.bin");
    const int nn=a.size()/33,n=raw.size()/19;
    if(a.size()!=33*size_t(nn)||raw.size()!=19*size_t(n)||own.size()!=2*size_t(nn)||off.size()!=size_t(n+1)||fr.size()!=size_t(n)||off[0]!=0||off[n]!=int64_t(inc.size()))return 12;
    for(int i=0;i<n;i++)if(off[i]>off[i+1] || (fr[i]!=0&&fr[i]!=1))return 12;
    for(auto code:inc)if(code==0 || code>nn || code< -int64_t(nn))return 12;
    std::vector<double> nf(19*nn),out(WIDTH*n),led(8);
    if(std::strcmp(argv[4],"cpu")==0){
        #pragma omp parallel for
        for(int k=0;k<nn;k++)node(a.data()+33*k,nf.data()+19*k);
        #pragma omp parallel for
        for(int i=0;i<n;i++)cell(i,a.data(),nf.data(),raw.data(),off.data(),inc.data(),fr.data(),out.data()+WIDTH*i);
        out[105]+=fixtures();ledger(n,nn,nf.data(),own.data(),fr.data(),out.data(),led.data());
    }else{
#ifdef __CUDACC__
        double *da=copy(a),*dr=copy(raw),*dnf=copy(nf),*dout=copy(out),*dled=copy(led);
        int64_t *down=copy(own),*doff=copy(off),*dinc=copy(inc),*dfr=copy(fr);
        nodes_kernel<<<(nn+127)/128,128>>>(da,dnf,nn);
        cells_kernel<<<(n+127)/128,128>>>(n,da,dnf,dr,doff,dinc,dfr,dout);
        ledger_kernel<<<1,1>>>(n,nn,dnf,down,dfr,dout,dled);
        if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 7;
        if(cudaMemcpy(out.data(),dout,out.size()*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess||cudaMemcpy(led.data(),dled,led.size()*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return 8;
        cudaFree(da);cudaFree(dr);cudaFree(dnf);cudaFree(dout);cudaFree(dled);cudaFree(down);cudaFree(doff);cudaFree(dinc);cudaFree(dfr);
#else
        return 9;
#endif
    }
    for(int i=0;i<n;i++)if(out[WIDTH*i+105]!=0.){std::fprintf(stderr,"cell%d status%g\n",i,out[WIDTH*i+105]);return 10;}
    write(prefix+".bin",out);write(prefix+"_ledger.bin",led);
    std::printf("{\"passed\":true,\"nodes\":%d,\"cells\":%d,\"backend\":\"%s\"}\n",nn,n,argv[4]);return 0;
}
#endif

/* Research bridge: same checked operators as the completed backend harnesses.
 * Outputs commit only after all cell/kernel status checks pass. No MPI here. */
#define LAGFORCE_HYDRO_NO_MAIN
#define LAGFORCE_STAGE_NO_MAIN
#include "moment_stage_harness.cu"
#include "moment_factor_update.h"
#include <omp.h>
#include <algorithm>

template<class T> struct Device {
    T *p=nullptr;size_t n=0;
    bool set(size_t count,const T *host=nullptr){
        n=count;if(!n)return true;
        if(cudaMalloc(&p,n*sizeof(T))!=cudaSuccess)return false;
        return !host||cudaMemcpy(p,host,n*sizeof(T),cudaMemcpyHostToDevice)==cudaSuccess;
    }
    bool get(T *host){return !n||cudaMemcpy(host,p,n*sizeof(T),cudaMemcpyDeviceToHost)==cudaSuccess;}
    ~Device(){if(p)cudaFree(p);}
};
__global__ void bridge_stage_kernel(int n,const double *a,const double *nf,const double *raw,
    const int64_t *off,const int64_t *inc,const int64_t *fr,const double *extra,const double *L,
    const double *gn,const int64_t *goff,double *out,double *qout,double *hydro){
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n)stage_cell(i,a,nf,raw,off,inc,fr,extra,L,gn,goff,out+SWIDTH*i,qout,hydro+WIDTH*i);
}
__global__ void bridge_ledger_kernel(int n,int nn,const double *nf,const int64_t *owners,
    const int64_t *fr,const double *stage,const double *hydro,double *result){
    if(threadIdx.x==0&&blockIdx.x==0){
        for(int k=0;k<12;k++)result[k]=0.;
        for(int i=0;i<n;i++)if(stage[SWIDTH*i+93]!=0.)return;
        ledger(n,nn,nf,owners,fr,hydro,result);gravity_ledger(n,fr,stage,result+8);
    }
}
/* backend0 CPU,1 CUDA. counts are exact flat-array sizes supplied by caller.
 * stage n*94,Q incidence*16,hydro n*106,ledgers12. No output overlap. */
extern "C" int lagmoment_stage(int backend,int threads,int n,int nn,int64_t ni,int64_t ng,
    const double *a,const double *raw,const int64_t *owners,const int64_t *off,
    const int64_t *inc,const int64_t *fr,const double *extra,const double *L,
    const double *gn,const int64_t *goff,double *stage_out,double *q_out,double *hydro_out,double *ledgers_out){
    if((backend!=0&&backend!=1)||threads<1||n<1||nn<1||ni<1||ng<1||!a||!raw||!owners||!off||!inc||!fr||!extra||!L||!gn||!goff||!stage_out||!q_out||!hydro_out||!ledgers_out)return 1;
    if(off[0]!=0||off[n]!=ni||goff[0]!=0||goff[n]!=ng)return 1;
    for(int i=0;i<n;i++)if(off[i]>off[i+1]||goff[i]>=goff[i+1]||(fr[i]!=0&&fr[i]!=1))return 1;
    int64_t expected_incidence=0;
    for(int i=0;i<nn;i++){
        if(owners[2*i]<0||owners[2*i]>=n||owners[2*i+1]>=n||owners[2*i+1]< -1)return 1;
        expected_incidence+=owners[2*i+1]>=0?2:1;
    }
    if(expected_incidence!=ni)return 1;
    for(int i=0;i<n;i++){
        int64_t previous=-1;
        for(int64_t t=off[i];t<off[i+1];t++){
            const int64_t code=inc[t];if(code==0||code>nn||code< -int64_t(nn))return 1;
            const int64_t k=code>0?code-1:-code-1;
            if(k<=previous||owners[2*k+(code>0?1:0)]!=i)return 1;
            previous=k;
        }
    }
    try{
        std::vector<double> nf(19*size_t(nn)),out(SWIDTH*size_t(n)),q(16*size_t(ni)),hydro(WIDTH*size_t(n)),led(12);
        if(backend==0){
            omp_set_num_threads(threads);
            #pragma omp parallel for
            for(int k=0;k<nn;k++)node(a+33*k,nf.data()+19*k);
            #pragma omp parallel for
            for(int i=0;i<n;i++)stage_cell(i,a,nf.data(),raw,off,inc,fr,extra,L,gn,goff,out.data()+SWIDTH*i,q.data(),hydro.data()+WIDTH*i);
            for(int i=0;i<n;i++)if(out[SWIDTH*i+93]!=0.)return 100+int(out[SWIDTH*i+93]);
            ledger(n,nn,nf.data(),owners,fr,hydro.data(),led.data());gravity_ledger(n,fr,out.data(),led.data()+8);
        }else{
            Device<double> da,dr,de,dl,dgn,dnf,dout,dq,dh,ds;
            Device<int64_t> down,doff,dinc,dfr,dgoff;
            if(!da.set(33*size_t(nn),a)||!dr.set(19*size_t(n),raw)||!de.set(20*size_t(nn),extra)||!dl.set(4*size_t(n),L)||!dgn.set(7*size_t(ng),gn)||!dnf.set(nf.size())||!dout.set(out.size())||!dq.set(q.size())||!dh.set(hydro.size())||!ds.set(12)||!down.set(2*size_t(nn),owners)||!doff.set(n+1,off)||!dinc.set(ni,inc)||!dfr.set(n,fr)||!dgoff.set(n+1,goff))return 2;
            nodes_kernel<<<(nn+127)/128,128>>>(da.p,dnf.p,nn);
            bridge_stage_kernel<<<(n+63)/64,64>>>(n,da.p,dnf.p,dr.p,doff.p,dinc.p,dfr.p,de.p,dl.p,dgn.p,dgoff.p,dout.p,dq.p,dh.p);
            bridge_ledger_kernel<<<1,1>>>(n,nn,dnf.p,down.p,dfr.p,dout.p,dh.p,ds.p);
            if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 3;
            if(!dout.get(out.data())||!dq.get(q.data())||!dh.get(hydro.data())||!ds.get(led.data()))return 4;
        }
        for(int i=0;i<n;i++)if(out[SWIDTH*i+93]!=0.)return 100+int(out[SWIDTH*i+93]);
        for(double v:out)if(!isfinite(v))return 5;
        for(double v:q)if(!isfinite(v))return 5;
        for(double v:hydro)if(!isfinite(v))return 5;
        for(double v:led)if(!isfinite(v))return 5;
        std::copy(out.begin(),out.end(),stage_out);std::copy(q.begin(),q.end(),q_out);
        std::copy(hydro.begin(),hydro.end(),hydro_out);std::copy(led.begin(),led.end(),ledgers_out);return 0;
    }catch(...){return 6;}
}
HD void bridge_factor_one(int i,const double *base,const double *linear,const double *A,
    const double *K,const double *q,const int64_t *off,double h,const double *raw,
    const int64_t *fr,double *out){
    const double *qi=q?q+2*off[i]:nullptr;
    out[62]=lmf_update(base+4*i,linear+4*i,A+4*i,K+4*i,qi,int(off[i+1]-off[i]),h,raw+18*i,int(fr[i]),out,out+4,out+22);
}
__global__ void bridge_factor_kernel(int n,const double *base,const double *linear,
    const double *A,const double *K,const double *q,const int64_t *off,double h,
    const double *raw,const int64_t *fr,double *out){
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n)bridge_factor_one(i,base,linear,A,K,q,off,h,raw,fr,out+63*i);
}
extern "C" int lagmoment_factor(int backend,int threads,int n,int64_t columns,
    const double *base,const double *linear,const double *A,const double *K,const double *q,
    const int64_t *off,double h,const double *raw,const int64_t *fr,double *output){
    if((backend!=0&&backend!=1)||threads<1||n<1||columns<0||!base||!linear||!A||!K||(!q&&columns)||!off||!raw||!fr||!output||!isfinite(h)||h<0.)return 1;
    if(off[0]!=0||off[n]!=columns)return 1;
    for(int i=0;i<n;i++)if(off[i]>off[i+1]||off[i+1]-off[i]>2147483647||(fr[i]!=0&&fr[i]!=1))return 1;
    try{
        std::vector<double> out(63*size_t(n));
        if(backend==0){
            omp_set_num_threads(threads);
            #pragma omp parallel for
            for(int i=0;i<n;i++)bridge_factor_one(i,base,linear,A,K,q,off,h,raw,fr,out.data()+63*i);
        }else{
            Device<double> db,dl,da,dk,dq,dr,dout;Device<int64_t> doff,dfr;
            if(!db.set(4*size_t(n),base)||!dl.set(4*size_t(n),linear)||!da.set(4*size_t(n),A)||!dk.set(4*size_t(n),K)||!dq.set(2*size_t(columns),q)||!dr.set(18*size_t(n),raw)||!dout.set(out.size())||!doff.set(n+1,off)||!dfr.set(n,fr))return 2;
            bridge_factor_kernel<<<(n+63)/64,64>>>(n,db.p,dl.p,da.p,dk.p,dq.p,doff.p,h,dr.p,dfr.p,dout.p);
            if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 3;
            if(!dout.get(out.data()))return 4;
        }
        for(int i=0;i<n;i++)if(out[63*i+62]!=0.)return 100+int(out[63*i+62]);
        for(double v:out)if(!isfinite(v))return 5;
        std::copy(out.begin(),out.end(),output);return 0;
    }catch(...){return 6;}
}

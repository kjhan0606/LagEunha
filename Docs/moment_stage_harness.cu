#define main hydro_reference_main
#include "moment_hydro_cell_harness.cu"
#undef main
#include "moment_covariance_stage.h"
enum {SWIDTH=94};
HD int stage_fixtures(){
    double nodeinput[7]={.7,-.2,.3,.4,.5,-.7,.2},g[18];
    if(lmc_gravity_node(nodeinput,g))return 1;
    if(fabs(g[17]-.5*(g[13]+g[16]))>1e-15||fabs(g[10]-g[11])>1e-15)return 2;
    double raw[18]={0.},L[4]={1.,.2,.3,.8},B[4],G[4],Ag[4],dg[3];
    raw[0]=1.;raw[3]=1.;raw[6]=2.;raw[9]=.2;raw[10]=.1;raw[11]=-.1;raw[12]=.3;
    double src[18]={0.};src[13]=.4;src[14]=src[15]=-.3;src[16]=.2;
    if(lmc_gravity_covariance(raw,src,L,B,G,Ag,dg))return 3;
    const double R[4]={1.04,.46,.46,.73};
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        double v=0.;for(int k=0;k<2;k++)v+=Ag[2*i+k]*R[2*k+j]+R[2*i+k]*Ag[2*k+j];
        if(fabs(v-G[2*i+j])>4e-15)return 4;
    }
    // Exact affine velocity field: gravity residual covariance rate is zero.
    const double dC[4]={.3,-.2,.1,.4};
    for(int k=0;k<4;k++)src[9+k]=dC[k];
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        src[13+2*i+j]=B[2*i]*dC[j]+B[2*i+1]*dC[2+j]+dC[i]*B[2*j]+dC[2+i]*B[2*j+1];
    }
    src[15]=src[14];
    if(lmc_gravity_covariance(raw,src,L,B,G,Ag,dg))return 5;
    for(int k=0;k<4;k++)if(fabs(G[k])>1e-15)return 6;
    for(int k=0;k<18;k++)src[k]=0.;
    if(lmc_gravity_covariance(raw,src,L,B,G,Ag,dg))return 7;
    for(int k=0;k<4;k++)if(G[k]!=0.||Ag[k]!=0.)return 8;
    for(int kind=0;kind<2;kind++){
        double badL[4]={1.,1.,1.,1.},badg[18]={0.},bo[4],go[4],ao[4],dd[3];
        for(int k=0;k<4;k++)bo[k]=go[k]=ao[k]=-73.;for(int k=0;k<3;k++)dd[k]=-73.;
        if(kind==1){for(int k=0;k<4;k++)badL[k]=L[k];badg[14]=1.;}
        if(lmc_gravity_covariance(raw,badg,badL,bo,go,ao,dd)==0)return 9;
        for(int k=0;k<4;k++)if(bo[k]!=-73.||go[k]!=-73.||ao[k]!=-73.)return 10;
        for(int k=0;k<3;k++)if(dd[k]!=-73.)return 11;
    }
    double a[33]={0.},extra[20]={0.},K[4],Q[8];a[4]=1.;
    if(lmc_transport_node(a,extra,0,1.,K,Q))return 12;
    for(int k=0;k<4;k++)if(K[k]!=0.)return 13;
    for(int k=0;k<8;k++)if(Q[k]!=0.)return 14;
    return 0;
}
HD void stage_cell(int i,const double *a,const double *nf,const double *raw,const int64_t *off,
                  const int64_t *inc,const int64_t *fr,const double *extra,const double *L,
                  const double *gn,const int64_t *goff,double *out,double *qout,double *hydro_capture=nullptr){
    double hydro[WIDTH]={0.};cell(i,a,nf,raw,off,inc,fr,hydro);
    if(hydro[105]){out[93]=10+hydro[105];return;}
    double K[4]={0.},incoming[4]={0.},grav[18]={0.};
    for(int64_t t=off[i];t<off[i+1];t++){
        const int64_t code=inc[t],k=code>0?code-1:-code-1;double kn[4],q[8];
        int status=lmc_transport_node(a+33*k,extra+20*k,code>0?1:0,raw[19*i],kn,q);
        if(status){out[93]=20+status;return;}
        for(int p=0;p<4;p++)K[p]+=kn[p];
        for(int s=0;s<4;s++){
            for(int p=0;p<2;p++)for(int d=0;d<2;d++)incoming[2*p+d]+=q[2*s+p]*q[2*s+d];
            for(int p=0;p<2;p++){qout[16*t+2*s+p]=q[2*s+p];qout[16*t+8+2*s+p]=fr[i]?0.:q[2*s+p];}
        }
    }
    for(int64_t t=goff[i];t<goff[i+1];t++){
        double g[18];int status=lmc_gravity_node(gn+7*t,g);if(status){out[93]=30+status;return;}
        for(int p=0;p<18;p++)grav[p]+=g[p];
    }
    double B[4],G[4],Ag[4],diag[3];int status=lmc_gravity_covariance(raw+19*i,grav,L+4*i,B,G,Ag,diag);
    if(status){out[93]=40+status;return;}
    for(int p=0;p<18;p++){out[p]=fr[i]?0.:hydro[18+p]+grav[p];out[30+p]=grav[p];out[75+p]=hydro[p]+grav[p];}
    for(int p=0;p<4;p++){
        out[18+p]=B[p];out[22+p]=fr[i]?0.:K[p];out[26+p]=fr[i]?0.:incoming[p];
        out[48+p]=G[p];out[52+p]=Ag[p];out[56+p]=fr[i]?0.:-B[p]+Ag[p];
        out[63+p]=K[p];out[67+p]=incoming[p];out[71+p]=-B[p]+Ag[p];
    }
    for(int p=0;p<3;p++)out[60+p]=diag[p];out[93]=0.;
    if(hydro_capture)for(int p=0;p<WIDTH;p++)hydro_capture[p]=hydro[p];
}
HD void gravity_ledger(int n,const int64_t *fr,const double *out,double *led){
    for(int k=0;k<4;k++)led[k]=0.;
    for(int i=0;i<n;i++)if(!fr[i]){const double *g=out+SWIDTH*i+30;led[0]+=g[7];led[1]+=g[8];led[2]+=g[17];led[3]+=g[10]-g[11];}
}
#ifdef __CUDACC__
__global__ void stage_kernel(int n,const double *a,const double *nf,const double *raw,const int64_t *off,const int64_t *inc,const int64_t *fr,const double *extra,const double *L,const double *gn,const int64_t *goff,double *out,double *qout){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){stage_cell(i,a,nf,raw,off,inc,fr,extra,L,gn,goff,out+SWIDTH*i,qout);if(i==0)out[93]+=stage_fixtures();}}
__global__ void gravity_ledger_kernel(int n,const int64_t *fr,const double *out,double *led){if(threadIdx.x==0&&blockIdx.x==0)gravity_ledger(n,fr,out,led);}
#endif
#ifndef LAGFORCE_STAGE_NO_MAIN
int main(int argc,char **argv){
    if(argc!=6)return 1;std::string capture=argv[1],prep=argv[2],stage=argv[3],prefix=argv[4];
    auto a=read<double>(capture+"/inputs.bin"),raw=read<double>(prep+"/cells.bin"),extra=read<double>(stage+"/node_extra.bin"),L=read<double>(stage+"/factors.bin"),gn=read<double>(stage+"/gravity_nodes.bin");
    auto off=read<int64_t>(prep+"/offsets.bin"),inc=read<int64_t>(prep+"/incidence.bin"),fr=read<int64_t>(prep+"/frozen.bin"),goff=read<int64_t>(stage+"/gravity_offsets.bin");
    const int n=raw.size()/19,nn=a.size()/33;
    if(extra.size()!=20*size_t(nn)||L.size()!=4*size_t(n)||goff.size()!=size_t(n+1)||goff[n]*7!=int64_t(gn.size()))return 2;
    std::vector<double> nf(19*nn),out(SWIDTH*n),qout(16*inc.size()),led(4);
    if(std::strcmp(argv[5],"cpu")==0){
        #pragma omp parallel for
        for(int k=0;k<nn;k++)node(a.data()+33*k,nf.data()+19*k);
        #pragma omp parallel for
        for(int i=0;i<n;i++)stage_cell(i,a.data(),nf.data(),raw.data(),off.data(),inc.data(),fr.data(),extra.data(),L.data(),gn.data(),goff.data(),out.data()+SWIDTH*i,qout.data());
        out[93]+=stage_fixtures();gravity_ledger(n,fr.data(),out.data(),led.data());
    }else{
#ifdef __CUDACC__
        double *da=copy(a),*dr=copy(raw),*de=copy(extra),*dl=copy(L),*dgn=copy(gn),*dnf=copy(nf),*dout=copy(out),*dqout=copy(qout),*dled=copy(led);
        int64_t *doff=copy(off),*dinc=copy(inc),*dfr=copy(fr),*dgoff=copy(goff);
        nodes_kernel<<<(nn+127)/128,128>>>(da,dnf,nn);
        stage_kernel<<<(n+63)/64,64>>>(n,da,dnf,dr,doff,dinc,dfr,de,dl,dgn,dgoff,dout,dqout);
        gravity_ledger_kernel<<<1,1>>>(n,dfr,dout,dled);
        if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 3;
        if(cudaMemcpy(out.data(),dout,out.size()*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess||cudaMemcpy(qout.data(),dqout,qout.size()*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess||cudaMemcpy(led.data(),dled,led.size()*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return 4;
        cudaFree(da);cudaFree(dr);cudaFree(de);cudaFree(dl);cudaFree(dgn);cudaFree(dnf);cudaFree(dout);cudaFree(dqout);cudaFree(dled);cudaFree(doff);cudaFree(dinc);cudaFree(dfr);cudaFree(dgoff);
#else
        return 5;
#endif
    }
    for(int i=0;i<n;i++)if(out[SWIDTH*i+93]){std::fprintf(stderr,"cell%d status%g\n",i,out[SWIDTH*i+93]);return 6;}
    write(prefix+".bin",out);write(prefix+"_Q.bin",qout);write(prefix+"_gravity_ledger.bin",led);
    std::printf("{\"passed\":true,\"cells\":%d,\"nodes\":%d,\"backend\":\"%s\"}\n",n,nn,argv[5]);return 0;
}
#endif

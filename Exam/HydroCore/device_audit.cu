#include <cstdio>
#include <cmath>
#include <cuda_runtime.h>
#include "../gfs_variational_core.h"
#include "../gfs_velocity_fit.h"
struct Result {double fit_error,heat_error,adiabat_error,shear_pressure;int valid;};
__global__ void audit(Result *out){
    int d=threadIdx.x+1;if(d>3)return;
    GfsVelocityFit fit;gfs_velocity_fit_init(&fit,d);
    double A[3][3]={{0,1.5,-.1},{-.7,0,.6},{.8,-.3,0}};
    for(int j=1;j<=12;j++){
        double r[3]={0},dv[3]={0};
        for(int a=0;a<d;a++)r[a]=sin(.7*j*(a+1))+.2*cos(.37*j*(a+2));
        for(int a=0;a<d;a++)for(int b=0;b<d;b++)dv[a]+=A[a][b]*r[b];
        if(!gfs_velocity_fit_add_face(&fit,r,dv,.2+j*.17)){out[d-1].valid=0;return;}
    }
    double G[3][3],condition;
    if(!gfs_velocity_fit_solve(&fit,G,&condition)){out[d-1].valid=0;return;}
    Result q={0,0,0,0,1};
    for(int a=0;a<d;a++)for(int b=0;b<d;b++)q.fit_error=fmax(q.fit_error,fabs(G[a][b]-A[a][b]));
    double r[3]={.7,d>1?.4:0,d>2?-.3:0},c[3],vi[3]={0},vj[3]={0},res[3]={0};
    for(int a=0;a<3;a++)c[a]=.5*r[a];
    for(int a=0;a<d;a++){
        for(int b=0;b<d;b++){vj[a]+=A[a][b]*r[b];res[a]-=G[a][b]*r[b];}
        res[a]+=vj[a];
    }
    GfsVariationalFace f;
    if(!gfs_variational_face(d,r,c,1,&f)){out[d-1].valid=0;return;}
    double trace=0;for(int a=0;a<d;a++)trace+=G[a][a];
    double sf=gfs_variational_compression_fraction(trace,1,trace,1);
    q.shear_pressure=gfs_variational_pair_pressure(&f,vi,vj,res,1,1,.125,sqrt(8.),sf,1,1);
    for(int a=0;a<d;a++){vi[a]=f.normal[a];vj[a]=-f.normal[a];res[a]=vj[a]-vi[a];}
    double p=gfs_variational_pair_pressure(&f,vi,vj,res,1,1,.125,sqrt(8.),1,1,1);
    double work=0,heat=0;
    for(int a=0;a<d;a++){
        work+=(-p*f.normal[a])*(vi[a]-vj[a]);
        heat+=p*f.normal[a]*(vi[a]-vj[a]);
    }
    q.heat_error=fabs(work+heat);
    double U;
    if(!gfs_variational_adiabatic(2,.3,.4,1.4,&U)){out[d-1].valid=0;return;}
    q.adiabat_error=fabs(U*pow(.4,.4)-2*pow(.3,.4));
    out[d-1]=q;
}
int main(){
    Result *device,host[3];
    cudaError_t e=cudaMalloc(&device,3*sizeof(Result));
    if(e!=cudaSuccess){fprintf(stderr,"%s\n",cudaGetErrorString(e));return 2;}
    audit<<<1,3>>>(device);
    e=cudaMemcpy(host,device,sizeof(host),cudaMemcpyDeviceToHost);
    if(e!=cudaSuccess){fprintf(stderr,"%s\n",cudaGetErrorString(e));return 2;}
    int pass=1;printf("{\"dimensions\":[");
    for(int i=0;i<3;i++){
        Result r=host[i];
        int ok=r.valid&&r.fit_error<1e-13&&r.heat_error<1e-13&&r.adiabat_error<1e-13&&r.shear_pressure<1e-13;
        pass&=ok;
        printf("%s{\"dim\":%d,\"fit_error\":%.17g,\"heat_error\":%.17g,\"adiabat_error\":%.17g,\"affine_shear_pressure\":%.17g,\"pass\":%s}",
            i?",":"",i+1,r.fit_error,r.heat_error,r.adiabat_error,r.shear_pressure,ok?"true":"false");
    }
    printf("],\"pass\":%s}\n",pass?"true":"false");cudaFree(device);return pass?0:1;
}

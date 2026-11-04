/* Same saved native geometry and velocity fit. No time evolution. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "../gfs_velocity_fit.h"
#include "../gfs_variational_core.h"
typedef struct {uint64_t magic,n,nf;double gamma;} Header;
typedef struct {int64_t id;double v[12];} Cell;
typedef struct {int i,j;double ax,ay,bx,by;} Face;
_Static_assert(sizeof(Header)==32&&sizeof(Cell)==104&&sizeof(Face)==40,"Audit ABI");
int main(int argc,char **argv){
    if(argc<2||argc>3)return 2;
    const char *mode=argc==3?argv[2]:"actual";
    int weighted=!strcmp(mode,"actual_area_fit")||!strcmp(mode,"circular_area_fit");
    int circular=!strcmp(mode,"circular_lsq")||!strcmp(mode,"circular_exact")||!strcmp(mode,"circular_area_fit");
    int exact=!strcmp(mode,"circular_exact");
    if(!circular&&!weighted&&strcmp(mode,"actual"))return 2;
    FILE *in=fopen(argv[1],"rb");Header h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=UINT64_C(0x554e494632444131)||h.n>2000000||h.nf>40000000)return 2;
    Cell *c=calloc(h.n,sizeof(*c));Face *f=calloc(h.nf,sizeof(*f));double *heat=calloc(5*h.n,sizeof(double));
    if(!c||!f||!heat||fread(c,sizeof(*c),h.n,in)!=h.n||fread(f,sizeof(*f),h.nf,in)!=h.nf)return 2;
    fclose(in);
    int *offset=calloc(h.n+1,sizeof(int));double *condition=calloc(h.n,sizeof(double));
    if(!offset||!condition)return 2;
    for(uint64_t k=0;k<h.nf;k++){
        if(f[k].i<0||f[k].i>=(int)h.n||(k&&f[k].i<f[k-1].i))return 2;
        offset[f[k].i+1]++;
    }
    for(uint64_t i=0;i<h.n;i++){
        offset[i+1]+=offset[i];
        if(circular){double x=c[i].v[3]-3,y=c[i].v[4]-3,w=pow(x*x+y*y+.01175*.01175,-.75);c[i].v[5]=-y*w;c[i].v[6]=x*w;}
    }
    for(uint64_t i=0;i<h.n;i++){
        GfsVelocityFit fit;gfs_velocity_fit_init(&fit,2);double G[3][3];
        for(int k=offset[i];k<offset[i+1];k++){
            int j=f[k].j;double r[3]={0},dv[3]={0};double *a=c[i].v;
            if(j>=0){
                double *b=c[j].v;r[0]=b[3]-a[3];r[1]=b[4]-a[4];dv[0]=b[5]-a[5];dv[1]=b[6]-a[6];
            }else if(j<=-1000){
                double sx=f[k].by-f[k].ay,sy=f[k].ax-f[k].bx,A=hypot(sx,sy),nx=sx/A,ny=sy/A;
                double d=.5*((f[k].ax+f[k].bx)*nx+(f[k].ay+f[k].by)*ny),vn=a[5]*nx+a[6]*ny;
                r[0]=2*d*nx;r[1]=2*d*ny;dv[0]=-2*vn*nx;dv[1]=-2*vn*ny;
            }else return 2;
            if(weighted){
                /* Diagnostic alternative only. Surface/distance weighting
                 * makes a vanishing interface contribute continuously. */
                double A=hypot(f[k].by-f[k].ay,f[k].ax-f[k].bx),w=A/hypot(r[0],r[1]);
                GfsVelocityFit sample;gfs_velocity_fit_init(&sample,2);
                if(!gfs_velocity_fit_add(&sample,r,dv))return 2;
                for(int a=0;a<2;a++)for(int b=0;b<2;b++){fit.M[a][b]+=w*sample.M[a][b];fit.R[a][b]+=w*sample.R[a][b];}
                fit.count++;
            }else if(!gfs_velocity_fit_add(&fit,r,dv))return 2;
        }
        if(!gfs_velocity_fit_solve(&fit,G,condition+i))return 2;
        if(circular||weighted){c[i].v[7]=G[0][0];c[i].v[8]=G[0][1];c[i].v[9]=G[1][0];c[i].v[10]=G[1][1];}
        if(exact){
            double x=c[i].v[3]-3,y=c[i].v[4]-3,t=x*x+y*y+.01175*.01175,w=pow(t,-.75),d=1.5*pow(t,-1.75);
            c[i].v[7]=d*x*y;c[i].v[8]=-w+d*y*y;c[i].v[9]=w-d*x*x;c[i].v[10]=-c[i].v[7];
        }
    }
    long double sum[3][5]={{0}};double max_condition[3]={0};int edges=0,residual_dominated=0;
    for(uint64_t k=0;k<h.nf;k++){
        int i=f[k].i,j=f[k].j;if(i<0||i>=(int)h.n||j>=(int)h.n)return 2;
        if(j>=0&&c[i].id>c[j].id)continue;
        double *a=c[i].v,vi[3]={a[5],a[6],0};
        double sx=f[k].by-f[k].ay,sy=f[k].ax-f[k].bx,A=hypot(sx,sy);
        double centre[3]={.5*(f[k].ax+f[k].bx),.5*(f[k].ay+f[k].by),0};
        double rhoi=a[0]/a[2],Pi=(h.gamma-1)*a[1]/a[2],ci=sqrt(h.gamma*Pi/rhoi);
        GfsVariationalFace face;double actual=0,q[2]={0};
        if(j>=0){
            double *b=c[j].v,r[3]={b[3]-a[3],b[4]-a[4],0},vj[3]={b[5],b[6],0};
            double jump[3]={vj[0]-vi[0]-.5*((a[7]+b[7])*r[0]+(a[8]+b[8])*r[1]),
                            vj[1]-vi[1]-.5*((a[9]+b[9])*r[0]+(a[10]+b[10])*r[1]),0};
            double zero[3]={0},rhoj=b[0]/b[2],Pj=(h.gamma-1)*b[1]/b[2],cj=sqrt(h.gamma*Pj/rhoj);
            if(!gfs_variational_face(2,r,centre,A,&face))return 2;
            double sf=gfs_variational_compression_fraction(a[7]+a[10],a[9]-a[8],b[7]+b[10],b[9]-b[8]);
            q[0]=gfs_variational_pair_pressure(&face,vi,vj,jump,rhoi,ci,rhoj,cj,sf,sqrt(a[2]),sqrt(b[2]));
            q[1]=gfs_variational_pair_pressure(&face,vi,vj,zero,rhoi,ci,rhoj,cj,sf,sqrt(a[2]),sqrt(b[2]));
            actual=(vi[0]-vj[0])*face.normal[0]+(vi[1]-vj[1])*face.normal[1];
        }else if(j<=-1000){
            double nx=sx/A,ny=sy/A,distance=centre[0]*nx+centre[1]*ny,vn=vi[0]*nx+vi[1]*ny;
            double derivative=nx*(a[7]*nx+a[8]*ny)+ny*(a[9]*nx+a[10]*ny);
            double r[3]={2*distance*nx,2*distance*ny,0},vj[3]={vi[0]-2*vn*nx,vi[1]-2*vn*ny,0};
            double jump[3]={-2*(vn+distance*derivative)*nx,-2*(vn+distance*derivative)*ny,0},zero[3]={0};
            if(!gfs_variational_face(2,r,centre,A,&face))return 2;
            double sf=gfs_variational_compression_fraction(a[7]+a[10],a[9]-a[8],a[7]+a[10],a[9]-a[8]);
            q[0]=gfs_variational_pair_pressure(&face,vi,vj,jump,rhoi,ci,rhoi,ci,sf,sqrt(a[2]),sqrt(a[2]));
            q[1]=gfs_variational_pair_pressure(&face,vi,vj,zero,rhoi,ci,rhoi,ci,sf,sqrt(a[2]),sqrt(a[2]));
            actual=vn;
        }else return 2;
        if(q[1]>q[0]*(1+1e-14)||q[1]<0)return 2;
        if(q[0]>0){edges++;if(q[0]>q[1])residual_dominated++;}
        double li=sqrt(c[i].v[2]),lj=j>=0?sqrt(c[j].v[2]):li;
        double scale=(.5*li+.5*lj)/face.separation;
        double mu=scale*fmax(actual,0)*(j>=0?1:2);
        double cj=j>=0?sqrt(h.gamma*(h.gamma-1)*c[j].v[1]/c[j].v[0]):ci;
        double fraction=(.5*ci+.5*cj)/(.5*ci+.5*cj+2*mu);
        for(int a=0;a<5;a++){
            double value=a<2?q[a]:a==2?q[0]*fraction:a==3?q[0]/scale:q[0]*(j>=0?.5*(condition[i]+condition[j]):condition[i]);
            double Q=value*A*fmax(actual,0);
            heat[a*h.n+i]+=j>=0?.5*Q:Q;
            if(j>=0)heat[a*h.n+j]+=.5*Q;
        }
    }
    for(uint64_t i=0;i<h.n;i++){
        int bin=c[i].v[11]<.5?0:c[i].v[11]<=2?1:2;
        for(int a=0;a<5;a++)sum[bin][a]+=heat[a*h.n+i];
        max_condition[bin]=fmax(max_condition[bin],condition[i]);
    }
    printf("{\"scope\":\"instantaneous native RHS channel audit only, no evolution\",\"velocity_gradient_mode\":\"%s\",\"active_pairs\":%d,\"pairs_with_extra_reconstruction_viscosity\":%d,\"cohorts\":[",mode,edges,residual_dominated);
    const char *names[]={"original_core","original_disk","original_ambient"};
    for(int b=0;b<3;b++)printf("%s{\"name\":\"%s\",\"total_heat_rate\":%.17Lg,\"compression_branch_rate\":%.17Lg,\"extra_reconstruction_rate\":%.17Lg,\"compression_fraction_of_heat\":%.17Lg,\"linear_viscosity_rate\":%.17Lg,\"quadratic_viscosity_rate\":%.17Lg,\"heat_weighted_d_over_h\":%.17Lg,\"heat_weighted_condition_bound\":%.17Lg,\"max_condition_bound\":%.17g}",b?",":"",names[b],sum[b][0],sum[b][1],sum[b][0]-sum[b][1],sum[b][0]>0?sum[b][1]/sum[b][0]:0.L,sum[b][2],sum[b][0]-sum[b][2],sum[b][0]>0?sum[b][3]/sum[b][0]:0.L,sum[b][0]>0?sum[b][4]/sum[b][0]:0.L,max_condition[b]);
    printf("]}\n");free(c);free(f);free(heat);free(offset);free(condition);return 0;
}

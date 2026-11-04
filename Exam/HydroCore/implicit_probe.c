/* Frozen-state sparse algebra probe. This is NOT a hydrodynamic timestep.
 * D is frozen from the native CUDA RHS. No pressure force or motion occurs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "../gfs_variational_core.h"
typedef struct {uint64_t magic,n,nf;double gamma;} Header;
typedef struct {int64_t id;double v[12];} Cell;
typedef struct {int i,j;double ax,ay,bx,by;} Face;
typedef struct {int i,j;double nx,ny,D;} Edge;
_Static_assert(sizeof(Header)==32&&sizeof(Cell)==104&&sizeof(Face)==40,"Probe input ABI");
static Cell *cells;static Edge *edges;static int *offset,*adj;
static double *inv00,*inv01,*inv11;static int n,ne;static double h;
static void *allocate(size_t count,size_t size){
    void *p=calloc(count,size);if(!p){fprintf(stderr,"allocation failed\n");exit(2);}return p;
}
static double dot(const double *a,const double *b){
    long double s=0;
#pragma omp parallel for reduction(+:s) schedule(static)
    for(int i=0;i<2*n;i++)s+=(long double)a[i]*b[i];
    return (double)s;
}
static void multiply(const double *x,double *y){
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){
        double a=cells[i].v[0]*x[2*i],b=cells[i].v[0]*x[2*i+1];
        for(int k=offset[i];k<offset[i+1];k++){
            Edge *e=edges+adj[k];int j=e->i==i?e->j:e->i;
            double dx=x[2*i]-(j>=0?x[2*j]:0),dy=x[2*i+1]-(j>=0?x[2*j+1]:0);
            double q=h*e->D*(e->nx*dx+e->ny*dy);a+=q*e->nx;b+=q*e->ny;
        }
        y[2*i]=a;y[2*i+1]=b;
    }
}
static void precondition(const double *r,double *z){
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){
        z[2*i]=inv00[i]*r[2*i]+inv01[i]*r[2*i+1];
        z[2*i+1]=inv01[i]*r[2*i]+inv11[i]*r[2*i+1];
    }
}
int main(int argc,char **argv){
    if(argc!=3)return 2;
    FILE *input=fopen(argv[1],"rb");Header head;
    if(!input||fread(&head,sizeof(head),1,input)!=1||head.magic!=UINT64_C(0x554e494632444131))return 2;
    if(head.n>2000000||head.nf>40000000)return 2;n=(int)head.n;
    h=atof(argv[2]);if(!(h>0))return 2;
    cells=allocate(n,sizeof(Cell));Face *faces=allocate(head.nf,sizeof(Face));
    if(fread(cells,sizeof(Cell),n,input)!=(size_t)n||
       fread(faces,sizeof(Face),head.nf,input)!=head.nf)return 2;
    fclose(input);edges=allocate(head.nf,sizeof(Edge));
    for(uint64_t k=0;k<head.nf;k++){
        Face *f=faces+k;int i=f->i,j=f->j;
        if(i<0||i>=n||j>=n)return 2;
        if(j>=0&&cells[i].id>cells[j].id)continue;
        double *a=cells[i].v,vi[3]={a[5],a[6],0};
        double S[3]={f->by-f->ay,f->ax-f->bx,0},area=hypot(S[0],S[1]);
        double c[3]={.5*(f->ax+f->bx),.5*(f->ay+f->by),0};
        double rhoi=a[0]/a[2],Pi=(head.gamma-1)*a[1]/a[2],ci=sqrt(head.gamma*Pi/rhoi);
        GfsVariationalFace b;double D=0;
        if(j>=0){
            double *v=cells[j].v,r[3]={v[3]-a[3],v[4]-a[4],0},vj[3]={v[5],v[6],0};
            double jump[3]={vj[0]-vi[0]-.5*((a[7]+v[7])*r[0]+(a[8]+v[8])*r[1]),
                            vj[1]-vi[1]-.5*((a[9]+v[9])*r[0]+(a[10]+v[10])*r[1]),0};
            if(!gfs_variational_face(2,r,c,area,&b))return 2;
            double rhoj=v[0]/v[2],Pj=(head.gamma-1)*v[1]/v[2],cj=sqrt(head.gamma*Pj/rhoj);
            double sf=gfs_variational_compression_fraction(a[7]+a[10],a[9]-a[8],v[7]+v[10],v[9]-v[8]);
            double q=gfs_variational_pair_pressure(&b,vi,vj,jump,rhoi,ci,rhoj,cj,sf,sqrt(a[2]),sqrt(v[2]));
            double closing=(vi[0]-vj[0])*b.normal[0]+(vi[1]-vj[1])*b.normal[1];
            if(q>0&&closing>0)D=q*area/closing;
        }else if(j<=-1000){
            double nx=S[0]/area,ny=S[1]/area,dist=c[0]*nx+c[1]*ny,vn=vi[0]*nx+vi[1]*ny;
            double derivative=nx*(a[7]*nx+a[8]*ny)+ny*(a[9]*nx+a[10]*ny);
            double r[3]={2*dist*nx,2*dist*ny,0},vg[3]={vi[0]-2*vn*nx,vi[1]-2*vn*ny,0};
            double jump[3]={-2*(vn+dist*derivative)*nx,-2*(vn+dist*derivative)*ny,0};
            if(!gfs_variational_face(2,r,c,area,&b))return 2;
            double sf=gfs_variational_compression_fraction(a[7]+a[10],a[9]-a[8],a[7]+a[10],a[9]-a[8]);
            double q=gfs_variational_pair_pressure(&b,vi,vg,jump,rhoi,ci,rhoi,ci,sf,sqrt(a[2]),sqrt(a[2]));
            if(q>0&&vn>0)D=q*area/vn;
            j=-1;
        }else return 2;
        if(D>0){if(!isfinite(D))return 2;edges[ne++]=(Edge){i,j,b.normal[0],b.normal[1],D};}
    }
    free(faces);
    offset=allocate(n+1,sizeof(int));
    for(int k=0;k<ne;k++){offset[edges[k].i+1]++;if(edges[k].j>=0)offset[edges[k].j+1]++;}
    for(int i=0;i<n;i++)offset[i+1]+=offset[i];
    adj=allocate(offset[n],sizeof(int));int *cursor=allocate(n,sizeof(int));memcpy(cursor,offset,n*sizeof(int));
    for(int k=0;k<ne;k++){adj[cursor[edges[k].i]++]=k;if(edges[k].j>=0)adj[cursor[edges[k].j]++]=k;}
    free(cursor);
    inv00=allocate(n,sizeof(double));inv01=allocate(n,sizeof(double));inv11=allocate(n,sizeof(double));
    for(int i=0;i<n;i++){
        double xx=cells[i].v[0],xy=0,yy=xx;
        for(int k=offset[i];k<offset[i+1];k++){
            Edge *e=edges+adj[k];xx+=h*e->D*e->nx*e->nx;xy+=h*e->D*e->nx*e->ny;yy+=h*e->D*e->ny*e->ny;
        }
        double det=xx*yy-xy*xy;if(!(det>0))return 2;
        inv00[i]=yy/det;inv01[i]=-xy/det;inv11[i]=xx/det;
    }
    double *v=allocate(2*n,sizeof(double)),*bvec=allocate(2*n,sizeof(double)),*r=allocate(2*n,sizeof(double));
    double *z=allocate(2*n,sizeof(double)),*p=allocate(2*n,sizeof(double)),*Ap=allocate(2*n,sizeof(double));
    for(int i=0;i<n;i++)for(int a=0;a<2;a++){v[2*i+a]=cells[i].v[5+a];bvec[2*i+a]=cells[i].v[0]*v[2*i+a];}
    multiply(v,Ap);for(int k=0;k<2*n;k++)r[k]=bvec[k]-Ap[k];
    precondition(r,z);memcpy(p,z,2*n*sizeof(double));double rz=dot(r,z),normb=sqrt(dot(bvec,bvec));
    double residual=sqrt(dot(r,r))/normb,start=omp_get_wtime();int iterations=0;
    while(residual>1e-12&&iterations<2000){
        multiply(p,Ap);double pap=dot(p,Ap);if(!(pap>0))return 2;
        double alpha=rz/pap;
#pragma omp parallel for schedule(static)
        for(int k=0;k<2*n;k++){v[k]+=alpha*p[k];r[k]-=alpha*Ap[k];}
        precondition(r,z);double next=dot(r,z),beta=next/rz;
#pragma omp parallel for schedule(static)
        for(int k=0;k<2*n;k++)p[k]=z[k]+beta*p[k];
        rz=next;residual=sqrt(dot(r,r))/normb;iterations++;
    }
    double seconds=omp_get_wtime()-start;
    multiply(v,Ap);for(int k=0;k<2*n;k++)r[k]=bvec[k]-Ap[k];
    double true_residual=sqrt(dot(r,r))/normb;
    double *heat=allocate(n,sizeof(double));long double K0=0,K1=0,H=0,oldheat=0,diskH=0;
    double max_heat_over_U=0,cohort_delta=0,cohort_max=0;int nc=0,flips=0;
    double flip_impulse=0,total_impulse=0;
    for(int i=0;i<n;i++){
        double *a=cells[i].v,dx=v[2*i]-a[5],dy=v[2*i+1]-a[6];
        heat[i]=.5*a[0]*(dx*dx+dy*dy);
        K0+=.5L*a[0]*(a[5]*a[5]+a[6]*a[6]);K1+=.5L*a[0]*(v[2*i]*v[2*i]+v[2*i+1]*v[2*i+1]);
        double x=a[3]-3,y=a[4]-3,R=hypot(x,y);
        if(R>=.8&&R<=1.2){
            double vk=R/pow(R*R+.01175*.01175,.75);
            double dvphi=fabs(x*dy-y*dx)/(R*vk);cohort_delta+=dvphi;cohort_max=fmax(cohort_max,dvphi);nc++;
        }
    }
    for(int k=0;k<ne;k++){
        Edge *e=edges+k;int i=e->i,j=e->j;
        double vn=e->nx*(v[2*i]-(j>=0?v[2*j]:0))+e->ny*(v[2*i+1]-(j>=0?v[2*j+1]:0));
        double old=e->nx*(cells[i].v[5]-(j>=0?cells[j].v[5]:0))+e->ny*(cells[i].v[6]-(j>=0?cells[j].v[6]:0));
        oldheat+=e->D*old*old;double Q=h*e->D*vn*vn;
        heat[i]+=j>=0?.5*Q:Q;if(j>=0)heat[j]+=.5*Q;
        total_impulse+=h*e->D*fabs(vn);
        if(vn<0){flips++;flip_impulse+=h*e->D*fabs(vn);}
    }
    for(int i=0;i<n;i++){
        H+=heat[i];max_heat_over_U=fmax(max_heat_over_U,heat[i]/cells[i].v[1]);
        if(cells[i].v[11]>=.5&&cells[i].v[11]<=2)diskH+=heat[i];
    }
    printf("{\"scope\":\"frozen-coefficient viscosity algebra only, no pressure or position update\",\"cells\":%d,\"active_edges\":%d,\"viscosity_interval\":%.17g,\"iterations\":%d,\"solve_seconds\":%.17g,\"relative_residual\":%.17g,\"true_relative_residual\":%.17g,\"initial_heat_rate\":%.17Lg,\"finite_heat\":%.17Lg,\"disk_heat\":%.17Lg,\"kinetic_loss\":%.17Lg,\"energy_identity_residual\":%.17Lg,\"max_heat_over_U\":%.17g,\"cohort_L1_velocity_change\":%.17g,\"cohort_max_velocity_change\":%.17g,\"final_expanding_edges\":%d,\"expanding_impulse_fraction\":%.17g,\"converged\":%s}\n",
        n,ne,h,iterations,seconds,residual,true_residual,oldheat,H,diskH,K0-K1,K1-K0+H,max_heat_over_U,cohort_delta/nc,cohort_max,flips,total_impulse>0?flip_impulse/total_impulse:0,true_residual<1e-10?"true":"false");
    return true_residual<1e-10?0:1;
}

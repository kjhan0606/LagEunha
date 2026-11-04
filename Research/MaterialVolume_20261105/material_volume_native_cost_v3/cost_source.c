/* Bounded native C/OpenMP coefficient screen on one actual saved geometry.
 * Jacobi PCG, maximum 256 iterations per RHS, one total extra-step time budget.
 * No PDE, geometry projection, timestep update, or production force changes.
 * Passing cannot qualify full weighted geometry or GIZMO throughput. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "../gfs_variational_core.h"
#include "../gfs_velocity_fit.h"
typedef struct {int i,j;double rx,ry,cx,cy,area;} ProbeFace;
_Static_assert(sizeof(ProbeFace)==48,"probe face ABI");
static int n,nf;static ProbeFace *face;static int *offset;
static double *coefficient,*diagonal,budget,start;
static void laplacian(const double *x,double *out){
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){
        double sum=0;for(int k=offset[i];k<offset[i+1];k++)sum+=coefficient[k]*(x[i]-x[face[k].j]);out[i]=sum;
    }
}
static double inner(const double *x,const double *y){
    double sum=0;
#pragma omp parallel for schedule(static) reduction(+:sum)
    for(int i=0;i<n;i++)sum+=x[i]*y[i];return sum;
}
static int solve(const double *rhs,double *x,int *iterations,double *relative){
    double *r=malloc(n*sizeof(double)),*z=malloc(n*sizeof(double)),*p=malloc(n*sizeof(double)),*Ap=malloc(n*sizeof(double));
    if(!r||!z||!p||!Ap)return 2;
    double norm=inner(rhs,rhs),rz=0,mean=0;*iterations=0;*relative=1;
    if(!(norm>0))return 2;
#pragma omp parallel for schedule(static) reduction(+:mean)
    for(int i=0;i<n;i++){x[i]=0;r[i]=rhs[i];z[i]=r[i]/diagonal[i];mean+=z[i];}
    mean/=n;
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++)p[i]=z[i]-mean;
    rz=inner(r,p);int status=3;
    for(int k=0;k<256;k++){
        if(omp_get_wtime()-start>budget)break;
        laplacian(p,Ap);double curvature=inner(p,Ap);if(!(curvature>0)||!isfinite(curvature)){status=2;break;}
        double alpha=rz/curvature;
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){x[i]+=alpha*p[i];r[i]-=alpha*Ap[i];}
        *iterations=k+1;*relative=sqrt(inner(r,r)/norm);
        if(*relative<=1e-10){status=0;break;}
        mean=0;
#pragma omp parallel for schedule(static) reduction(+:mean)
        for(int i=0;i<n;i++){z[i]=r[i]/diagonal[i];mean+=z[i];}
        mean/=n;
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++)z[i]-=mean;
        double next=inner(r,z),beta=next/rz;rz=next;
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++)p[i]=z[i]+beta*p[i];
    }
    mean=0;for(int i=0;i<n;i++)mean+=x[i]/n;for(int i=0;i<n;i++)x[i]-=mean;
    if(!status){
        laplacian(x,Ap);for(int i=0;i<n;i++)Ap[i]-=rhs[i];
        *relative=sqrt(inner(Ap,Ap)/norm);if(*relative>1e-10)status=3;
    }
    free(r);free(z);free(p);free(Ap);return status;
}
int main(int argc,char **argv){
    if(argc!=3)return 2;budget=atof(argv[2]);if(!(budget>0)||!isfinite(budget))return 2;
    FILE *in=fopen(argv[1],"rb");uint64_t header[3];
    if(!in||fread(header,sizeof(header),1,in)!=1||header[0]!=UINT64_C(0x4d56434f53543131)||header[1]>1000000||header[2]>20000000)return 2;
    n=header[1];nf=header[2];double *cell=malloc(6*n*sizeof(double));face=malloc(nf*sizeof(*face));
    if(!cell||!face||fread(cell,6*sizeof(double),n,in)!=(size_t)n||fread(face,sizeof(*face),nf,in)!=(size_t)nf||fgetc(in)!=EOF||fclose(in))return 2;
    offset=calloc(n+1,sizeof(int));coefficient=calloc(nf,sizeof(double));diagonal=calloc(n,sizeof(double));
    double *rhs[2]={calloc(n,sizeof(double)),calloc(n,sizeof(double))},*psi[2]={calloc(n,sizeof(double)),calloc(n,sizeof(double))};
    double *inverse=calloc(4*n,sizeof(double)),*kvec=calloc(2*nf,sizeof(double)),*Rv=calloc(n,sizeof(double)),*q=calloc(n,sizeof(double)),*F=calloc(2*n,sizeof(double));
    int *incoming=calloc(nf,sizeof(int)),*inoffset=calloc(n+1,sizeof(int));
    if(!offset||!coefficient||!diagonal||!rhs[0]||!rhs[1]||!psi[0]||!psi[1]||!inverse||!kvec||!Rv||!q||!F||!incoming||!inoffset)return 2;
    int actual_threads=1;
#pragma omp parallel
    {
#pragma omp single
        actual_threads=omp_get_num_threads();
    }
    start=omp_get_wtime();
    for(int k=0;k<nf;k++){
        if(face[k].i<0||face[k].i>=n||face[k].j<0||face[k].j>=n||!(face[k].area>0))return 2;
        if(k&&face[k].i<face[k-1].i)return 2;
        offset[face[k].i+1]++;inoffset[face[k].j+1]++;
    }
    for(int i=0;i<n;i++){offset[i+1]+=offset[i];inoffset[i+1]+=inoffset[i];}
    int *cursor=malloc(n*sizeof(int));if(!cursor)return 2;for(int i=0;i<n;i++)cursor[i]=inoffset[i];
    for(int k=0;k<nf;k++)incoming[cursor[face[k].j]++]=k;free(cursor);
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){
        GfsVelocityFit fit;gfs_velocity_fit_init(&fit,2);double K[2][2]={{0}},V=cell[6*i+2];
        for(int k=offset[i];k<offset[i+1];k++){
            ProbeFace *f=face+k;double r[3]={f->rx,f->ry,0},c[3]={f->cx,f->cy,0},zero[3]={0};GfsVariationalFace b;
            if(!gfs_variational_face(2,r,c,f->area,&b)||!gfs_velocity_fit_add_face(&fit,r,zero,f->area))exit(2);
            coefficient[k]=f->area/(2*b.separation);diagonal[i]+=coefficient[k];
            for(int a=0;a<2;a++)for(int d=0;d<2;d++)K[a][d]+=b.bj[a]*r[d]/V;
        }
        rhs[0][i]=V*(K[1][1]-K[0][0])/sqrt(2.);rhs[1][i]=-V*(K[0][1]+K[1][0])/sqrt(2.);
        for(int a=0;a<2;a++)for(int d=0;d<2;d++)fit.R[a][d]=a==d;
        double G[3][3];if(!gfs_velocity_fit_solve(&fit,G,NULL))exit(2);
        for(int a=0;a<2;a++)for(int d=0;d<2;d++)inverse[4*i+2*a+d]=G[a][d];
    }
    double compatibility[2];int it[2]={0,0};double rel[2]={1,1};int rc=0;
    for(int b=0;b<2;b++){
        double sum=0,absolute=0;for(int i=0;i<n;i++){sum+=rhs[b][i];absolute+=fabs(rhs[b][i]);}compatibility[b]=sum;
        if(fabs(sum)>4096*DBL_EPSILON*fmax(1,absolute))return 2;
        for(int i=0;i<n;i++)rhs[b][i]-=sum/n;
    }
    for(int b=0;b<2;b++){rc=solve(rhs[b],psi[b],it+b,rel+b);if(rc)break;}
    double solve_seconds=omp_get_wtime()-start,identity=NAN;
    if(!rc){
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++)for(int k=offset[i];k<offset[i+1];k++){
            ProbeFace *f=face+k;int j=f->j;double d2=f->rx*f->rx+f->ry*f->ry,weight=f->area/(d2*sqrt(d2));
            double z0=inverse[4*i]*f->rx+inverse[4*i+1]*f->ry,z1=inverse[4*i+2]*f->rx+inverse[4*i+3]*f->ry;
            kvec[2*k]=weight*(psi[0][i]*z0+psi[1][i]*z1)/sqrt(2.);
            kvec[2*k+1]=weight*(psi[1][i]*z0-psi[0][i]*z1)/sqrt(2.);
            Rv[i]+=kvec[2*k]*(cell[6*j+4]-cell[6*i+4])+kvec[2*k+1]*(cell[6*j+5]-cell[6*i+5]);
            q[i]+=coefficient[k]*(cell[6*i+3]-cell[6*j+3]);
        }
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){
            for(int k=offset[i];k<offset[i+1];k++)for(int a=0;a<2;a++)F[2*i+a]-=q[i]*kvec[2*k+a];
            for(int s=inoffset[i];s<inoffset[i+1];s++){int k=incoming[s];for(int a=0;a<2;a++)F[2*i+a]+=q[face[k].i]*kvec[2*k+a];}
        }
        long double left=0,right=0;
        for(int i=0;i<n;i++){left+=(long double)cell[6*i+4]*F[2*i]+(long double)cell[6*i+5]*F[2*i+1];right+=(long double)q[i]*Rv[i];}
        identity=(double)(fabsl(left-right)/fmaxl(1,fabsl(left)+fabsl(right)));
        if(identity>1e-10)rc=2;
    }
    double elapsed=omp_get_wtime()-start;if(elapsed>budget&&!rc)rc=3;
    printf("{\"scope\":\"native CPU bounded coefficient/R/transpose screen only; weighted geometry and PDE cost not included\",\"cells\":%d,\"faces\":%d,\"actual_threads\":%d,\"maximum_CG_iterations\":256,\"relative_tolerance\":1e-10,\"extra_time_budget_seconds\":%.9g,\"candidate_seconds\":%.9g,\"preparation_and_coefficient_seconds\":%.9g,\"iterations\":[%d,%d],\"relative_residuals\":[%.17g,%.17g],\"rhs_compatibility\":[%.17g,%.17g],\"R_transpose_work_residual\":",n,nf,actual_threads,budget,elapsed,solve_seconds,it[0],it[1],rel[0],rel[1],compatibility[0],compatibility[1]);
    if(isfinite(identity))printf("%.17g",identity);else printf("null");
    printf(",\"status\":%d,\"accepted_steps\":0}\n",rc);return rc;
}

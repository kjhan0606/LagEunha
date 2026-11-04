#include "gfs_viscosity_implicit.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
struct GfsViscosityWorkspace {
    int capacity_n,capacity_e,n,dim,ne;
    int *offset,*adj,*cursor;
    unsigned char *active,*next;
    double *v0,*dv,*rhs,*r,*z,*p,*Ap,*trial,*heat,*inverse;
    const double *mass;
    const GfsViscosityEdge *edges;
    double h;
};
static double clock_seconds(void){
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock()/CLOCKS_PER_SEC;
#endif
}
GfsViscosityWorkspace *gfs_viscosity_workspace_create(void){return calloc(1,sizeof(GfsViscosityWorkspace));}
void gfs_viscosity_workspace_free(GfsViscosityWorkspace *w){
    if(!w)return;
    free(w->offset);free(w->adj);free(w->cursor);free(w->active);free(w->next);
    free(w->v0);free(w->dv);free(w->rhs);free(w->r);free(w->z);free(w->p);
    free(w->Ap);free(w->trial);free(w->heat);free(w->inverse);free(w);
}
static int resize_array(void **p,size_t count,size_t size){
    void *q=realloc(*p,count*size);if(!q&&count)return 0;*p=q;return 1;
}
static int reserve(GfsViscosityWorkspace *w,int n,int ne){
#define RESIZE(field,count,type) if(!resize_array((void**)&w->field,(count),sizeof(type)))return 0
    if(n>w->capacity_n){
        RESIZE(offset,n+1,int);RESIZE(cursor,n,int);
        RESIZE(v0,3*(size_t)n,double);RESIZE(dv,3*(size_t)n,double);
        RESIZE(rhs,3*(size_t)n,double);RESIZE(r,3*(size_t)n,double);
        RESIZE(z,3*(size_t)n,double);RESIZE(p,3*(size_t)n,double);
        RESIZE(Ap,3*(size_t)n,double);RESIZE(trial,3*(size_t)n,double);
        RESIZE(heat,n,double);RESIZE(inverse,9*(size_t)n,double);
        w->capacity_n=n;
    }
    if(ne>w->capacity_e){
        RESIZE(adj,2*(size_t)ne,int);RESIZE(active,ne,unsigned char);RESIZE(next,ne,unsigned char);
        w->capacity_e=ne;
    }
#undef RESIZE
    return 1;
}
static double edge_closing(const GfsViscosityWorkspace *w,int k,const double *v){
    const GfsViscosityEdge *e=w->edges+k;double q=0;
    for(int a=0;a<w->dim;a++)q+=e->normal[a]*(v[w->dim*e->i+a]-(e->j>=0?v[w->dim*e->j+a]:0));
    return q;
}
static void multiply(GfsViscosityWorkspace *w,const double *v,double *out){
    int d=w->dim;
#pragma omp parallel for schedule(static)
    for(int i=0;i<w->n;i++){
        double y[3]={0};
        for(int a=0;a<d;a++)y[a]=w->mass[i]*v[d*i+a];
        for(int l=w->offset[i];l<w->offset[i+1];l++){
            int k=w->adj[l];if(!w->active[k])continue;
            const GfsViscosityEdge *e=w->edges+k;int j=e->i==i?e->j:e->i;
            double g=0;
            for(int a=0;a<d;a++)g+=e->normal[a]*(v[d*i+a]-(j>=0?v[d*j+a]:0));
            for(int a=0;a<d;a++)y[a]+=w->h*e->coefficient*g*e->normal[a];
        }
        for(int a=0;a<d;a++)out[d*i+a]=y[a];
    }
}
static int preconditioner(GfsViscosityWorkspace *w){
    int bad=0,d=w->dim;
#pragma omp parallel for schedule(static) reduction(max:bad)
    for(int i=0;i<w->n;i++){
        double M[3][3]={{0}},L[3][3]={{0}};
        for(int a=0;a<d;a++)M[a][a]=w->mass[i];
        for(int l=w->offset[i];l<w->offset[i+1];l++){
            int k=w->adj[l];if(!w->active[k])continue;
            const GfsViscosityEdge *e=w->edges+k;
            for(int a=0;a<d;a++)for(int b=0;b<d;b++)M[a][b]+=w->h*e->coefficient*e->normal[a]*e->normal[b];
        }
        for(int a=0;a<d;a++)for(int b=0;b<=a;b++){
            double q=M[a][b];for(int c=0;c<b;c++)q-=L[a][c]*L[b][c];
            if(a==b){if(!(q>0)||!isfinite(q)){bad=1;q=1;}L[a][b]=sqrt(q);}
            else L[a][b]=q/L[b][b];
        }
        for(int b=0;b<d;b++){
            double y[3]={0},v[3]={0};
            for(int a=0;a<d;a++){
                double q=a==b?1:0;for(int c=0;c<a;c++)q-=L[a][c]*y[c];y[a]=q/L[a][a];
            }
            for(int a=d-1;a>=0;a--){
                double q=y[a];for(int c=a+1;c<d;c++)q-=L[c][a]*v[c];v[a]=q/L[a][a];
                w->inverse[9*i+3*a+b]=v[a];
            }
        }
    }
    return !bad;
}
static double scalar_product(const double *a,const double *b,int count){
    long double s=0;
#pragma omp parallel for schedule(static) reduction(+:s)
    for(int k=0;k<count;k++)s+=(long double)a[k]*b[k];
    return (double)s;
}
static double residual_and_precondition(GfsViscosityWorkspace *w,double *weighted_norm2){
    int d=w->dim;double rz=0,norm=0;
#pragma omp parallel for schedule(static) reduction(+:rz,norm)
    for(int i=0;i<w->n;i++){
        for(int a=0;a<d;a++){
            double z=0;for(int b=0;b<d;b++)z+=w->inverse[9*i+3*a+b]*w->r[d*i+b];
            w->z[d*i+a]=z;rz+=w->r[d*i+a]*z;norm+=w->r[d*i+a]*w->r[d*i+a]/w->mass[i];
        }
    }
    *weighted_norm2=norm;return rz;
}
static void cross(const double *x,const double *p,int dim,double out[3]){
    double a[3]={0},b[3]={0};for(int k=0;k<dim;k++){a[k]=x[k];b[k]=p[k];}
    out[0]=a[1]*b[2]-a[2]*b[1];out[1]=a[2]*b[0]-a[0]*b[2];out[2]=a[0]*b[1]-a[1]*b[0];
}
int gfs_viscosity_be(GfsViscosityWorkspace *w,int dim,int n,int ne,
        const GfsViscosityEdge *edges,const double *mass,const double *x,
        double *v,double *U,double h,GfsViscosityStats *stats){
    double start=clock_seconds();GfsViscosityStats s={0};
    if(!w||!stats||dim<1||dim>3||n<1||ne<0||!mass||!x||!v||!U||!(h>0)||!isfinite(h)||(!edges&&ne))return 1;
    if(!reserve(w,n,ne))return 1;
    w->n=n;w->dim=dim;w->ne=ne;w->mass=mass;w->edges=edges;w->h=h;
    int count=dim*n;memset(w->offset,0,(n+1)*sizeof(int));
    memcpy(w->v0,v,count*sizeof(double));memset(w->dv,0,count*sizeof(double));
    for(int i=0;i<n;i++){
        if(!(mass[i]>0)||!(U[i]>0)||!isfinite(mass[i])||!isfinite(U[i]))return 1;
        for(int a=0;a<dim;a++)if(!isfinite(v[dim*i+a])||!isfinite(x[dim*i+a]))return 1;
    }
    for(int k=0;k<ne;k++){
        const GfsViscosityEdge *e=edges+k;
        if(e->i<0||e->i>=n||e->j>=n||e->j < -1||e->j==e->i||!(e->coefficient>=0)||!isfinite(e->coefficient))return 1;
        double norm=0,dist2=0,projection=0;
        for(int a=0;a<dim;a++){
            norm+=e->normal[a]*e->normal[a];
            if(e->j>=0){double r=x[dim*e->j+a]-x[dim*e->i+a];dist2+=r*r;projection+=r*e->normal[a];}
        }
        if(fabs(norm-1)>256*DBL_EPSILON||!isfinite(norm))return 1;
        if(e->j>=0&&(!(projection>0)||fabs(dist2-projection*projection)>1024*DBL_EPSILON*dist2))return 1;
        w->offset[e->i+1]++;if(e->j>=0)w->offset[e->j+1]++;
        w->active[k]=e->coefficient>0&&edge_closing(w,k,w->v0)>0;
    }
    for(int i=0;i<n;i++)w->offset[i+1]+=w->offset[i];
    memcpy(w->cursor,w->offset,n*sizeof(int));
    for(int k=0;k<ne;k++){
        w->adj[w->cursor[edges[k].i]++]=k;
        if(edges[k].j>=0)w->adj[w->cursor[edges[k].j]++]=k;
    }
    int converged=0;
    for(int outer=0;outer<32;outer++){
        s.active_iterations=outer+1;
        if(!preconditioner(w))return 2;
        /* RHS=-h L v0, formed without subtracting orbital M v0. */
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){
            for(int a=0;a<dim;a++)w->rhs[dim*i+a]=0;
            for(int l=w->offset[i];l<w->offset[i+1];l++){
                int k=w->adj[l];if(!w->active[k])continue;
                const GfsViscosityEdge *e=edges+k;int j=e->i==i?e->j:e->i;double g=0;
                for(int a=0;a<dim;a++)g+=e->normal[a]*(w->v0[dim*i+a]-(j>=0?w->v0[dim*j+a]:0));
                for(int a=0;a<dim;a++)w->rhs[dim*i+a]-=h*e->coefficient*g*e->normal[a];
            }
        }
        double normb2=0;
#pragma omp parallel for schedule(static) reduction(+:normb2)
        for(int i=0;i<n;i++)for(int a=0;a<dim;a++)normb2+=w->rhs[dim*i+a]*w->rhs[dim*i+a]/mass[i];
        multiply(w,w->dv,w->Ap);
        for(int k=0;k<count;k++)w->r[k]=w->rhs[k]-w->Ap[k];
        double normr2,rz=residual_and_precondition(w,&normr2);
        memcpy(w->p,w->z,count*sizeof(double));
        int iteration=0;
        while(normr2>1e-28*normb2&&iteration<2000){
            multiply(w,w->p,w->Ap);double pap=scalar_product(w->p,w->Ap,count);
            if(!(pap>0)||!isfinite(pap)||!(rz>0))return 2;
            double alpha=rz/pap;
#pragma omp parallel for schedule(static)
            for(int k=0;k<count;k++){w->dv[k]+=alpha*w->p[k];w->r[k]-=alpha*w->Ap[k];}
            double next=residual_and_precondition(w,&normr2),beta=next/rz;
#pragma omp parallel for schedule(static)
            for(int k=0;k<count;k++)w->p[k]=w->z[k]+beta*w->p[k];
            rz=next;iteration++;
        }
        s.iterations+=iteration;
        if(iteration>=2000||!isfinite(normr2))return 2;
        for(int k=0;k<count;k++)w->trial[k]=w->v0[k]+w->dv[k];
        int changed=0;
        for(int k=0;k<ne;k++){
            w->next[k]=edges[k].coefficient>0&&edge_closing(w,k,w->trial)>0;
            if(w->next[k]!=w->active[k])changed++;
        }
        if(!changed){converged=1;break;}
        memcpy(w->active,w->next,ne);
    }
    if(!converged)return 3;
    /* Measure the equation after storing v0+dv, including addition roundoff. */
    for(int k=0;k<count;k++)w->dv[k]=w->trial[k]-w->v0[k];
    multiply(w,w->dv,w->Ap);
    double residual2=0,rhs2=0;
    for(int i=0;i<n;i++){
        w->heat[i]=0;
        for(int a=0;a<dim;a++){
            double dv=w->dv[dim*i+a],r=w->Ap[dim*i+a]-w->rhs[dim*i+a];
            residual2+=r*r/mass[i];rhs2+=w->rhs[dim*i+a]*w->rhs[dim*i+a]/mass[i];
            w->heat[i]+=.5*mass[i]*dv*dv;
        }
    }
    s.absolute_equation_residual=sqrt(residual2);
    s.relative_equation_residual=rhs2>0?sqrt(residual2/rhs2):sqrt(residual2);
    for(int k=0;k<ne;k++)if(w->active[k]){
        const GfsViscosityEdge *e=edges+k;double g=edge_closing(w,k,w->trial),Q=h*e->coefficient*g*g;
        s.active_edges++;w->heat[e->i]+=e->j>=0?.5*Q:Q;
        if(e->j>=0)w->heat[e->j]+=.5*Q;
        else{
            double impulse[3]={0},torque[3];
            for(int a=0;a<dim;a++){impulse[a]=-h*e->coefficient*g*e->normal[a];s.wall_impulse[a]+=impulse[a];}
            cross(x+dim*e->i,impulse,dim,torque);
            for(int a=0;a<3;a++)s.wall_torque[a]+=torque[a];
        }
    }
    long double energy_defect=0,heat=0,dp[3]={0},dJ[3]={0};s.min_heat=HUGE_VAL;
    for(int i=0;i<n;i++){
        if(!(w->heat[i]>=0)||!isfinite(w->heat[i])||!isfinite(U[i]+w->heat[i]))return 2;
        heat+=w->heat[i];s.min_heat=fmin(s.min_heat,w->heat[i]);s.max_heat_over_U=fmax(s.max_heat_over_U,w->heat[i]/U[i]);
        energy_defect+=w->heat[i];double momentum[3]={0},torque[3];
        for(int a=0;a<dim;a++){
            double dv=w->dv[dim*i+a];momentum[a]=mass[i]*dv;dp[a]+=momentum[a];
            energy_defect+=.5L*mass[i]*dv*(w->trial[dim*i+a]+w->v0[dim*i+a]);
        }
        cross(x+dim*i,momentum,dim,torque);for(int a=0;a<3;a++)dJ[a]+=torque[a];
    }
    s.heat=(double)heat;s.energy_defect=(double)energy_defect;
    for(int a=0;a<3;a++){s.momentum_defect[a]=(double)dp[a]-s.wall_impulse[a];s.angular_defect[a]=(double)dJ[a]-s.wall_torque[a];}
    for(int i=0;i<n;i++)U[i]+=w->heat[i];
    memcpy(v,w->trial,count*sizeof(double));s.solve_seconds=clock_seconds()-start;*stats=s;
    return 0;
}

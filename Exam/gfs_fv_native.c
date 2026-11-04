/* Shared 1D/2D/3D FV reconstruction, complete face flux and screened
 * weight-rate law. CPU/OpenMP backend; no GPU-to-host Krylov callbacks.
 * Geometry and SSPRK state updates are separate caller-owned operations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "gfs_fv_native.h"
#include "gfs_velocity_fit.h"
#ifdef GFS_NATIVE_KEP
#include "gfs_kep_ale.h"
typedef struct {
    double weight,wq,inv2d,density,compression,pressure,delta;
    GfsFvPrimitive left,right;
} FvFanPoint;
static void kep_entropy(int d,GfsFvPrimitive p,double gamma,double v[5]){
    double v2=0,beta=.5*p.rho/p.p;for(int c=0;c<5;c++)v[c]=0;
    for(int a=0;a<d;a++){v2+=p.u[a]*p.u[a];v[a+1]=2*beta*p.u[a];}
    v[0]=(gamma-log(p.p)+gamma*log(p.rho))/(gamma-1)-beta*v2;v[4]=-2*beta;
}
#else
typedef struct {GfsFvFan fan;double weight,wq,inv2d;} FvFanPoint;
#endif
typedef struct {int i,j,face,offset,count;double jacobian,flux[5],Vdot,wave,viscous,pressure_work,stress_heat,entropy_dissipation,entropy_boundary;} FvEdge;
struct GfsFvWorkspace {
    int n,ne,np,edge_capacity,point_capacity;GfsFvPrimitive *primitive;double (*gradient)[5][3],(*velocity)[3];
    FvEdge *edge;FvFanPoint *point;int *offset,*incident;
    double *z,*R,*diagonal,*scale,*direction,*r,*p,*Ap,*pre,*trial;
};
static double dot3(const double a[3],const double b[3],int d){double v=0;for(int k=0;k<d;k++)v+=a[k]*b[k];return v;}
static void *memory(size_t n){void *p=calloc(n?n:1,1);if(!p){fprintf(stderr,"FV native allocation failed\n");exit(12);}return p;}
GfsFvWorkspace *gfs_fv_workspace_create(void){return memory(sizeof(GfsFvWorkspace));}
void gfs_fv_workspace_free(GfsFvWorkspace *w){
    if(!w)return;free(w->primitive);free(w->gradient);free(w->velocity);free(w->edge);free(w->point);free(w->offset);free(w->incident);
    free(w->z);free(w->R);free(w->diagonal);free(w->scale);free(w->direction);free(w->r);free(w->p);free(w->Ap);free(w->pre);free(w->trial);free(w);
}
static void allocate(GfsFvWorkspace *w,int n){
    if(w->n==n)return;
    free(w->primitive);free(w->gradient);free(w->velocity);free(w->offset);
    w->primitive=memory(n*sizeof(*w->primitive));w->gradient=memory(n*sizeof(*w->gradient));w->velocity=memory(n*sizeof(*w->velocity));w->offset=memory((n+1)*sizeof(int));
#define VECTOR(field) free(w->field);w->field=memory(n*sizeof(double))
    VECTOR(z);VECTOR(R);VECTOR(diagonal);VECTOR(scale);VECTOR(direction);VECTOR(r);VECTOR(p);VECTOR(Ap);VECTOR(pre);VECTOR(trial);
#undef VECTOR
    w->n=n;
}
static double component(GfsFvPrimitive p,int c){return c==0?p.rho:c==4?p.p:p.u[c-1];}
static int mirrored(int flags,int a){return (flags&(1<<a))||(flags&(1<<(a+3)));}
static void neighbor(const GfsFvMesh *m,const GfsFvWorkspace *w,int i,const GfsFvNativeFace *f,
        GfsFvPrimitive *p,double r[3]){
    *p=w->primitive[f->j];
    for(int a=0;a<3;a++){
        double sign=mirrored(f->mirror,a)?-1:1;
        r[a]=f->r[a]+sign*m->geometry[f->j].centroid[a]-m->geometry[i].centroid[a];
        if(sign<0&&m->bc[a]==GFS_FV_REFLECT)p->u[a]*=-1;
    }
}
static int reconstruct(GfsFvWorkspace *w,GfsFvMesh *m){
    int n=m->n,d=m->dim,bad=0;
#pragma omp parallel for if(n>1024) schedule(static) reduction(max:bad)
    for(int i=0;i<n;i++)if(!gfs_fv_primitive(d,m->cell[i].Q,m->geometry[i].V,m->gamma,w->primitive+i))bad=1;
    if(bad)return 1;
#pragma omp parallel for if(n>1024) schedule(static) reduction(max:bad)
    for(int i=0;i<n;i++){
        GfsFvGeometry *g=m->geometry+i;GfsVelocityFit fit;gfs_velocity_fit_init(&fit,d);
        double rhs[5][3]={{0}},lo[5]={0},hi[5]={0},velocity_radius2=0,zero[3]={0};
        for(int k=0;k<g->faces;k++){
            GfsFvNativeFace *f=g->face+k;GfsFvPrimitive other;double r[3];neighbor(m,w,i,f,&other,r);
            if(!gfs_velocity_fit_add_face(&fit,r,zero,f->area)){bad=1;continue;}
            double distance2=dot3(r,r,d),angular=f->area;
            if(d==2)angular/=sqrt(distance2);if(d==3)angular/=distance2;
            double delta[5]={other.rho-w->primitive[i].rho,0,0,0,other.p-w->primitive[i].p};
            for(int a=0;a<d;a++)delta[1+a]=other.u[a]-w->primitive[i].u[a];
            double norm=0;for(int a=0;a<d;a++)norm+=delta[1+a]*delta[1+a];velocity_radius2=fmax(velocity_radius2,norm);
            for(int c=0;c<5;c++){
                lo[c]=fmin(lo[c],delta[c]);hi[c]=fmax(hi[c],delta[c]);
                for(int a=0;a<d;a++)rhs[c][a]+=angular*delta[c]*r[a]/distance2;
            }
        }
        /* Native Cholesky/conditioning checks, with R=I to obtain M^-1.
         * All five RHS rows then use the same single factorization. */
        for(int a=0;a<3;a++)for(int b=0;b<3;b++)fit.R[a][b]=a==b?1:0;
        double inverse[3][3],condition;
        if(!gfs_velocity_fit_solve(&fit,inverse,&condition)){bad=1;continue;}
        for(int c=0;c<5;c++)for(int a=0;a<3;a++){
            w->gradient[i][c][a]=0;
            for(int b=0;b<d;b++)if(a<d)w->gradient[i][c][a]+=rhs[c][b]*inverse[b][a];
        }
        double alpha[5]={1,1,1,1,1},velocity_alpha=1;
        for(int k=0;k<g->faces;k++){
            GfsFvNativeFace *f=g->face+k;GfsFvPrimitive other;double r[3];neighbor(m,w,i,f,&other,r);
            for(int c=0;c<5;c+=4){
                double prediction=.5*dot3(w->gradient[i][c],r,d);
                if(prediction>hi[c])alpha[c]=fmin(alpha[c],hi[c]/prediction);
                if(prediction<lo[c])alpha[c]=fmin(alpha[c],lo[c]/prediction);
            }
            double norm=0;
            for(int a=0;a<d;a++){double du=.5*dot3(w->gradient[i][1+a],r,d);norm+=du*du;}
            if(norm>velocity_radius2)velocity_alpha=fmin(velocity_alpha,sqrt(velocity_radius2/norm));
        }
        for(int c=0;c<5;c+=4){
            for(int q=0;q<g->points;q++){
                double r[3];for(int a=0;a<3;a++)r[a]=g->point[q].x[a]-g->centroid[a];
                double change=alpha[c]*dot3(w->gradient[i][c],r,d),base=component(w->primitive[i],c);
                if(change<0&&base+change<=.1*base)alpha[c]=fmin(alpha[c],.9*base/(-dot3(w->gradient[i][c],r,d)));
            }
            for(int a=0;a<d;a++)w->gradient[i][c][a]*=alpha[c];
        }
        for(int c=1;c<=d;c++)for(int a=0;a<d;a++)w->gradient[i][c][a]*=velocity_alpha;
        for(int a=0;a<3;a++){
            w->velocity[i][a]=a<d?w->primitive[i].u[a]-dot3(w->gradient[i][a+1],g->centroid,d):0;
            if(!isfinite(w->velocity[i][a]))bad=1;
        }
    }
    return bad;
}
static GfsFvPrimitive trace(const GfsFvWorkspace *w,const GfsFvMesh *m,int i,const double rel[3]){
    double r[3];for(int a=0;a<3;a++)r[a]=rel[a]-m->geometry[i].centroid[a];
    GfsFvPrimitive p=w->primitive[i];p.rho+=dot3(w->gradient[i][0],r,m->dim);p.p+=dot3(w->gradient[i][4],r,m->dim);
    for(int a=0;a<m->dim;a++)p.u[a]+=dot3(w->gradient[i][a+1],r,m->dim);return p;
}
static int build_points(GfsFvWorkspace *w,GfsFvMesh *m,GfsFvRhsStats *stats){
    int n=m->n,d=m->dim,ne=0,np=0,bad=0;
    for(int i=0;i<n;i++)for(int k=0;k<m->geometry[i].faces;k++){
        GfsFvNativeFace *f=m->geometry[i].face+k;
        if(f->j==i&&!f->mirror)return 1;
        if(f->mirror||i<f->j){ne++;np+=f->count;}
    }
    if(ne>w->edge_capacity){
        free(w->edge);free(w->incident);w->edge_capacity=2*ne;w->edge=memory(w->edge_capacity*sizeof(*w->edge));w->incident=memory(2*w->edge_capacity*sizeof(int));
    }
    if(np>w->point_capacity){free(w->point);w->point_capacity=2*np;w->point=memory(w->point_capacity*sizeof(*w->point));}
    w->ne=ne;w->np=np;
    memset(w->offset,0,(n+1)*sizeof(int));ne=np=0;
    for(int i=0;i<n;i++)for(int k=0;k<m->geometry[i].faces;k++){
        GfsFvNativeFace *f=m->geometry[i].face+k;
        if(f->mirror||i<f->j){
            FvEdge *e=w->edge+ne;e->i=i;e->j=f->mirror?-1:f->j;e->face=k;e->offset=np;e->count=f->count;
            w->offset[i+1]++;if(e->j>=0)w->offset[e->j+1]++;np+=f->count;ne++;
        }
    }
    for(int i=1;i<=n;i++)w->offset[i]+=w->offset[i-1];int *cursor=memory(n*sizeof(int));memcpy(cursor,w->offset,n*sizeof(int));
    for(int e=0;e<ne;e++){w->incident[cursor[w->edge[e].i]++]=e;if(w->edge[e].j>=0)w->incident[cursor[w->edge[e].j]++]=e;}free(cursor);
    int fallbacks=0;
#pragma omp parallel for if(n>1024) schedule(static) reduction(max:bad) reduction(+:fallbacks)
    for(int k=0;k<ne;k++){
        FvEdge *e=w->edge+k;int i=e->i,j=e->j;GfsFvNativeFace *f=m->geometry[i].face+e->face;
        if(j<0&&f->j!=i){fprintf(stderr,"FV reflected nonowner face dim=%d i=%d image_j=%d mirror=%d area=%.17g r=(%.17g %.17g %.17g)\n",d,i,f->j,f->mirror,f->area,f->r[0],f->r[1],f->r[2]);bad=1;continue;}
        for(int q=0;q<e->count;q++){
            FvFanPoint *point=w->point+e->offset+q;GfsFvPoint *g=m->geometry[i].point+f->point+q;
            GfsFvPrimitive left=trace(w,m,i,g->x),right;
            if(j>=0){double r[3];for(int a=0;a<3;a++)r[a]=g->x[a]-f->r[a];right=trace(w,m,j,r);}
            else{right=left;for(int a=0;a<d;a++)if(mirrored(f->mirror,a)&&m->bc[a]==GFS_FV_REFLECT)right.u[a]*=-1;}
#ifdef GFS_NATIVE_KEP
            /* Reconstruction determines physical qdot and the compression
             * stress; the EC core itself uses the saved cell means. */
            double reconstructed=0,delta=0;
            point->left=w->primitive[i];point->right=j>=0?w->primitive[j]:w->primitive[i];
            if(j<0)for(int a=0;a<d;a++)if(mirrored(f->mirror,a)&&m->bc[a]==GFS_FV_REFLECT)point->right.u[a]*=-1;
            for(int a=0;a<d;a++){
                reconstructed+=(right.u[a]-left.u[a])*f->normal[a];
                delta+=(point->right.u[a]-point->left.u[a])*f->normal[a];
            }
            double compression=delta<0&&reconstructed<0?fmin(-delta,-reconstructed):0;
            double hi=pow(m->geometry[i].V,1./d),hj=j>=0?pow(m->geometry[j].V,1./d):hi;
            compression*=.5*(hi+hj)/f->distance;
            double rhobar=.5*(point->left.rho+point->right.rho);
            double sound=.5*(sqrt(m->gamma*point->left.p/point->left.rho)+sqrt(m->gamma*point->right.p/point->right.rho));
            point->compression=rhobar*(sound*compression+2*compression*compression);point->delta=delta;
            point->density=gfs_kep_logmean(point->left.rho,point->right.rho);
            point->pressure=rhobar/(point->left.rho/(2*point->left.p)+point->right.rho/(2*point->right.p));
            if(!(point->density>0)||!isfinite(point->compression)){bad=1;continue;}
#else
            point->fan=gfs_fv_fan(d,left,right,f->normal,m->gamma);
            if(!point->fan.valid){fprintf(stderr,"FV invalid fan dim=%d i=%d j=%d rho=(%.17g %.17g) P=(%.17g %.17g)\n",d,i,j,left.rho,right.rho,left.p,right.p);bad=1;continue;}fallbacks+=point->fan.hll;
#endif
            point->weight=g->weight;point->inv2d=j>=0?1/(2*f->distance):0;
            if(j<0)point->wq=0;
            else{
                double mean[3],difference[3],offset[3];for(int a=0;a<3;a++){
                    mean[a]=.5*(w->velocity[i][a]+w->velocity[j][a]);difference[a]=w->velocity[j][a]-w->velocity[i][a];offset[a]=g->x[a]-.5*f->r[a];}
                point->wq=dot3(mean,f->normal,d)-dot3(difference,offset,d)/f->distance;
            }
        }
    }
    stats->hll_points=fallbacks;stats->points=np;return bad;
}
static int residual(GfsFvWorkspace *w,const GfsFvMesh *m,const double *z,double *R,int coefficients){
    int bad=0,n=m->n;
#pragma omp parallel for if(n>1024) schedule(static) reduction(max:bad)
    for(int e=0;e<w->ne;e++){
        FvEdge *edge=w->edge+e;double mass=0,jac=0;
        for(int q=0;q<edge->count;q++){
            FvFanPoint *p=w->point+edge->offset+q;double wn=p->wq+(edge->j>=0?p->inv2d*(z[edge->i]-z[edge->j]):0);
#ifdef GFS_NATIVE_KEP
            GfsFvNativeFace *face=m->geometry[edge->i].face+edge->face;
            double un=0;for(int a=0;a<m->dim;a++)un+=.5*(p->left.u[a]+p->right.u[a])*face->normal[a];
            mass+=p->weight*p->density*(un-wn);jac+=p->weight*p->inv2d*p->density;
#else
            GfsFvSample sample=gfs_fv_sample(&p->fan,wn);
            if(!sample.valid){bad=1;continue;}mass+=p->weight*sample.flux[0];jac+=p->weight*p->inv2d*sample.density;
#endif
        }
        edge->flux[0]=mass;if(coefficients)edge->jacobian=jac;
    }
#pragma omp parallel for if(n>1024) schedule(static)
    for(int i=0;i<n;i++){
        R[i]=m->cell[i].mu*z[i];if(coefficients)w->diagonal[i]=m->cell[i].mu;
        for(int k=w->offset[i];k<w->offset[i+1];k++){
            FvEdge *e=w->edge+w->incident[k];R[i]+=(e->i==i?-1:1)*e->flux[0];
            if(coefficients&&e->j>=0)w->diagonal[i]+=e->jacobian;
        }
    }
    return bad;
}
static double norm(const GfsFvWorkspace *w,const double *r){double value=0;for(int i=0;i<w->n;i++)value=fmax(value,fabs(r[i])/w->scale[i]);return value;}
static double product_sum(const double *a,const double *b,int n){
    double sum=0;
#pragma omp parallel for if(n>1024) schedule(static) reduction(+:sum)
    for(int i=0;i<n;i++)sum+=a[i]*b[i];return sum;
}
static void multiply(GfsFvWorkspace *w,const GfsFvMesh *m,const double *v,double *out){
#pragma omp parallel for if(m->n>1024) schedule(static)
    for(int i=0;i<m->n;i++){
        double value=m->cell[i].mu*v[i];
        for(int k=w->offset[i];k<w->offset[i+1];k++){
            FvEdge *e=w->edge+w->incident[k];if(e->j>=0){int j=e->i==i?e->j:e->i;value+=e->jacobian*(v[i]-v[j]);}
        }out[i]=value;
    }
}
static int pcg(GfsFvWorkspace *w,const GfsFvMesh *m){
    int n=m->n;for(int i=0;i<n;i++){w->direction[i]=0;w->r[i]=-w->R[i];w->pre[i]=w->r[i]/w->diagonal[i];w->p[i]=w->pre[i];}
    double rz=product_sum(w->r,w->pre,n),initial=rz;if(initial==0)return 0;
    for(int step=1;step<=64;step++){
        multiply(w,m,w->p,w->Ap);double den=product_sum(w->p,w->Ap,n);if(!(den>0))return -1;double a=rz/den;
        for(int i=0;i<n;i++){w->direction[i]+=a*w->p[i];w->r[i]-=a*w->Ap[i];w->pre[i]=w->r[i]/w->diagonal[i];}
        double next=product_sum(w->r,w->pre,n);if(next<=1e-24*initial)return step;
        for(int i=0;i<n;i++)w->p[i]=w->pre[i]+(next/rz)*w->p[i];rz=next;
    }
    return -1;
}
#ifndef GFS_NATIVE_KEP
static double merit(const GfsFvMesh *m,const double *R){double value=0;for(int i=0;i<m->n;i++)value+=R[i]*R[i]/m->cell[i].mu;return value;}
#endif
static int solve(GfsFvWorkspace *w,const GfsFvMesh *m,GfsFvRhsStats *s){
    memset(w->z,0,m->n*sizeof(double));
    for(int i=0;i<m->n;i++){
        double h0=sqrt(m->cell[i].m0/m->cell[i].mu),speed=sqrt(dot3(w->primitive[i].u,w->primitive[i].u,m->dim))+sqrt(m->gamma*w->primitive[i].p/w->primitive[i].rho);
        w->scale[i]=m->cell[i].mu*h0*speed;if(!(w->scale[i]>0))return 1;
    }
#ifdef GFS_NATIVE_KEP
    if(residual(w,m,w->z,w->R,1))return 1;s->root_residual=norm(w,w->R);
    if(s->root_residual<=1e-12)return 0;
    int count=pcg(w,m);if(count<0)return 3;s->cg=count;
    memcpy(w->z,w->direction,m->n*sizeof(double));
    if(residual(w,m,w->z,w->R,0))return 1;s->root_residual=norm(w,w->R);
    return s->root_residual<=1e-12?0:2;
#else
    for(int iteration=0;iteration<=8;iteration++){
        if(residual(w,m,w->z,w->R,1))return 1;s->root_residual=norm(w,w->R);s->newton=iteration;
        if(s->root_residual<=1e-12)return 0;if(iteration==8)return 2;
        int count=pcg(w,m);if(count<0)return 3;s->cg+=count;
        double old=merit(m,w->R),factor=1;int accepted=0;
        for(int retry=0;retry<16;retry++){
            for(int i=0;i<m->n;i++)w->trial[i]=w->z[i]+factor*w->direction[i];
            if(residual(w,m,w->trial,w->pre,0))return 1;
            if(merit(m,w->pre)<old*(1-1e-4*factor)){memcpy(w->z,w->trial,m->n*sizeof(double));accepted=1;break;}
            factor*=.5;s->backtracks++;
        }
        if(!accepted)return 4;
    }return 2;
#endif
}
int gfs_fv_native_rhs(GfsFvWorkspace *w,GfsFvMesh *m,GfsFvRate *rate,GfsFvRhsStats *s){
    memset(s,0,sizeof(*s));allocate(w,m->n);double mark=omp_get_wtime();
    if(reconstruct(w,m))return 1;s->reconstruction_seconds=omp_get_wtime()-mark;mark=omp_get_wtime();
    if(build_points(w,m,s))return 2;s->fan_seconds=omp_get_wtime()-mark;mark=omp_get_wtime();
    int code=solve(w,m,s);s->solve_seconds=omp_get_wtime()-mark;if(code){fprintf(stderr,"FV screened solve failed code=%d Newton=%d CG=%d residual=%.17g\n",code,s->newton,s->cg,s->root_residual);return 10+code;}
    mark=omp_get_wtime();int bad=0;
#pragma omp parallel for if(m->n>1024) schedule(static) reduction(max:bad)
    for(int e=0;e<w->ne;e++){
        FvEdge *f=w->edge+e;memset(f->flux,0,sizeof(f->flux));f->Vdot=f->wave=f->viscous=f->pressure_work=f->stress_heat=f->entropy_dissipation=f->entropy_boundary=0;
        for(int q=0;q<f->count;q++){
            FvFanPoint *p=w->point+f->offset+q;double wn=p->wq+(f->j>=0?p->inv2d*(w->z[f->i]-w->z[f->j]):0);
#ifdef GFS_NATIVE_KEP
            GfsFvNativeFace *face=m->geometry[f->i].face+f->face;
            double flux[5],rho;
            if(!gfs_kep_flux(m->dim,p->left,p->right,face->normal,m->gamma,wn,p->compression,flux,&rho)){bad=1;continue;}
            for(int c=0;c<5;c++)f->flux[c]+=p->weight*flux[c];
            double left=dot3(p->left.u,face->normal,m->dim),right=dot3(p->right.u,face->normal,m->dim);
            double signal=fmax(fabs(left-wn)+sqrt(m->gamma*p->left.p/p->left.rho),fabs(right-wn)+sqrt(m->gamma*p->right.p/p->right.rho));
            f->Vdot+=p->weight*wn;f->wave+=p->weight*signal;
            if(p->compression>0)f->viscous+=p->weight*p->compression/(-p->delta);
            double fraction=f->j>=0?1:.5;
            f->pressure_work+=fraction*p->weight*(p->pressure+p->compression)*p->delta;
            f->stress_heat-=fraction*p->weight*p->compression*p->delta;
            f->entropy_dissipation+=fraction*p->weight*.5*(p->left.rho/p->left.p+p->right.rho/p->right.p)*p->compression*p->delta;
            if(f->j<0){
                double vl[5],vr[5],entropy_flux=0;kep_entropy(m->dim,p->left,m->gamma,vl);kep_entropy(m->dim,p->right,m->gamma,vr);
                for(int c=0;c<5;c++)entropy_flux+=.5*(vl[c]+vr[c])*flux[c];
                entropy_flux-=.5*(p->left.rho*left+p->right.rho*right);entropy_flux+=wn*.5*(p->left.rho+p->right.rho);
                f->entropy_boundary+=p->weight*entropy_flux;
            }
#else
            GfsFvSample value=gfs_fv_sample(&p->fan,wn);if(!value.valid){bad=1;continue;}
            for(int c=0;c<5;c++)f->flux[c]+=p->weight*value.flux[c];
            f->Vdot+=p->weight*wn;f->wave+=p->weight*fmax(fabs(p->fan.speed[0]-wn),fabs(p->fan.speed[2]-wn));
#endif
        }
    }
    double wave=HUGE_VAL,limit=HUGE_VAL,link=0,viscous_dt=HUGE_VAL,kinetic_rate=0,entropy_rate=0,entropy_absolute_rate=0;
#pragma omp parallel for if(m->n>1024) schedule(static) reduction(max:bad,link) reduction(min:wave,limit,viscous_dt) reduction(+:kinetic_rate,entropy_rate,entropy_absolute_rate)
    for(int i=0;i<m->n;i++){
        GfsFvRate *r=rate+i;memset(r,0,sizeof(*r));memcpy(r->qdot,w->velocity[i],sizeof(r->qdot));r->z=w->z[i];r->residual=w->R[i];double signals=0,damping=0;
        for(int k=w->offset[i];k<w->offset[i+1];k++){
            FvEdge *f=w->edge+w->incident[k];double sign=f->i==i?1:-1;
            for(int c=0;c<5;c++)r->dQ[c]-=sign*f->flux[c];r->Vdot+=sign*f->Vdot;signals+=f->wave;damping+=f->viscous;
        }
        if(!(signals>0)){bad=1;continue;}
        wave=fmin(wave,.4*2*m->geometry[i].V/signals);
        if(damping>0)viscous_dt=fmin(viscous_dt,.4*m->cell[i].Q[0]/damping);
        double u2=0;for(int a=0;a<m->dim;a++){double u=w->primitive[i].u[a];kinetic_rate+=u*r->dQ[a+1];u2+=u*u;}kinetic_rate-=.5*u2*r->dQ[0];
#ifdef GFS_NATIVE_KEP
        double entropy_vector[5],entropy_cell=-w->primitive[i].rho*r->Vdot;kep_entropy(m->dim,w->primitive[i],m->gamma,entropy_vector);
        for(int c=0;c<5;c++)entropy_cell+=entropy_vector[c]*r->dQ[c];entropy_rate+=entropy_cell;entropy_absolute_rate+=fabs(entropy_cell);
#endif
        double ray=gfs_fv_admissible_interval(m->dim,m->cell[i].Q,r->dQ);if(!(ray>0))bad=1;limit=fmin(limit,.9*ray);
        link=fmax(link,fabs(r->dQ[0]+m->cell[i].mu*r->z+r->residual)/w->scale[i]);
    }
    for(int e=0;e<w->ne;e++){
        if(w->edge[e].j<0)for(int c=0;c<5;c++)s->boundary_flux[c]+=w->edge[e].flux[c];
        s->kinetic_pressure_work+=w->edge[e].pressure_work;s->stress_thermal_rate+=w->edge[e].stress_heat;
        s->entropy_boundary_flux+=w->edge[e].entropy_boundary;s->entropy_dissipation+=w->edge[e].entropy_dissipation;
    }
    s->entropy_rate=entropy_rate;s->entropy_balance_defect=entropy_rate+s->entropy_boundary_flux-s->entropy_dissipation;
#ifdef GFS_NATIVE_KEP
    if(fabs(s->entropy_balance_defect)>1e-10*(1+entropy_absolute_rate+fabs(s->entropy_boundary_flux)+fabs(s->entropy_dissipation))){fprintf(stderr,"KEP native entropy ledger failed defect=%.17g\n",s->entropy_balance_defect);bad=1;}
#endif
    s->kinetic_rate=kinetic_rate;s->viscous_dt=viscous_dt;
    s->wave_dt=fmin(wave,viscous_dt);s->admissible_dt=limit;s->mass_link_defect=link;s->flux_seconds=omp_get_wtime()-mark;return bad?3:0;
}
int gfs_fv_native_velocity(GfsFvWorkspace *w,GfsFvMesh *m,double (*velocity)[3]){
    allocate(w,m->n);if(reconstruct(w,m))return 1;memcpy(velocity,w->velocity,m->n*sizeof(*velocity));return 0;
}

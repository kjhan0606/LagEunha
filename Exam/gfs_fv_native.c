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
#if defined(GFS_NATIVE_GLOBAL_VELOCITY) && (!defined(GFS_NATIVE_KEP) || !defined(GFS_NATIVE_KEP_FILTER))
#error Global physical velocity is registered with the KEP/EC filtered-stress candidate
#endif
#if defined(GFS_NATIVE_BLOCK_JACOBI) && !defined(GFS_NATIVE_GLOBAL_VELOCITY)
#error The conditional block preconditioner is registered for global physical velocity
#endif
#if defined(GFS_NATIVE_IC0) && (defined(GFS_NATIVE_BLOCK_JACOBI) || !defined(GFS_NATIVE_GLOBAL_VELOCITY))
#error IC0 is a separate preconditioner for the common KEP physical field
#endif
#if defined(GFS_NATIVE_WALL_VELOCITY) && !defined(GFS_NATIVE_GLOBAL_VELOCITY)
#error The constrained wall field is registered with global physical velocity
#endif
#if (defined(GFS_NATIVE_SCALAR_VELOCITY) || defined(GFS_NATIVE_CENTROID_STRAIN)) && !defined(GFS_NATIVE_WALL_VELOCITY)
#error The common limiter/strain candidates retain the constrained physical wall field
#endif
#ifdef GFS_NATIVE_KEP
#include "gfs_kep_ale.h"
#ifdef GFS_NATIVE_KEP_FILTER
#include "gfs_kep_stress.h"
#endif
typedef struct {
    double weight,wq,inv2d,density,compression,pressure,delta;
    GfsFvPrimitive left,right;
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
    double original_compression,fraction;
#endif
} FvFanPoint;
static void kep_entropy(int d,GfsFvPrimitive p,double gamma,double v[5]){
    double v2=0,beta=.5*p.rho/p.p;for(int c=0;c<5;c++)v[c]=0;
    for(int a=0;a<d;a++){v2+=p.u[a]*p.u[a];v[a+1]=2*beta*p.u[a];}
    v[0]=(gamma-log(p.p)+gamma*log(p.rho))/(gamma-1)-beta*v2;v[4]=-2*beta;
}
#else
typedef struct {GfsFvFan fan;double weight,wq,inv2d;} FvFanPoint;
#endif
typedef struct {int i,j,face,offset,count;double jacobian,flux[5],Vdot,wave,viscous,pressure_work,kinetic_boundary,stress_heat,entropy_dissipation,entropy_boundary;} FvEdge;
struct GfsFvWorkspace {
    int n,ne,np,edge_capacity,point_capacity;GfsFvPrimitive *primitive;double (*gradient)[5][3],(*velocity)[3];
    int pcg_iterations,pcg_failure,pcg_stop_reason,pcg_true_replacements;
    double pcg_last_relative,pcg_first_candidate_true;
#ifdef GFS_PCG_TRUE_HISTORY
    double pcg_true_history[65],pcg_recursive_history[65],pcg_relative_history[65];
#endif
#ifdef GFS_NATIVE_WALL_VELOCITY
    unsigned char *wall_mask;double (*wall_relative)[3],(*wall_distance)[3];
#endif
#ifdef GFS_NATIVE_BLOCK_JACOBI
    int *strongest,*partner;double *block_a,*block_b;
#endif
#ifdef GFS_NATIVE_IC0
    int *ic_offset,*ic_column,ic_capacity,ic_nnz,ic_applications,ic_failure_row;
    double *ic_matrix,*ic_lower,*ic_pivot;
    double ic_assembly_seconds,ic_factor_seconds,ic_apply_seconds,ic_min_pivot;
#endif
#ifdef GFS_NATIVE_GLOBAL_VELOCITY
    int *velocity_owner;
    int owner_changed,owner_switches,query_max_hops;
    long long query_hops;
    double query_seconds;
#endif
    FvEdge *edge;FvFanPoint *point;int *offset,*incident;
    double *z,*R,*diagonal,*scale,*direction,*r,*p,*Ap,*pre,*trial;
#ifdef GFS_NATIVE_KEP_FILTER
    double *stress_fraction;
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
    double *parallel_energy,*perpendicular_energy,*data_energy,*bulk_fraction,*cross_energy;
    double (*raw_velocity_gradient)[3][3];
#endif
#endif
};
static double dot3(const double a[3],const double b[3],int d){double v=0;for(int k=0;k<d;k++)v+=a[k]*b[k];return v;}
static void *memory(size_t n){void *p=calloc(n?n:1,1);if(!p){fprintf(stderr,"FV native allocation failed\n");exit(12);}return p;}
GfsFvWorkspace *gfs_fv_workspace_create(void){return memory(sizeof(GfsFvWorkspace));}
void gfs_fv_workspace_free(GfsFvWorkspace *w){
    if(!w)return;free(w->primitive);free(w->gradient);free(w->velocity);free(w->edge);free(w->point);free(w->offset);free(w->incident);
    free(w->z);free(w->R);free(w->diagonal);free(w->scale);free(w->direction);free(w->r);free(w->p);free(w->Ap);free(w->pre);free(w->trial);
#ifdef GFS_NATIVE_WALL_VELOCITY
    free(w->wall_mask);free(w->wall_relative);free(w->wall_distance);
#endif
#ifdef GFS_NATIVE_BLOCK_JACOBI
    free(w->strongest);free(w->partner);free(w->block_a);free(w->block_b);
#endif
#ifdef GFS_NATIVE_IC0
    free(w->ic_offset);free(w->ic_column);free(w->ic_matrix);free(w->ic_lower);free(w->ic_pivot);
#endif
#ifdef GFS_NATIVE_GLOBAL_VELOCITY
    free(w->velocity_owner);
#endif
#ifdef GFS_NATIVE_KEP_FILTER
    free(w->stress_fraction);
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
    free(w->parallel_energy);free(w->perpendicular_energy);free(w->data_energy);free(w->bulk_fraction);free(w->cross_energy);free(w->raw_velocity_gradient);
#endif
#endif
    free(w);
}
static void allocate(GfsFvWorkspace *w,int n){
    if(w->n==n)return;
    free(w->primitive);free(w->gradient);free(w->velocity);free(w->offset);
    w->primitive=memory(n*sizeof(*w->primitive));w->gradient=memory(n*sizeof(*w->gradient));w->velocity=memory(n*sizeof(*w->velocity));w->offset=memory((n+1)*sizeof(int));
#ifdef GFS_NATIVE_IC0
    free(w->ic_offset);free(w->ic_pivot);w->ic_offset=memory((n+1)*sizeof(int));w->ic_pivot=memory(n*sizeof(double));
#endif
#ifdef GFS_NATIVE_WALL_VELOCITY
    free(w->wall_mask);free(w->wall_relative);free(w->wall_distance);
    w->wall_mask=memory(n*sizeof(*w->wall_mask));w->wall_relative=memory(n*sizeof(*w->wall_relative));w->wall_distance=memory(n*sizeof(*w->wall_distance));
#endif
#ifdef GFS_NATIVE_GLOBAL_VELOCITY
    free(w->velocity_owner);w->velocity_owner=memory(n*sizeof(int));for(int i=0;i<n;i++)w->velocity_owner[i]=-1;
#endif
#define VECTOR(field) free(w->field);w->field=memory(n*sizeof(double))
    VECTOR(z);VECTOR(R);VECTOR(diagonal);VECTOR(scale);VECTOR(direction);VECTOR(r);VECTOR(p);VECTOR(Ap);VECTOR(pre);VECTOR(trial);
#ifdef GFS_NATIVE_BLOCK_JACOBI
    free(w->strongest);free(w->partner);w->strongest=memory(n*sizeof(int));w->partner=memory(n*sizeof(int));VECTOR(block_a);VECTOR(block_b);
#endif
#ifdef GFS_NATIVE_KEP_FILTER
    VECTOR(stress_fraction);
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
    VECTOR(parallel_energy);VECTOR(perpendicular_energy);VECTOR(data_energy);VECTOR(bulk_fraction);VECTOR(cross_energy);
    free(w->raw_velocity_gradient);w->raw_velocity_gradient=memory(n*sizeof(*w->raw_velocity_gradient));
#endif
#endif
#undef VECTOR
    w->n=n;
}
static double component(GfsFvPrimitive p,int c){return c==0?p.rho:c==4?p.p:p.u[c-1];}
static int mirrored(int flags,int a){return (flags&(1<<a))||(flags&(1<<(a+3)));}
static int face_bc(const GfsFvMesh *m,const GfsFvNativeFace *f,int a){return gfs_fv_boundary_type(m->bc[a],!!(f->mirror&(1<<(a+3))));}
#ifdef GFS_NATIVE_WALL_VELOCITY
static int constrain_wall(GfsFvWorkspace *w,const GfsFvMesh *m,int i){
    const GfsFvGeometry *g=m->geometry+i;unsigned char mask=0;
    for(int k=0;k<g->faces;k++){
        const GfsFvNativeFace *f=g->face+k;
        for(int a=0;a<m->dim;a++)if(mirrored(f->mirror,a)&&face_bc(m,f,a)==GFS_FV_REFLECT){
            if(f->j!=i)return 1;
            for(int b=0;b<m->dim;b++)if(b!=a&&fabs(f->normal[b])>4096*DBL_EPSILON)return 1;
            mask|=f->mirror&((1<<a)|(1<<(a+3)));
        }
    }
    w->wall_mask[i]=mask;
    for(int a=0;a<m->dim;a++)if(mirrored(mask,a)){
        int both=(mask&(1<<a))&&(mask&(1<<(a+3)));
        if(both&&w->primitive[i].u[a]!=0){
            fprintf(stderr,"Native affine wall field infeasible cell=%d axis=%d normal_mean=%.17g: opposite walls require zero normal mean\n",i,a,w->primitive[i].u[a]);return 1;
        }
        double wall=m->lo[a]+((mask&(1<<(a+3)))?m->L[a]:0);
        w->wall_relative[i][a]=wall-m->cell[i].q[a];
        double distance=g->centroid[a]-w->wall_relative[i][a];
        if(!isfinite(distance)||distance==0)return 1;
        if(!both&&((mask&(1<<a))?!(distance>0):!(distance<0)))return 1;
        w->wall_distance[i][a]=distance;
        for(int b=0;b<m->dim;b++)w->gradient[i][a+1][b]=0;
        if(!both)w->gradient[i][a+1][a]=w->primitive[i].u[a]/distance;
    }
    return 0;
}
#endif
/* physical is available for generator queries, so x-wall is evaluated
 * directly rather than by subtracting two large affine contributions. */
static double physical_velocity(const GfsFvWorkspace *w,const GfsFvMesh *m,int i,int a,
        const double relative[3],const double physical[3]){
#ifdef GFS_NATIVE_WALL_VELOCITY
    if(mirrored(w->wall_mask[i],a)){
        if((w->wall_mask[i]&(1<<a))&&(w->wall_mask[i]&(1<<(a+3))))return 0;
        double wall=m->lo[a]+((w->wall_mask[i]&(1<<(a+3)))?m->L[a]:0);
        double distance=physical?physical[a]-wall:relative[a]-w->wall_relative[i][a];
        return w->primitive[i].u[a]*(distance/w->wall_distance[i][a]);
    }
#else
    (void)physical;
#endif
    double offset[3];for(int b=0;b<3;b++)offset[b]=relative[b]-m->geometry[i].centroid[b];
    return w->primitive[i].u[a]+dot3(w->gradient[i][a+1],offset,m->dim);
}
static void neighbor(const GfsFvMesh *m,const GfsFvWorkspace *w,int i,const GfsFvNativeFace *f,
        GfsFvPrimitive *p,double r[3]){
    *p=w->primitive[f->j];
    for(int a=0;a<3;a++){
        double sign=mirrored(f->mirror,a)?-1:1;
        r[a]=f->r[a]+sign*m->geometry[f->j].centroid[a]-m->geometry[i].centroid[a];
        if(sign<0&&face_bc(m,f,a)==GFS_FV_REFLECT)p->u[a]*=-1;
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
#ifdef GFS_NATIVE_SCALAR_VELOCITY
        (void)velocity_radius2;
#endif
        for(int k=0;k<g->faces;k++){
            GfsFvNativeFace *f=g->face+k;GfsFvPrimitive other;double r[3];neighbor(m,w,i,f,&other,r);
            if(!gfs_velocity_fit_add_face(&fit,r,zero,f->area)){bad=1;continue;}
            double distance2=dot3(r,r,d),angular=f->area;
            if(d==2)angular/=sqrt(distance2);if(d==3)angular/=distance2;
            double delta[5]={other.rho-w->primitive[i].rho,0,0,0,other.p-w->primitive[i].p};
            for(int a=0;a<d;a++)delta[1+a]=other.u[a]-w->primitive[i].u[a];
#ifndef GFS_NATIVE_SCALAR_VELOCITY
            double norm=0;for(int a=0;a<d;a++)norm+=delta[1+a]*delta[1+a];velocity_radius2=fmax(velocity_radius2,norm);
#endif
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
#ifdef GFS_NATIVE_WALL_VELOCITY
        if(constrain_wall(w,m,i)){bad=1;continue;}
#endif
        double alpha[5]={1,1,1,1,1},velocity_alpha=1;
#ifdef GFS_NATIVE_SCALAR_VELOCITY
        (void)velocity_alpha;
#endif
#ifdef GFS_NATIVE_KEP_FILTER
        double Eparallel=0,Eperpendicular=0,G[3][3]={{0}};
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
        double Edata=0,Ecross=0;
#endif
        for(int a=0;a<d;a++)for(int b=0;b<d;b++)G[a][b]=w->gradient[i][a+1][b];
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
        memcpy(w->raw_velocity_gradient[i],G,sizeof(G));
#endif
#endif
        for(int k=0;k<g->faces;k++){
            GfsFvNativeFace *f=g->face+k;GfsFvPrimitive other;double r[3];neighbor(m,w,i,f,&other,r);
#ifdef GFS_NATIVE_KEP_FILTER
            double r2=dot3(r,r,d),omega=f->area/r2;
            if(d==2)omega/=sqrt(r2);if(d==3)omega/=r2;
            for(int a=0;a<d;a++){
                double predicted=dot3(G[a],r,d),delta=other.u[a]-w->primitive[i].u[a];
                Eparallel+=omega*predicted*predicted;Eperpendicular+=omega*(delta-predicted)*(delta-predicted);
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
                Edata+=omega*delta*delta;Ecross+=2*omega*(delta-predicted)*predicted;
#endif
            }
#endif
#ifdef GFS_NATIVE_SCALAR_VELOCITY
            for(int c=0;c<5;c++){
                if(c>0&&c<4&&(c>d||mirrored(w->wall_mask[i],c-1)))continue;
#else
            for(int c=0;c<5;c+=4){
#endif
                double prediction=.5*dot3(w->gradient[i][c],r,d);
                if(prediction>hi[c])alpha[c]=fmin(alpha[c],hi[c]/prediction);
                if(prediction<lo[c])alpha[c]=fmin(alpha[c],lo[c]/prediction);
            }
#ifndef GFS_NATIVE_SCALAR_VELOCITY
            double norm=0;
            for(int a=0;a<d;a++){double du=.5*dot3(w->gradient[i][1+a],r,d);norm+=du*du;}
            if(norm>velocity_radius2)velocity_alpha=fmin(velocity_alpha,sqrt(velocity_radius2/norm));
#endif
        }
#ifdef GFS_NATIVE_KEP_FILTER
        double bulk;w->stress_fraction[i]=gfs_kep_stress_fraction(d,G,Eparallel,Eperpendicular,&bulk);
        if(!(w->stress_fraction[i]>=0)||w->stress_fraction[i]>1||!isfinite(w->stress_fraction[i]))bad=1;
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
        w->parallel_energy[i]=Eparallel;w->perpendicular_energy[i]=Eperpendicular;w->data_energy[i]=Edata;w->bulk_fraction[i]=bulk;w->cross_energy[i]=Ecross;
#endif
#endif
        for(int c=0;c<5;c+=4){
            for(int q=0;q<g->points;q++){
                double r[3];for(int a=0;a<3;a++)r[a]=g->point[q].x[a]-g->centroid[a];
                double change=alpha[c]*dot3(w->gradient[i][c],r,d),base=component(w->primitive[i],c);
                if(change<0&&base+change<=.1*base)alpha[c]=fmin(alpha[c],.9*base/(-dot3(w->gradient[i][c],r,d)));
            }
            for(int a=0;a<d;a++)w->gradient[i][c][a]*=alpha[c];
        }
        for(int c=1;c<=d;c++){
#ifdef GFS_NATIVE_WALL_VELOCITY
            if(mirrored(w->wall_mask[i],c-1))continue;
#endif
#ifdef GFS_NATIVE_SCALAR_VELOCITY
            for(int a=0;a<d;a++)w->gradient[i][c][a]*=alpha[c];
#else
            for(int a=0;a<d;a++)w->gradient[i][c][a]*=velocity_alpha;
#endif
        }
        for(int a=0;a<3;a++){
#ifdef GFS_NATIVE_WALL_VELOCITY
            w->velocity[i][a]=a<d?physical_velocity(w,m,i,a,zero,m->cell[i].q):0;
#else
            w->velocity[i][a]=a<d?w->primitive[i].u[a]-dot3(w->gradient[i][a+1],g->centroid,d):0;
#endif
            if(!isfinite(w->velocity[i][a]))bad=1;
        }
    }
#ifdef GFS_NATIVE_GLOBAL_VELOCITY
    if(bad)return bad;
    double start=omp_get_wtime();int changed=0,switches=0,max_hops=0;long long hops_total=0;
#pragma omp parallel for if(n>1024) schedule(static) reduction(max:bad,max_hops) reduction(+:changed,switches,hops_total)
    for(int i=0;i<n;i++){
        int owner,hops;double relative[3];int code=gfs_fv_point_owner(m,i,m->cell[i].q,&owner,relative,&hops);
        if(code){fprintf(stderr,"Native physical velocity owner query failed cell=%d code=%d\n",i,code);bad=1;continue;}
        hops_total+=hops;if(hops>max_hops)max_hops=hops;
        changed+=owner!=i;switches+=w->velocity_owner[i]>=0&&w->velocity_owner[i]!=owner;w->velocity_owner[i]=owner;
        if(owner!=i){
            for(int a=0;a<d;a++)w->velocity[i][a]=physical_velocity(w,m,owner,a,relative,m->cell[i].q);
        }
    }
    w->query_seconds=omp_get_wtime()-start;w->owner_changed=changed;w->owner_switches=switches;w->query_hops=hops_total;w->query_max_hops=max_hops;
#endif
    return bad;
}
static GfsFvPrimitive trace(const GfsFvWorkspace *w,const GfsFvMesh *m,int i,const double rel[3]){
    double r[3];for(int a=0;a<3;a++)r[a]=rel[a]-m->geometry[i].centroid[a];
    GfsFvPrimitive p=w->primitive[i];p.rho+=dot3(w->gradient[i][0],r,m->dim);p.p+=dot3(w->gradient[i][4],r,m->dim);
    for(int a=0;a<m->dim;a++){
#ifdef GFS_NATIVE_WALL_VELOCITY
        p.u[a]=physical_velocity(w,m,i,a,rel,NULL);
#else
        /* Preserve the historical operation order outside the opt-in field. */
        p.u[a]+=dot3(w->gradient[i][a+1],r,m->dim);
#endif
    }return p;
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
#ifdef GFS_NATIVE_CENTROID_STRAIN
        double physical_separation[3];
        for(int a=0;a<3;a++)physical_separation[a]=f->r[a]+(mirrored(f->mirror,a)?-1:1)*m->geometry[f->j].centroid[a]-m->geometry[i].centroid[a];
        double strain_distance=dot3(physical_separation,f->normal,d);
        if(!(strain_distance>0)||!isfinite(strain_distance)){fprintf(stderr,"FV nonpositive centroid-normal strain distance i=%d j=%d value=%.17g\n",i,j,strain_distance);bad=1;continue;}
        double hi=pow(m->geometry[i].V,1./d),hj=j>=0?pow(m->geometry[j].V,1./d):hi;
        double strain_scale=.5*(hi+hj)/strain_distance;
#endif
        for(int q=0;q<e->count;q++){
            FvFanPoint *point=w->point+e->offset+q;GfsFvPoint *g=m->geometry[i].point+f->point+q;
            GfsFvPrimitive left=trace(w,m,i,g->x),right;
            if(j>=0){double r[3];for(int a=0;a<3;a++)r[a]=g->x[a]-f->r[a];right=trace(w,m,j,r);}
            else{right=left;for(int a=0;a<d;a++)if(mirrored(f->mirror,a)&&face_bc(m,f,a)==GFS_FV_REFLECT)right.u[a]*=-1;}
#ifdef GFS_NATIVE_KEP
            /* Reconstruction determines physical qdot and the compression
             * stress; the EC core itself uses the saved cell means. */
            double reconstructed=0,delta=0;
            point->left=w->primitive[i];point->right=j>=0?w->primitive[j]:w->primitive[i];
            if(j<0)for(int a=0;a<d;a++)if(mirrored(f->mirror,a)&&face_bc(m,f,a)==GFS_FV_REFLECT)point->right.u[a]*=-1;
            for(int a=0;a<d;a++){
                reconstructed+=(right.u[a]-left.u[a])*f->normal[a];
                delta+=(point->right.u[a]-point->left.u[a])*f->normal[a];
            }
            double compression=delta<0&&reconstructed<0?fmin(-delta,-reconstructed):0;
#ifdef GFS_NATIVE_CENTROID_STRAIN
            compression*=strain_scale;
#else
            double hi=pow(m->geometry[i].V,1./d),hj=j>=0?pow(m->geometry[j].V,1./d):hi;
            compression*=.5*(hi+hj)/f->distance;
#endif
            double rhobar=.5*(point->left.rho+point->right.rho);
            double sound=.5*(sqrt(m->gamma*point->left.p/point->left.rho)+sqrt(m->gamma*point->right.p/point->right.rho));
            point->compression=rhobar*(sound*compression+2*compression*compression);point->delta=delta;
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
            point->original_compression=point->compression;
#endif
#ifdef GFS_NATIVE_KEP_FILTER
            double fraction=j>=0?fmax(w->stress_fraction[i],w->stress_fraction[j]):w->stress_fraction[i];
            point->compression*=fraction;
#ifdef GFS_KEP_STRESS_DIAGNOSTICS
            point->fraction=fraction;
#endif
#endif
            point->density=gfs_kep_logmean(point->left.rho,point->right.rho);
            point->pressure=rhobar/(point->left.rho/(2*point->left.p)+point->right.rho/(2*point->right.p));
            if(!(point->density>0)||!isfinite(point->compression)){bad=1;continue;}
#else
            point->fan=gfs_fv_fan(d,left,right,f->normal,m->gamma);
            if(!point->fan.valid){fprintf(stderr,"FV invalid fan dim=%d i=%d j=%d rho=(%.17g %.17g) P=(%.17g %.17g)\n",d,i,j,left.rho,right.rho,left.p,right.p);bad=1;continue;}fallbacks+=point->fan.hll;
#endif
            point->weight=g->weight;point->inv2d=j>=0?1/(2*f->distance):0;
            if(j<0){
                point->wq=0;
                if(d==1&&face_bc(m,f,0)==GFS_FV_MOVING_TRANSMIT){
                    int inner=i==0?1:n-2;
                    point->wq=(1.5*w->velocity[i][0]-.5*w->velocity[inner][0])*f->normal[0];
                }
            }
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
static double norm(const GfsFvWorkspace *w,const double *r){
    double value=0;for(int i=0;i<w->n;i++){
        if(!isfinite(r[i])||!isfinite(w->scale[i])||!(w->scale[i]>0))return HUGE_VAL;
        value=fmax(value,fabs(r[i])/w->scale[i]);
    }return value;
}
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
#ifdef GFS_NATIVE_IC0
static int ic0_error(GfsFvWorkspace *w,const char *stage,int row,int column,double value){
    w->ic_failure_row=row;
    fprintf(stderr,"Native unshifted IC0 failure stage=%s row=%d column=%d value=%.17g; no pivot replacement or fallback\n",stage,row,column,value);return 1;
}
static int prepare_ic0(GfsFvWorkspace *w,const GfsFvMesh *m){
    double mark=omp_get_wtime();w->ic_nnz=0;
    if(w->ne>w->ic_capacity){
        w->ic_capacity=2*w->ne+16;free(w->ic_column);free(w->ic_matrix);free(w->ic_lower);
        w->ic_column=memory(w->ic_capacity*sizeof(int));w->ic_matrix=memory(w->ic_capacity*sizeof(double));w->ic_lower=memory(w->ic_capacity*sizeof(double));
    }
    /* Native row order, original lower sparsity only. Parallel/image faces
     * for the same owned pair are summed, without changing edge fluxes. */
    for(int i=0;i<m->n;i++){
        int start=w->ic_nnz;w->ic_offset[i]=start;
        for(int k=w->offset[i];k<w->offset[i+1];k++){
            FvEdge *e=w->edge+w->incident[k];if(e->j<0)continue;int j=e->i==i?e->j:e->i;
            if(j>=i)continue;
            double a=e->jacobian;if(!isfinite(a)||a<0){w->ic_assembly_seconds=omp_get_wtime()-mark;return ic0_error(w,"assembly",i,j,a);}
            if(a==0)continue;
            int found=start;while(found<w->ic_nnz&&w->ic_column[found]!=j)found++;
            if(found==w->ic_nnz){w->ic_column[found]=j;w->ic_matrix[found]=0;w->ic_nnz++;}
            w->ic_matrix[found]-=a;
        }
        for(int k=start+1;k<w->ic_nnz;k++){
            int column=w->ic_column[k],p=k;double value=w->ic_matrix[k];
            while(p>start&&w->ic_column[p-1]>column){w->ic_column[p]=w->ic_column[p-1];w->ic_matrix[p]=w->ic_matrix[p-1];p--;}
            w->ic_column[p]=column;w->ic_matrix[p]=value;
        }
    }
    w->ic_offset[m->n]=w->ic_nnz;w->ic_assembly_seconds=omp_get_wtime()-mark;mark=omp_get_wtime();w->ic_min_pivot=HUGE_VAL;
    for(int i=0;i<m->n;i++){
        double diagonal=w->diagonal[i];
        for(int p=w->ic_offset[i];p<w->ic_offset[i+1];p++){
            int j=w->ic_column[p],a=w->ic_offset[i],b=w->ic_offset[j];double value=w->ic_matrix[p];
            while(a<p&&b<w->ic_offset[j+1]){
                int ka=w->ic_column[a],kb=w->ic_column[b];
                if(ka<kb)a++;else if(kb<ka)b++;else{value-=w->ic_lower[a]*w->ic_pivot[ka]*w->ic_lower[b];a++;b++;}
            }
            value/=w->ic_pivot[j];
            if(!isfinite(value)){w->ic_factor_seconds=omp_get_wtime()-mark;return ic0_error(w,"lower",i,j,value);}
            w->ic_lower[p]=value;diagonal-=value*value*w->ic_pivot[j];
        }
        if(!(diagonal>0)||!isfinite(diagonal)){w->ic_factor_seconds=omp_get_wtime()-mark;return ic0_error(w,"pivot",i,i,diagonal);}
        w->ic_pivot[i]=diagonal;w->ic_min_pivot=fmin(w->ic_min_pivot,diagonal);
    }
    w->ic_factor_seconds=omp_get_wtime()-mark;return 0;
}
static int ic0_precondition(GfsFvWorkspace *w,const GfsFvMesh *m){
    double mark=omp_get_wtime();int bad=-1;w->ic_applications++;
    for(int i=0;i<m->n;i++){
        double value=w->r[i];for(int p=w->ic_offset[i];p<w->ic_offset[i+1];p++)value-=w->ic_lower[p]*w->pre[w->ic_column[p]];
        w->pre[i]=value;if(!isfinite(value)){bad=i;break;}
    }
    if(bad<0)for(int i=0;i<m->n;i++){w->pre[i]/=w->ic_pivot[i];if(!isfinite(w->pre[i])){bad=i;break;}}
    if(bad<0)for(int i=m->n-1;i>=0;i--){
        double value=w->pre[i];if(!isfinite(value)){bad=i;break;}
        for(int p=w->ic_offset[i];p<w->ic_offset[i+1];p++)w->pre[w->ic_column[p]]-=w->ic_lower[p]*value;
    }
    w->ic_apply_seconds+=omp_get_wtime()-mark;
    return bad>=0?ic0_error(w,"triangular",bad,bad,w->pre[bad]):0;
}
static void ic0_stats(const GfsFvWorkspace *w,GfsFvRhsStats *s){
    s->preconditioner_assembly_seconds=w->ic_assembly_seconds;s->preconditioner_factor_seconds=w->ic_factor_seconds;s->preconditioner_apply_seconds=w->ic_apply_seconds;
    s->minimum_preconditioner_pivot=isfinite(w->ic_min_pivot)?w->ic_min_pivot:0;s->preconditioner_nonzeros=w->ic_nnz+w->n;s->preconditioner_applications=w->ic_applications;s->preconditioner_failure_row=w->ic_failure_row;
}
#endif
#ifdef GFS_NATIVE_BLOCK_JACOBI
static void prepare_blocks(GfsFvWorkspace *w,const GfsFvMesh *m){
#pragma omp parallel for if(m->n>1024) schedule(static)
    for(int i=0;i<m->n;i++){
        int best=-1;double largest=0;
        for(int k=w->offset[i];k<w->offset[i+1];k++){
            FvEdge *e=w->edge+w->incident[k];if(e->j<0)continue;int j=e->i==i?e->j:e->i,seen=0;
            for(int t=w->offset[i];t<k;t++){FvEdge *old=w->edge+w->incident[t];if(old->j>=0&&(old->i==i?old->j:old->i)==j){seen=1;break;}}
            if(seen)continue;double a=0;
            for(int t=k;t<w->offset[i+1];t++){FvEdge *f=w->edge+w->incident[t];if(f->j>=0&&(f->i==i?f->j:f->i)==j)a+=f->jacobian;}
            if(a>largest||(a==largest&&a>0&&(best<0||m->cell[j].id<m->cell[best].id))){largest=a;best=j;}
        }
        w->strongest[i]=best;w->block_a[i]=largest;
    }
#pragma omp parallel for if(m->n>1024) schedule(static)
    for(int i=0;i<m->n;i++){
        int j=w->strongest[i];w->partner[i]=j>=0&&w->strongest[j]==i?j:-1;
        double b=m->cell[i].mu;
        for(int k=w->offset[i];k<w->offset[i+1];k++){
            FvEdge *f=w->edge+w->incident[k];if(f->j<0)continue;int neighbor=f->i==i?f->j:f->i;
            if(neighbor!=w->partner[i])b+=f->jacobian;
        }w->block_b[i]=b;
    }
}
static int block_precondition(GfsFvWorkspace *w,const GfsFvMesh *m){
    int bad=0;
#pragma omp parallel for if(m->n>1024) schedule(static) reduction(max:bad)
    for(int i=0;i<m->n;i++){
        int j=w->partner[i];
        if(j<0)w->pre[i]=w->r[i]/w->diagonal[i];
        else{
            double a=w->block_a[i],bi=w->block_b[i],bj=w->block_b[j],scale=fmax(a,fmax(bi,bj));
            double A=a/scale,B=bi/scale,C=bj/scale,det=A*(B+C)+B*C;
            /* a*(bi+bj)+bi*bj, never (a+bi)*(a+bj)-a*a. */
            w->pre[i]=(A*(w->r[i]+w->r[j])+C*w->r[i])/(scale*det);
        }
        if(!isfinite(w->pre[i]))bad=1;
    }return bad;
}
#endif
static int pcg(GfsFvWorkspace *w,const GfsFvMesh *m){
    w->pcg_iterations=0;w->pcg_failure=0;w->pcg_last_relative=1;w->pcg_stop_reason=0;w->pcg_true_replacements=0;w->pcg_first_candidate_true=NAN;
#ifdef GFS_PCG_TRUE_HISTORY
    for(int k=0;k<=64;k++)w->pcg_true_history[k]=w->pcg_recursive_history[k]=w->pcg_relative_history[k]=NAN;
    w->pcg_true_history[0]=w->pcg_recursive_history[0]=norm(w,w->R);w->pcg_relative_history[0]=1;
#endif
    int n=m->n;for(int i=0;i<n;i++){w->direction[i]=0;w->r[i]=-w->R[i];}
#ifdef GFS_NATIVE_IC0
    if(prepare_ic0(w,m)||ic0_precondition(w,m)){w->pcg_failure=5;return -1;}
#elif defined(GFS_NATIVE_BLOCK_JACOBI)
    prepare_blocks(w,m);if(block_precondition(w,m)){w->pcg_failure=5;return -1;}
#else
    for(int i=0;i<n;i++)w->pre[i]=w->r[i]/w->diagonal[i];
#endif
    for(int i=0;i<n;i++)w->p[i]=w->pre[i];
    double rz=product_sum(w->r,w->pre,n),initial=rz;if(initial==0)return 0;
    for(int step=1;step<=64;step++){
        w->pcg_iterations=step;
        multiply(w,m,w->p,w->Ap);double den=product_sum(w->p,w->Ap,n);if(!(den>0)){w->pcg_failure=1;return -1;}double a=rz/den;
        for(int i=0;i<n;i++){w->direction[i]+=a*w->p[i];w->r[i]-=a*w->Ap[i];}
#ifdef GFS_NATIVE_IC0
        if(ic0_precondition(w,m)){w->pcg_failure=5;return -1;}
#elif defined(GFS_NATIVE_BLOCK_JACOBI)
        if(block_precondition(w,m)){w->pcg_failure=5;return -1;}
#else
        for(int i=0;i<n;i++)w->pre[i]=w->r[i]/w->diagonal[i];
#endif
        double next=product_sum(w->r,w->pre,n);
        w->pcg_last_relative=sqrt(fmax(0,next/initial));
#ifdef GFS_PCG_TRUE_HISTORY
        /* Passive saved-system audit only. Flux[] is scratch; this does
         * not alter A, b, the Krylov recurrence, or stopping decisions. */
        if(!residual(w,m,w->direction,w->trial,0))w->pcg_true_history[step]=norm(w,w->trial);
        w->pcg_recursive_history[step]=norm(w,w->r);w->pcg_relative_history[step]=w->pcg_last_relative;
#endif
        if(next<=1e-24*initial){
#ifdef GFS_NATIVE_KEP
            /* An L2 recursive residual is not the registered scaled true
             * flux residual. Strong pressure contrasts exposed premature
             * convergence. Check that residual and reliably replace it
             * within the SAME 64-iteration cap and same SPD system. */
            if(residual(w,m,w->direction,w->trial,0)){w->pcg_failure=2;return -1;}
            double true_norm=norm(w,w->trial);
            if(isnan(w->pcg_first_candidate_true))w->pcg_first_candidate_true=true_norm;
            if(true_norm<=1e-12){w->pcg_stop_reason=1;return step;}
            w->pcg_true_replacements++;
            for(int i=0;i<n;i++)w->r[i]=-w->trial[i];
#ifdef GFS_NATIVE_IC0
            if(ic0_precondition(w,m)){w->pcg_failure=5;return -1;}
#elif defined(GFS_NATIVE_BLOCK_JACOBI)
            if(block_precondition(w,m)){w->pcg_failure=5;return -1;}
#else
            for(int i=0;i<n;i++)w->pre[i]=w->r[i]/w->diagonal[i];
#endif
            for(int i=0;i<n;i++)w->p[i]=w->pre[i];
            rz=product_sum(w->r,w->pre,n);if(!(rz>0)){w->pcg_failure=3;return -1;}continue;
#else
            return step;
#endif
        }
        for(int i=0;i<n;i++)w->p[i]=w->pre[i]+(next/rz)*w->p[i];rz=next;
    }
#ifdef GFS_NATIVE_KEP
    /* The recursive preconditioned norm is only a trigger. The registered
     * physical accuracy criterion also decides the last iterate at cap64. */
    if(residual(w,m,w->direction,w->trial,0)){w->pcg_failure=2;return -1;}
    if(norm(w,w->trial)<=1e-12){w->pcg_stop_reason=2;return 64;}
#endif
    w->pcg_failure=4;return -1;
}
#ifndef GFS_NATIVE_KEP
static double merit(const GfsFvMesh *m,const double *R){double value=0;for(int i=0;i<m->n;i++)value+=R[i]*R[i]/m->cell[i].mu;return value;}
#endif
static int solve(GfsFvWorkspace *w,const GfsFvMesh *m,GfsFvRhsStats *s){
    w->pcg_stop_reason=0;w->pcg_true_replacements=0;w->pcg_first_candidate_true=NAN;
    s->preconditioner_failure_row=-1;
#ifdef GFS_NATIVE_IC0
    w->ic_assembly_seconds=w->ic_factor_seconds=w->ic_apply_seconds=w->ic_min_pivot=0;w->ic_nnz=w->ic_applications=0;w->ic_failure_row=-1;
#endif
    memset(w->z,0,m->n*sizeof(double));
    for(int i=0;i<m->n;i++){
        double h0=sqrt(m->cell[i].m0/m->cell[i].mu),speed=sqrt(dot3(w->primitive[i].u,w->primitive[i].u,m->dim))+sqrt(m->gamma*w->primitive[i].p/w->primitive[i].rho);
        w->scale[i]=m->cell[i].mu*h0*speed;if(!(w->scale[i]>0))return 1;
    }
#ifdef GFS_NATIVE_KEP
    if(residual(w,m,w->z,w->R,1))return 1;s->root_residual=norm(w,w->R);
    if(s->root_residual<=1e-12){s->pcg_stop_reason=0;return 0;}
    int count=pcg(w,m);s->cg=w->pcg_iterations;s->pcg_stop_reason=w->pcg_stop_reason;
#ifdef GFS_NATIVE_IC0
    ic0_stats(w,s);
#endif
    if(count<0){
        double initial_residual=s->root_residual;
        if(!residual(w,m,w->direction,w->trial,0))s->root_residual=norm(w,w->trial);
        fprintf(stderr,"PCG failure iterations=%d reason=%d initial_scaled_residual=%.17g final_true_scaled_residual=%.17g (1=curvature,2=flux,3=replacement,4=cap,5=preconditioner)\n",w->pcg_iterations,w->pcg_failure,initial_residual,s->root_residual);return 3;
    }
    memcpy(w->z,w->direction,m->n*sizeof(double));
    if(residual(w,m,w->z,w->R,0))return 1;s->root_residual=norm(w,w->R);
    return s->root_residual<=1e-12?0:2;
#else
    for(int iteration=0;iteration<=8;iteration++){
        if(residual(w,m,w->z,w->R,1))return 1;s->root_residual=norm(w,w->R);s->newton=iteration;
        if(s->root_residual<=1e-12)return 0;if(iteration==8)return 2;
        int count=pcg(w,m);s->cg+=w->pcg_iterations;if(count<0)return 3;
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
    memset(s,0,sizeof(*s));s->pcg_stop_reason=-1;s->preconditioner_failure_row=-1;allocate(w,m->n);double mark=omp_get_wtime();
    if(reconstruct(w,m))return 1;s->reconstruction_seconds=omp_get_wtime()-mark;mark=omp_get_wtime();
#ifdef GFS_NATIVE_GLOBAL_VELOCITY
    s->velocity_query_seconds=w->query_seconds;s->velocity_owner_changed=w->owner_changed;s->velocity_owner_switches=w->owner_switches;s->maximum_query_hops=w->query_max_hops;s->total_query_hops=w->query_hops;
#endif
    if(build_points(w,m,s))return 2;s->fan_seconds=omp_get_wtime()-mark;mark=omp_get_wtime();
    int code=solve(w,m,s);s->solve_seconds=omp_get_wtime()-mark;if(code){fprintf(stderr,"FV screened solve failed code=%d Newton=%d CG=%d residual=%.17g\n",code,s->newton,s->cg,s->root_residual);return 10+code;}
    mark=omp_get_wtime();int bad=0;
#pragma omp parallel for if(m->n>1024) schedule(static) reduction(max:bad)
    for(int e=0;e<w->ne;e++){
        FvEdge *f=w->edge+e;memset(f->flux,0,sizeof(f->flux));f->Vdot=f->wave=f->viscous=f->pressure_work=f->kinetic_boundary=f->stress_heat=f->entropy_dissipation=f->entropy_boundary=0;
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
            if(f->j>=0)f->pressure_work+=p->weight*(p->pressure+p->compression)*p->delta;
            else f->kinetic_boundary+=p->weight*(.5*dot3(p->left.u,p->left.u,m->dim)*flux[0]+(p->pressure+p->compression)*left);
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
        s->kinetic_boundary_flux+=w->edge[e].kinetic_boundary;
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

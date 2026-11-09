#ifndef GFS_MATERIAL_TRACE_CSR2D_H
#define GFS_MATERIAL_TRACE_CSR2D_H
#if defined(GFS_MATERIAL_DIRECT_WALLS) && !defined(GFS_MATERIAL_DIRECT_FV)
#error "GFS_MATERIAL_DIRECT_WALLS requires GFS_MATERIAL_DIRECT_FV"
#endif
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
/* Directed owner-sorted periodic faces with canonical local cell indices.
 * d and c must retain the actual periodic image separation and owner-relative
 * face centroid. This function does not infer MPI ownership or halo aliases.
 */
typedef struct {int i,j;double dx,dy,cx,cy,area;} GfsTraceDirected2D;
/* In 2D a power face is a segment perpendicular to d. Its endpoints are
 * c +/- (A/2)(-dy,dx)/|d|. A nearly cancelled midpoint does not bound the
 * rounding scale of the endpoints from which it was formed. These are local
 * coordinate lengths, not an absolute floor or a 3D face-area rule. */
static inline int gfs_trace_endpoint_scales2d(const GfsTraceDirected2D*e,
 double d,double scale[2]){
 if(!e||!scale||!(d>0)||!isfinite(d))return 0;
 scale[0]=fabs(e->cx)+(.5*e->area)*(fabs(e->dy)/d);
 scale[1]=fabs(e->cy)+(.5*e->area)*(fabs(e->dx)/d);
 return isfinite(scale[0])&&isfinite(scale[1]);
}
#include "gfs_material_plane.h"
/* General stationary planes are supported by the direct one-sided wall
 * flux. The legacy reflected-entry/adjoint representation remains y-only. */
static inline int gfs_material_stationary_wall2d(const GfsTraceDirected2D*e){
 if(!e||e->i!=e->j)return 0;
#ifdef GFS_MATERIAL_DIRECT_WALLS
 double d[2]={e->dx,e->dy},c[2]={e->cx,e->cy};
 return gfs_material_plane_bisector(2,d,c);
#else
 double d=hypot(e->dx,e->dy);
 return d>0&&isfinite(d)&&fabs(e->dx)<=256*DBL_EPSILON*d&&fabs(e->cy-.5*e->dy)<=256*DBL_EPSILON*d;
#endif
}

#include "gfs_material_centroid2d.h"
typedef struct {int cell,pad;double weight;
#ifdef GFS_MATERIAL_FACE_NORMAL_LIMIT
 double velocity_alpha; /* Prepared once per current-state face side. */
#endif
#ifdef GFS_MATERIAL_POSITIVE_FAN
 double reversible_weight;
#endif
} GfsTraceEntry;
typedef struct {int i,j,sl,nl,sr,nr;double area,nx,ny;
#ifdef GFS_MATERIAL_GRADIENT_PROJECT
 double thermo[4]; /* rhoL,pL,rhoR,pR prepared once per immutable RHS */
#endif
} GfsTraceFace2D;
typedef struct {int cells,faces,entries;GfsTraceFace2D *face;GfsTraceEntry *entry;} GfsTraceCsr2D;
static inline void gfs_trace_csr2d_free(GfsTraceCsr2D *t){if(!t)return;free(t->face);free(t->entry);memset(t,0,sizeof(*t));}
/* Unrestricted inverse-distance-squared least squares, identical to the
 * controlled residual-adjoint screen. A limiter/positivity policy is NOT
 * silently introduced. Invalid geometry fails rather than replacing a fit.
 */
/* Candidate only: face-length / generator-distance measures the angular
 * support of a 2D face. Its contribution vanishes with a disappearing face.
 * Use the same metric in the normal matrix and trace coefficients. */
static inline
#ifdef __CUDACC__
__host__ __device__
#endif
double gfs_trace_metric2d(const GfsTraceDirected2D *e,double d2){
#ifdef GFS_MATERIAL_FACE_MEASURE_LS
 return e->area/sqrt(d2);
#else
 (void)e;(void)d2;return 1;
#endif
}
/* n owned reconstruction rows, nt physical samples, ns total data slots.
 * External faces use a prescribed face-state slot (>=nt), supplied through
 * boundary_slot[k]. This slot is only a boundary datum and exchange ledger;
 * it is never advanced as an owned fluid cell. Interior reconstruction and
 * its C/CUDA trace representation are unchanged. */
static inline int gfs_trace_csr2d_build_coordinates_domain_planes(int n,int nt,int ns,int nf,const GfsTraceDirected2D *f,GfsTraceCsr2D *out,const unsigned char *reflected_y,const GfsTraceDirected2D *fit,const int*boundary_slot,int physical_planes){
 if(!out)return 0;memset(out,0,sizeof(*out));if(n<1||nt<n||ns<nt||nf<1||!f||!fit||(nt>n&&!boundary_slot))return 0;
 int *off=(int*)calloc((size_t)n+1,sizeof(int));double *inv=(double*)calloc((size_t)n*3,sizeof(double));
 GfsTraceCsr2D t={0};int ok=0,m=0,row=0,emission_failed=0;size_t count=0;
 int *face_key=NULL;size_t *entry_start=NULL;
 if(!off||!inv)goto done;
 for(int k=0;k<nf;k++){
  const GfsTraceDirected2D *a=f+k;double d2=a->dx*a->dx+a->dy*a->dy;
  int wall=reflected_y?reflected_y[k]:0;
  if(wall>1)goto done;
  if(a->i<0||a->i>=n||a->j<0||a->j>=nt||(!wall&&a->i==a->j)||
     (k&&f[k-1].i>a->i)||!(d2>0)||!isfinite(d2)||!(a->area>0)||
     !isfinite(a->area)||!isfinite(a->cx)||!isfinite(a->cy))goto done;
  if(a->j>=n&&(wall||boundary_slot[k]<nt||boundary_slot[k]>=ns))goto done;
  /* A supplied physical wall plane has a normal and centroid independent
   * of a mirror-site bisector. Its underlying power geometry must already
   * have been certified by the caller. Legacy reflected walls retain their
   * bisector check. No cell/face geometry is changed by reconstruction. */
#ifndef GFS_MATERIAL_DIRECT_WALLS
  if(physical_planes)goto done;
#endif
  if(wall&&(a->i!=a->j||(!physical_planes&&!gfs_material_stationary_wall2d(a))))goto done;
  const GfsTraceDirected2D *z=fit+k;double dx=z->dx,dy=z->dy,fd2=dx*dx+dy*dy;
  if(!(fd2>0)||!isfinite(fd2))goto done;
  double metric=gfs_trace_metric2d(z,fd2);
#ifdef GFS_MATERIAL_DIRECT_WALLS
  /* A physical wall is a flux boundary, not an independent scalar sample.
   * Fit the interior field one-sided; the wall Riemann state is reflected
   * only after extrapolation. Otherwise an affine pressure is flattened. */
  if(wall)metric=0;
#endif
  off[a->i+1]++;inv[3*a->i]+=metric*dx*dx/fd2;
  inv[3*a->i+1]+=metric*dx*dy/fd2;inv[3*a->i+2]+=metric*dy*dy/fd2;
  if(a->i<a->j)m++;
#ifdef GFS_MATERIAL_DIRECT_WALLS
  if(wall)m++;
#endif
 }
 for(int i=0;i<n;i++){
  off[i+1]+=off[i];double a=inv[3*i],b=inv[3*i+1],c=inv[3*i+2],det=a*c-b*b;
  if(!(det>64*DBL_EPSILON*a*c)||!isfinite(det))goto done;
  inv[3*i]=c/det;inv[3*i+1]=-b/det;inv[3*i+2]=a/det;
 }
 /* Integer-only canonical offsets make each face output independent.
  * All geometric validation and per-face floating-point sums retain
  * their original order. This storage is local to this one build. */
 if(m<1)goto done;
 face_key=(int*)malloc((size_t)m*sizeof(int));
 entry_start=(size_t*)malloc(((size_t)m+1)*sizeof(size_t));
 if(!face_key||!entry_start)goto done;
 /* Require exactly one reverse directed image for every undirected face. */
 for(int k=0;k<nf;k++){
  const GfsTraceDirected2D *a=f+k;int reverse=0;
  if(reflected_y&&reflected_y[k]){
#ifdef GFS_MATERIAL_DIRECT_WALLS
   face_key[row]=k;entry_start[row++]=count;
   count+=(size_t)(off[a->i+1]-off[a->i]+1);
#endif
   continue;
  }
  if(a->j>=n){face_key[row]=k;entry_start[row++]=count;
   count+=(size_t)(off[a->i+1]-off[a->i]+2);continue;}
  for(int h=off[a->j];h<off[a->j+1];h++){
   const GfsTraceDirected2D *b=f+h;
   double scale=sqrt(a->dx*a->dx+a->dy*a->dy),tol=256*DBL_EPSILON*scale;
   if(b->j==a->i&&fabs(b->dx+a->dx)<=tol&&fabs(b->dy+a->dy)<=tol){
    /* In 2D area is a face LENGTH. Its rounding scale includes
     * the face itself, not only the possibly much shorter site distance.
     * Keep image/centroid matching unchanged. This is not a 3D area rule. */
    double length_tol=256*DBL_EPSILON*fmax(scale,fmax(a->area,b->area));
    /* Preserve the epsilon multiplier and separation/centroid lower bound,
     * while accounting for cancellation in the endpoint midpoint. */
    double sa[2],sb[2];
    if(!gfs_trace_endpoint_scales2d(a,scale,sa)||!gfs_trace_endpoint_scales2d(b,scale,sb))goto done;
    double cx_tol=256*DBL_EPSILON*fmax(scale,fmax(sa[0],sb[0]));
    double cy_tol=256*DBL_EPSILON*fmax(scale,fmax(sa[1],sb[1]));
    if(fabs(b->cx-(a->cx-a->dx))>cx_tol||fabs(b->cy-(a->cy-a->dy))>cy_tol||fabs(b->area-a->area)>length_tol){
#ifdef GFS_MATERIAL_FIXED_SOLID
     continue;
#else
     goto done;
#endif
    }
    reverse++;
   }
  }
  if(reverse!=1)goto done;
  if(a->i<a->j){face_key[row]=k;entry_start[row++]=count;
   count+=(size_t)(off[a->i+1]-off[a->i]+off[a->j+1]-off[a->j]+2);}
 }
 if(count>INT_MAX||row!=m)goto done;entry_start[m]=count;
 t.face=(GfsTraceFace2D*)malloc((size_t)m*sizeof(*t.face));t.entry=(GfsTraceEntry*)malloc(count*sizeof(*t.entry));if(!t.face||!t.entry)goto done;
 t.cells=ns;t.faces=m;t.entries=(int)count;
#pragma omp parallel for schedule(static) reduction(|:emission_failed)
 for(int face_row=0;face_row<m;face_row++){
  int k=face_key[face_row],slot=(int)entry_start[face_row];
  const GfsTraceDirected2D *a=f+k;double d=hypot(a->dx,a->dy);
  GfsTraceFace2D *b=t.face+face_row;b->i=a->i;b->j=a->j>=n?boundary_slot[k]:a->j;b->area=a->area;b->nx=a->dx/d;b->ny=a->dy/d;
  int sides=a->j>=n?1:2;
#ifdef GFS_MATERIAL_DIRECT_WALLS
  if(reflected_y&&reflected_y[k])sides=1;
#endif
  for(int side=0;side<sides;side++){
   int owner=side?a->j:a->i,start=slot++;double x=fit[k].cx-(side?fit[k].dx:0),y=fit[k].cy-(side?fit[k].dy:0),sum=0;
   double gx=x*inv[3*owner]+y*inv[3*owner+1],gy=x*inv[3*owner+1]+y*inv[3*owner+2];
   for(int h=off[owner];h<off[owner+1];h++){
    const GfsTraceDirected2D *z=fit+h;double d2=z->dx*z->dx+z->dy*z->dy,w=gfs_trace_metric2d(z,d2)*(gx*z->dx+gy*z->dy)/d2;
#ifdef GFS_MATERIAL_DIRECT_WALLS
    if(reflected_y&&reflected_y[h])w=0;
#endif
    if(!isfinite(w))emission_failed=1;
    t.entry[slot++]=(GfsTraceEntry){z->j,(reflected_y&&reflected_y[h])?(2|(z->dy>0?4:0)):0,w};sum+=w;
   }
   t.entry[start]=(GfsTraceEntry){owner,0,1-sum};
   if(side){b->sr=start;b->nr=slot-start;}else{b->sl=start;b->nl=slot-start;}
  }
  if(a->j>=n){
   b->sr=slot;b->nr=1;t.entry[slot++]=(GfsTraceEntry){b->j,0,1};
  }
#ifdef GFS_MATERIAL_DIRECT_WALLS
  /* One interior trace; the wall primitive constructs its reflected state. */
  if(sides==1&&a->j<n){b->sr=b->sl;b->nr=b->nl;}
#endif
  if(slot!=(int)entry_start[face_row+1])emission_failed=1;
 }
 if(emission_failed)goto done;
#ifdef GFS_MATERIAL_FACE_NORMAL_LIMIT
 for(int j=0;j<t.entries;j++)t.entry[j].velocity_alpha=NAN;
#endif
#ifdef GFS_MATERIAL_POSITIVE_FAN
 for(int j=0;j<t.entries;j++)t.entry[j].reversible_weight=NAN;
#endif
 *out=t;memset(&t,0,sizeof(t));ok=1;
 done:free(face_key);free(entry_start);free(off);free(inv);gfs_trace_csr2d_free(&t);return ok;
}
/* Existing APIs keep the stricter reflected-bisector contract. */
static inline int gfs_trace_csr2d_build_coordinates_domain(int n,int nt,int ns,int nf,
 const GfsTraceDirected2D*f,GfsTraceCsr2D*out,const unsigned char*reflected,
 const GfsTraceDirected2D*fit,const int*boundary_slot){
 return gfs_trace_csr2d_build_coordinates_domain_planes(n,nt,ns,nf,f,out,reflected,fit,boundary_slot,0);
}
/* Original closed-domain API and arithmetic. */
static inline int gfs_trace_csr2d_build_coordinates(int n,int nf,const GfsTraceDirected2D*f,GfsTraceCsr2D*out,const unsigned char*reflected,const GfsTraceDirected2D*fit){
 return gfs_trace_csr2d_build_coordinates_domain(n,n,n,nf,f,out,reflected,fit,NULL);
}
static inline int gfs_trace_csr2d_build_boundaries(int n,int nf,const GfsTraceDirected2D*f,GfsTraceCsr2D*out,const unsigned char*reflected){
#ifdef GFS_MATERIAL_CENTROID_TRACE
 GfsTraceDirected2D *fit=gfs_trace_centroid_coordinates2d(n,nf,f,reflected);
 if(!fit){if(out)memset(out,0,sizeof(*out));return 0;}
 int ok=gfs_trace_csr2d_build_coordinates(n,nf,f,out,reflected,fit);free(fit);return ok;
#else
 return gfs_trace_csr2d_build_coordinates(n,nf,f,out,reflected,f);
#endif
}
/* pad bit2 denotes a y-reflected image; serialized periodic entries remain0. */
#ifdef __CUDACC__
#define GFS_TRACE_HD static __host__ __device__ inline
#else
#define GFS_TRACE_HD static inline
#endif
GFS_TRACE_HD double gfs_trace_velocity_sign(const GfsTraceEntry *e,int axis){
 return axis==1&&(e->pad&2)?-1.:1.;
}
GFS_TRACE_HD int gfs_trace_owner_entry(const GfsTraceEntry *e,int owner){
 return e->cell==owner&&!(e->pad&2);
}
#undef GFS_TRACE_HD
static inline int gfs_trace_csr2d_build(int n,int nf,const GfsTraceDirected2D *f,GfsTraceCsr2D *out){
 return gfs_trace_csr2d_build_boundaries(n,nf,f,out,NULL);
}
#endif

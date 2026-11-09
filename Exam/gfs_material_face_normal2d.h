#ifndef GFS_MATERIAL_FACE_NORMAL2D_H
#define GFS_MATERIAL_FACE_NORMAL2D_H
#include "gfs_material_trace_csr2d.h"
#include "gfs_material_face_response.h"
#include "gfs_material_ordered_jump.h"
#if defined(GFS_MATERIAL_POSITIVE_FAN) && !defined(GFS_MATERIAL_FACE_NORMAL_LIMIT)
#error "Positive fan candidate requires face-normal candidate ABI"
#endif
#ifdef GFS_MATERIAL_FACE_NORMAL_LIMIT
/* Direct material FV uses the reconstructed normal states only; it has no
 * separate reversible transpose that must remain unlimited. */
#if !defined(GFS_MATERIAL_DIRECT_FV) && (!defined(GFS_MATERIAL_SPLIT_REVERSIBLE) || GFS_MATERIAL_SPLIT_REVERSIBLE != 1)
#error "Face-normal candidate requires the unlimited reversible trace"
#endif
/* Frozen face-side map, shared by gather, transpose force, constraint actions,
 * and damping bound. Recompute at each RHS; never differentiate/re-limit a
 * constraint vector. Reflection follows the same ghost map as reconstruction. */
/* wall_normal is optional prescribed normal speed per trace face. Only
 * one-sided direct-wall faces (i==j) read it. Limit their reflected states
 * in the boundary frame, as used by the moving-wall Riemann flux.
 * NULL preserves the existing stationary-wall API and arithmetic. */
static inline int gfs_material_face_normal_prepare2d_selected(
 GfsTraceCsr2D*t,const GfsFvPrimitive*c,const double*wall_normal,int walls_only){
 if(!t||!c)return 0;
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr};
  if(walls_only && f->i!=f->j)continue;
#ifdef GFS_MATERIAL_DIRECT_WALLS
  if(wall_normal&&f->i==f->j&&!isfinite(wall_normal[k]))return 0;
#endif
#ifdef GFS_MATERIAL_ORDERED_JUMP
  double center[2],increment[2]={0};
  for(int s=0;s<2;s++){
   int owner=s?f->j:f->i;center[s]=f->nx*c[owner].u[0]+f->ny*c[owner].u[1];
#ifdef GFS_MATERIAL_DIRECT_WALLS
   if(f->i==f->j){
    if(wall_normal)center[s]-=wall_normal[k];
    if(s==1)center[s]=-center[s];
   }
#endif
   if(!isfinite(center[s]))return 0;
   for(int j=starts[s];j<starts[s]+counts[s];j++){
    const GfsTraceEntry*e=t->entry+j;
    double v=f->nx*c[e->cell].u[0]+f->ny*gfs_trace_velocity_sign(e,1)*c[e->cell].u[1];
#ifdef GFS_MATERIAL_DIRECT_WALLS
    if(f->i==f->j){
     if(wall_normal)v-=wall_normal[k];
     if(s==1)v=-v;
    }
#endif
    if(!isfinite(v))return 0;
#ifdef GFS_MATERIAL_COMMON_POSITIVE
    increment[s]+=e->reversible_weight*(v-center[s]);
#else
    increment[s]+=e->weight*(v-center[s]);
#endif
   }
  }
  double mean_jump=center[1]-center[0],change=increment[1]-increment[0],a=1;
  if(!gfs_material_ordered_jump_alpha(mean_jump,change,&a))return 0;
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++)t->entry[j].velocity_alpha=a;
#else
  for(int s=0;s<2;s++){
   int owner=s?f->j:f->i;double center=f->nx*c[owner].u[0]+f->ny*c[owner].u[1];
   double lo=center,hi=center,delta=0;if(!isfinite(center))return 0;
   for(int j=starts[s];j<starts[s]+counts[s];j++){
    const GfsTraceEntry*e=t->entry+j;
    double v=f->nx*c[e->cell].u[0]+f->ny*gfs_trace_velocity_sign(e,1)*c[e->cell].u[1];
    if(!isfinite(v))return 0;lo=fmin(lo,v);hi=fmax(hi,v);delta+=e->weight*(v-center);
   }
   if(!isfinite(delta))return 0;
   double a=1;if(delta>0)a=fmin(1,(hi-center)/delta);else if(delta<0)a=fmin(1,(lo-center)/delta);
   a=fmax(0,a);if(!isfinite(a))return 0;
   for(int j=starts[s];j<starts[s]+counts[s];j++)t->entry[j].velocity_alpha=a;
  }
#endif
 }
 return 1;
}
static inline int gfs_material_face_normal_prepare2d_wall_velocity(
 GfsTraceCsr2D*t,const GfsFvPrimitive*c,const double*wall_normal){
 return gfs_material_face_normal_prepare2d_selected(t,c,wall_normal,0);
}
/* Within one immutable-state RHS, only the recovered wall speeds changed.
 * Requires an earlier successful full preparation with unchanged trace and
 * cell velocities. Interior/exterior entries keep those prepared values. */
static inline int gfs_material_face_normal_refresh_walls2d(
 GfsTraceCsr2D*t,const GfsFvPrimitive*c,const double*wall_normal){
 return gfs_material_face_normal_prepare2d_selected(t,c,wall_normal,1);
}
static inline int gfs_material_face_normal_prepare2d(GfsTraceCsr2D*t,const GfsFvPrimitive*c){
 return gfs_material_face_normal_prepare2d_wall_velocity(t,c,NULL);
}
#endif
#ifdef __CUDACC__
static __host__ __device__ inline
#else
static inline
#endif
GfsMaterialTraceWeights gfs_material_entry_weights(const GfsTraceEntry*e,int owner,const double*alpha){
#ifdef GFS_MATERIAL_FACE_NORMAL_LIMIT
 double a[3]={alpha?alpha[0]:1,alpha?alpha[1]:1,e->velocity_alpha};
#ifdef GFS_MATERIAL_COMMON_POSITIVE
 GfsMaterialTraceWeights w=gfs_material_trace_weights(e->reversible_weight,gfs_trace_owner_entry(e,owner),a);
#else
 GfsMaterialTraceWeights w=gfs_material_trace_weights(e->weight,gfs_trace_owner_entry(e,owner),a);
#endif
#ifdef GFS_MATERIAL_POSITIVE_FAN
 w.reversible=e->reversible_weight;
#endif
 return w;
#else
 return gfs_material_trace_weights(e->weight,gfs_trace_owner_entry(e,owner),alpha);
#endif
}
#endif

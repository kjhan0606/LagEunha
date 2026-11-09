#ifndef GFS_MATERIAL_DAMPING_BOUND2D_H
#define GFS_MATERIAL_DAMPING_BOUND2D_H
#include "gfs_material_constraint_ops2d.h"
#include "gfs_material_damping_bound.h"
/* Frozen velocity damping: M^-1/2 K M^-1/2 = sum_f b_f b_f^T.
 * Group each b_f by cell vector, AFTER combining duplicate trace entries.
 * Block Gershgorin: lambda_max <= max_i sum_f |b_fi| sum_j |b_fj|.
 * Vector norms make this bound invariant under a common spatial rotation.
 * No global matrix/eigensolve. This is not the full hydro Jacobian bound.
 */
static inline int gfs_material_damping_bound2d(const GfsTraceCsr2D*t,
 const GfsMaterialFaceResponse*cache,const double*alpha,const double*mass,double*bound){
 if(!t||!cache||!mass||!bound||t->cells<1)return 0;
 int n=t->cells,ok=0;double *row=calloc(n,sizeof(double)),*v=calloc((size_t)n*2,sizeof(double));
 int *stamp=malloc(n*sizeof(int)),*ids=malloc(n*sizeof(int));
 if(!row||!v||!stamp||!ids)goto done;
 for(int i=0;i<n;i++){stamp[i]=-1;if(!(mass[i]>0)||!isfinite(mass[i]))goto done;}
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;const GfsMaterialFaceResponse*r=cache+k;
  double coeff=f->area*r->zleft*(r->zright/(r->zleft+r->zright));
  if(!(coeff>0)||!isfinite(coeff))goto done;
  int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr},used=0;
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
   const GfsTraceEntry*e=t->entry+j;int i=e->cell,owner=s?f->j:f->i;
   if(i<0||i>=n)goto done;
   if(stamp[i]!=k){stamp[i]=k;ids[used++]=i;v[2*i]=v[2*i+1]=0;}
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha?alpha+3*owner:NULL);
   v[2*i]+=(s?1:-1)*w.velocity*f->nx;
   v[2*i+1]+=(s?1:-1)*w.velocity*f->ny*gfs_trace_velocity_sign(e,1);
  }
  if(!gfs_material_damping_row(2,used,ids,v,mass,coeff,row))goto done;
 }
 *bound=0;for(int i=0;i<n;i++){if(!isfinite(row[i]))goto done;*bound=fmax(*bound,row[i]);}ok=1;
 done:free(row);free(v);free(stamp);free(ids);return ok;
}
#endif

/* Opt-in inverse after existing native GPU geometry. No force, U or x update. */
#include "../gfs_material_volume_closed.h"
int gfs_material_native_volume_correct(void*opaque,int stage,double dt,int iteration){
 SimParameters*s=opaque;GPUContext*c=gpu_get_context();
 GfsMaterialVolumeAudit*a=&gfs_material_native_volume_audit;
 if(!s||!c||!c->initialized||!c->tess_faces_on_device||NID(s)!=1||
    (stage!=1&&stage!=2)||a->phase!=stage||
    !a->row||a->n!=VORO_NP(s)||!(dt>0)||!isfinite(dt))return 0;
 int n=a->n,nt=c->h_parts.n_total;ParticleSoA*p=&c->h_parts;FaceCSR*f=&c->h_faces;
 int cylinder=gfs_material_cylinder_mode(),exterior=gfs_material_noh_native_mode();
 int anchored=cylinder||exterior;
 const GfsCylinderMaterialSites*b=cylinder?gfs_material_cylinder_sites(n,nt,p->indx):NULL;
 const int*aliases=NULL;int*rt_alias=NULL;unsigned char*reflection=NULL;
 double lx=Xmax_HydroExam(s)-Xmin_HydroExam(s),ly=Ymax_HydroExam(s)-Ymin_HydroExam(s);
 const char*rt=getenv("GFS_MATERIAL_RT_NATIVE");int wall=rt&&!strcmp(rt,"1");
 if(!anchored){
  if(wall){
   int count=0;const GfsMaterialRtOwner*pad=gfs_material_rt_owner_map(&count);
   rt_alias=malloc((size_t)nt*sizeof(int));reflection=malloc((size_t)nt);
   if(!rt_alias||!reflection||!gfs_material_rt_resolve_owners(n,nt,p->indx,p->x,p->y,count,pad,lx,Ymin_HydroExam(s),Ymax_HydroExam(s),rt_alias,reflection)){free(rt_alias);free(reflection);return 0;}
   aliases=rt_alias;ly=0;
  }else aliases=gpu_unified_owner_aliases(c,n,nt,lx,ly);
  if(!aliases){free(rt_alias);free(reflection);return 0;}
 }
 treevorostressrk4particletype*cell=(void*)VORORK4_TBP(s);
 if((cylinder&&!b)||!cell){free(rt_alias);free(reflection);return 0;}
 double*rhs=malloc((size_t)n*sizeof(double)),*delta=malloc((size_t)n*sizeof(double));
 int*owner=NULL,*resolved_offset=NULL;double*coefficient=NULL;
 GfsTraceDirected2D*resolved_edge=NULL;unsigned char*resolved_wall=NULL;
 GfsMaterialExterior2D exterior_geometry={0};
 const int*offset=f->face_offset;
 double*target_volume=malloc((size_t)n*sizeof(double));
 const double volume_tolerance=1e-11;
 int answer=0;double max_relative=0,start=MPI_Wtime();long long worst=-1;
 if(!rhs||!delta||!target_volume)goto done;
 for(int i=0;i<n;i++){
  if((long long)PINDX(cell+i)!=p->indx[i]||cell[i].weight.W!=p->weight_W[i])goto done;
  GfsMaterialVolumeRow key={0};key.id=p->indx[i];
  const GfsMaterialVolumeRow*r=bsearch(&key,a->row,n,sizeof(*a->row),gfs_material_volume_order);
  if(!r||(stage==2&&!r->second))goto done;
  double target=stage==1?r->volume+dt*r->rate[0]:
                        r->volume+.5*dt*(r->rate[0]+r->rate[1]);
  double actual=p->weight_volume64[i];
  if(!(target>0)||!isfinite(target)||!(actual>0)||!isfinite(actual))goto done;
  target_volume[i]=target;rhs[i]=target-actual;double relative=fabs(rhs[i])/target;
  if(relative>max_relative){max_relative=relative;worst=p->indx[i];}
 }
 {
  double t;int phase;if(!gfs_material_native_evaluation_time(&t,&phase))goto done;
  fprintf(stderr,"[MATERIAL_VOLUME_INVERSE] t=%.17g target_stage=%d evaluation_phase=%d iteration=%d max_relative=%.17g worst_id=%lld\n",
          t,stage,phase,iteration,max_relative,worst);
 }
 if(max_relative<=volume_tolerance){answer=1;goto done;}
 if(iteration>=8)goto done;
 gpu_download_face_csr(c,n,c->d_faces.n_faces_total,f->face_offset,f->c1x,f->c1y,
                      f->c2x,f->c2y,f->neighbor_idx,f->kp_idx,f->km_idx,f->is_ghost);
 int nf=f->face_offset[n];
 if(nf<1)goto done;
 owner=malloc((size_t)nf*sizeof(int));coefficient=malloc((size_t)nf*sizeof(double));
 if(!owner||!coefficient)goto done;
 if(cylinder){
 for(int i=0;i<n;i++)for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
  int j=f->neighbor_idx[k];
#ifdef GFS_MATERIAL_FIXED_SOLID
  if(j<=GFS_SOLID_TAG_BASE){owner[k]=i;coefficient[k]=0;continue;}
#endif
  if(j<0||j>=nt)goto done;
  owner[k]=cylinder?(j<n?j:b->alias[j]>=0?b->alias[j]:b->image[j]>=0?b->image[j]:-1):aliases[j];
  if(!cylinder&&(owner[k]<0||owner[k]>=n))goto done;
  if(!cylinder&&reflection&&reflection[j]){
   if(owner[k]!=i)goto done;coefficient[k]=0;continue;
  }
  double dx=cylinder?gfs_power_coordinate_delta(p->x,b->lowx,i,j):material_canonical_delta(p->x[i],p->x[aliases[j]],p->x[j],lx);
  double dy=cylinder?gfs_power_coordinate_delta(p->y,b->lowy,i,j):(ly>0?material_canonical_delta(p->y[i],p->y[aliases[j]],p->y[j],ly):p->y[j]-p->y[i]),d=hypot(dx,dy);
  double area=hypot(f->c2x[k]-f->c1x[k],f->c2y[k]-f->c1y[k]);
  if(!(d>0)||!isfinite(d)||!(area>=0)||!isfinite(area))goto done;
  coefficient[k]=area/(2*d);
 }
 }else if(exterior){
  /* Identical resolved exterior faces to the pressure-work adapter. All
   * prescribed sample weights remain fixed during this Newton correction. */
  if(!gfs_material_exterior2d_build(n,nt,nf,f->face_offset,f->neighbor_idx,
      p->x,p->y,p->weight_W,f->c1x,f->c1y,f->c2x,f->c2y,
      p->weight_volume64,&exterior_geometry))goto done;
  resolved_offset=calloc((size_t)n+1,sizeof(int));if(!resolved_offset)goto done;
  for(int k=0;k<exterior_geometry.directed_count;k++){
   const GfsTraceDirected2D*e=exterior_geometry.directed+k;
   if(e->i<0||e->i>=n||e->j<0||e->j>=nt||e->i==e->j||
      (k&&e->i<exterior_geometry.directed[k-1].i))goto done;
   resolved_offset[e->i+1]++;owner[k]=e->j<n?e->j:-1;
   coefficient[k]=e->area/(2*hypot(e->dx,e->dy));
  }
  for(int i=0;i<n;i++)resolved_offset[i+1]+=resolved_offset[i];
  offset=resolved_offset;
 }else{
  resolved_edge=malloc((size_t)nf*sizeof(*resolved_edge));
  resolved_wall=reflection?calloc((size_t)nf,1):NULL;
  resolved_offset=calloc((size_t)n+1,sizeof(int));
  if(!resolved_edge||!resolved_offset||(reflection&&!resolved_wall))goto done;
  int count=0;
  if(!gfs_material_resolved_faces2d(c,n,nt,lx,ly,aliases,reflection,resolved_edge,resolved_wall,&count))goto done;
  for(int k=0;k<count;k++){
   const GfsTraceDirected2D*e=resolved_edge+k;
   if(e->i<0||e->i>=n||e->j<0||e->j>=n||(k&&e->i<resolved_edge[k-1].i))goto done;
   resolved_offset[e->i+1]++;owner[k]=e->j;
   coefficient[k]=(resolved_wall&&resolved_wall[k])?0:e->area/(2*hypot(e->dx,e->dy));
  }
  for(int i=0;i<n;i++)resolved_offset[i+1]+=resolved_offset[i];
  offset=resolved_offset;
 }
 if(anchored&&!gfs_lu_prepare_basis(n,p->indx,offset,owner,coefficient))goto done;
 GfsVolumeLinearStats linear;
 if(!(anchored?gfs_material_volume_linear_correction(n,offset,owner,coefficient,rhs,delta,&linear):gfs_material_volume_closed_correction_scaled(n,p->indx,offset,owner,coefficient,rhs,delta,&linear,target_volume,volume_tolerance))){
  fprintf(stderr,"[MATERIAL_VOLUME_INVERSE_LINEAR] failed error=%d pivots=%d residual=%.17g\n",
          linear.error,linear.perturbed_pivots,linear.relative_residual);goto done;
 }
 for(int i=0;i<n;i++)if(!isfinite(cell[i].weight.W+delta[i]))goto done;
 for(int i=0;i<n;i++)cell[i].weight.W+=delta[i];
 fprintf(stderr,"[MATERIAL_VOLUME_INVERSE_LINEAR] iteration=%d rows=%d entries=%d residual=%.17g seconds=%.9g\n",
         iteration,n,linear.entries,linear.relative_residual,MPI_Wtime()-start);
 answer=2;
done:
 gfs_material_exterior2d_free(&exterior_geometry);
 free(target_volume);free(resolved_wall);free(resolved_edge);free(resolved_offset);
 free(reflection);free(rt_alias);free(coefficient);free(owner);free(delta);free(rhs);return answer;
}


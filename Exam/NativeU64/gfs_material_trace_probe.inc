/* Canonical periodic displacement: subtract owner coordinates before
 * adding the integer image shift. Padding coordinates were already checked
 * by the existing owner resolver. Do not recompute a minimum image: multiple
 * periodic images may be distinct faces on small meshes. */
static inline double material_canonical_delta(double xi,double xj,double image,double period){
 return (xj-xi)+period*nearbyint((image-xj)/period);
}
#ifdef GFS_MATERIAL_POSITIVE_FAN
#include "../gfs_material_positive_fan2d.h"
#endif
/* Single-rank bridge. Periodic evolution and reflected-wall diagnostics.
 * apply=0 exports a diagnostic; apply=1 writes material rates to host staging
 * arrays. The public RT wrapper rejects apply until gravity is connected. */
#include "../gfs_material_constraint_project2d.h"
#ifdef GFS_MATERIAL_DIRECT_PAIR
#include "../gfs_material_direct_pair2d.h"
#endif
#if defined(GFS_MATERIAL_STAR_ENTROPY) && !defined(GFS_MATERIAL_DIRECT_PAIR)
#error Limited material star requires its differentiated direct pair constraint
#endif
#include "../gfs_material_free_budget2d.h"
#ifdef GFS_MATERIAL_SWEPT_CFL
#include "gfs_material_stage_time.h"
#endif
#include "../gfs_material_limit_maps2d.h"
#ifdef GFS_MATERIAL_GRADIENT_PROJECT
#include "../gfs_material_project_shared2d.h"
#endif
#include "../gfs_material_rt_owner.h"
#include "../gfs_material_gravity2d.h"
#include "../gfs_material_point_gravity.h"
#if defined(GFS_MATERIAL_POINT_GRAVITY) && (!defined(GFS_MATERIAL_DIRECT_FV) || !defined(GFS_MATERIAL_DIRECT_WALLS))
#error "Point-gravity experiment requires direct FV and direct walls"
#endif
#include "../gfs_material_damping_bound2d.h"
#ifdef GFS_MATERIAL_CUDA_RHS
#include "../gfs_material_cuda_rhs.h"
#endif
extern const GfsMaterialRtOwner *gfs_material_rt_owner_map(int *count);
typedef struct {long long id;double x,y;} GfsMaterialSupportSite;
static GfsMaterialSupportSite *gfs_material_support_sites;
static int gfs_material_support_n;
static double gfs_material_support_period[2];
static int gfs_material_support_cmp(const void*a,const void*b){long long x=((const GfsMaterialSupportSite*)a)->id,y=((const GfsMaterialSupportSite*)b)->id;return (x>y)-(x<y);}
void gfs_material_step_support_end(void){free(gfs_material_support_sites);gfs_material_support_sites=NULL;gfs_material_support_n=0;}
int gfs_material_step_support_begin(int n,const long long*ids,const double*xy,double lx,double ly){
 gfs_material_step_support_end();if(n<1||!ids||!xy||!(lx>0)||ly<0)return 0;
 gfs_material_support_sites=malloc((size_t)n*sizeof(*gfs_material_support_sites));if(!gfs_material_support_sites)return 0;
 for(int i=0;i<n;i++)gfs_material_support_sites[i]=(GfsMaterialSupportSite){ids[i],xy[2*i],xy[2*i+1]};
 qsort(gfs_material_support_sites,n,sizeof(*gfs_material_support_sites),gfs_material_support_cmp);
 for(int i=1;i<n;i++)if(gfs_material_support_sites[i-1].id==gfs_material_support_sites[i].id){gfs_material_step_support_end();return 0;}
 gfs_material_support_n=n;gfs_material_support_period[0]=lx;gfs_material_support_period[1]=ly;return 1;
}
static GfsTraceDirected2D *gfs_material_support_edges(int n,int count,const GfsTraceDirected2D*edge,const ParticleSoA*p){
 if(n!=gfs_material_support_n)return NULL;
 double *xy=malloc((size_t)2*n*sizeof(double));GfsTraceDirected2D*out=malloc((size_t)count*sizeof(*out));
 if(!xy||!out){free(xy);free(out);return NULL;}
 for(int i=0;i<n;i++){GfsMaterialSupportSite key={(long long)p->indx[i],0,0};
  GfsMaterialSupportSite *v=bsearch(&key,gfs_material_support_sites,n,sizeof(key),gfs_material_support_cmp);
  if(!v){free(xy);free(out);return NULL;}xy[2*i]=v->x;xy[2*i+1]=v->y;
 }
 for(int k=0;k<count;k++){out[k]=edge[k];int i=edge[k].i,j=edge[k].j;
  double dx=xy[2*j]-xy[2*i],dy=xy[2*j+1]-xy[2*i+1];
  if(gfs_material_support_period[0]>0)dx+=gfs_material_support_period[0]*nearbyint((edge[k].dx-dx)/gfs_material_support_period[0]);
  if(gfs_material_support_period[1]>0)dy+=gfs_material_support_period[1]*nearbyint((edge[k].dy-dy)/gfs_material_support_period[1]);
  out[k].dx=dx;out[k].dy=dy;
 }free(xy);return out;
}
static double gfs_material_constraint_dt;
static GfsMaterialFreePair2D *gfs_material_free_pairs;
static GfsMaterialFreeSupport2D *gfs_material_free_support;
static int gfs_material_free_count,gfs_material_free_valid;
static double gfs_material_free_budget;
int gfs_material_free_rhs_valid_for_dt(double dt){
 int valid=0;if(gfs_material_mpi_rank()<=0)valid=gfs_material_free_valid&&dt<=gfs_material_free_budget&&gfs_material_free_pairs_support_valid2d(gfs_material_free_count,gfs_material_free_pairs,gfs_material_free_support,dt,.9);
 gfs_material_mpi_bcast(&valid,sizeof(valid));return valid;
}
static double gfs_material_last_face_crossing;
static double gfs_material_last_geometry_bound;
static double gfs_material_last_damping_bound;
static double gfs_material_face_cfl_bound(double courant){
 double dt=courant>0&&isfinite(courant)&&gfs_material_last_face_crossing>0?courant*gfs_material_last_face_crossing:0;
 /* SSPRK2 frozen negative-real interval; geometric CFL remains independent. */
 if(dt>0&&gfs_material_last_damping_bound>0)dt=fmin(dt,2/gfs_material_last_damping_bound);
 /* geometry_bound already includes its own Kepler sweep coefficient. */
 if(dt>0&&gfs_material_last_geometry_bound>0)dt=fmin(dt,gfs_material_last_geometry_bound);
 return dt;
}
void gfs_material_constraint_step_budget(double dt){gfs_material_constraint_dt=dt;gfs_material_free_valid=0;}
static double gfs_material_last_ledger[3];
static double gfs_material_last_potential[2];
static int gfs_material_last_potential_valid;
/* Valid only for the most recent successful applying RHS evaluation.
 * A geometry refresh or RK stage still needs its own matching evaluation. */
int gfs_material_potential_ledger(double *energy,double *power){
 if(!energy||!power||!gfs_material_last_potential_valid)return 0;
 *energy=gfs_material_last_potential[0];*power=gfs_material_last_potential[1];return 1;
}
void gfs_material_rhs_ledger(double out[3]){memcpy(out,gfs_material_last_ledger,sizeof(gfs_material_last_ledger));}
#include "gfs_material_resolved_faces_native.inc"
static int gfs_material_trace_probe_gravity_impl(GPUContext *ctx,int n,int nt,double lx,double ly,double gamma,int apply,const int *owner,const unsigned char *image_reflection,const double *potential_gradient){
 if(apply){gfs_material_last_potential_valid=0;gfs_material_last_face_crossing=0;gfs_material_last_damping_bound=0;gfs_material_last_geometry_bound=0;
  gfs_material_free_valid=0;gfs_material_free_count=0;free(gfs_material_free_pairs);gfs_material_free_pairs=NULL;free(gfs_material_free_support);gfs_material_free_support=NULL;
 }
 /* This adapter is periodic in x: a nonzero linear x-potential is not single-valued. */
 if(potential_gradient&&(!image_reflection||potential_gradient[0]!=0||!isfinite(potential_gradient[1])))return 0;
 int damping_ready=0;double damping_ms=0;
 const char *damping_cfl=getenv("GFS_MATERIAL_DAMPING_CFL");
 int damping_enabled=damping_cfl&&!strcmp(damping_cfl,"1");
#ifdef GFS_MATERIAL_DIRECT_FV
 int unsupported=damping_enabled;
#ifndef GFS_MATERIAL_DIRECT_WALLS
 unsupported|=image_reflection!=NULL;
#endif
#ifndef GFS_MATERIAL_POINT_GRAVITY
 unsupported|=potential_gradient!=NULL;
#else
 const char *point_constraints=getenv("GFS_MATERIAL_CONSTRAINTS");
 if(potential_gradient&&point_constraints&&!strcmp(point_constraints,"1"))unsupported=1;
#endif
 if(unsupported){fprintf(stderr,"[MATERIAL_DIRECT_UNSUPPORTED] gravity_wall_or_adjoint_damping=1\n");return 0;}
#endif
 double profile_start=MPI_Wtime();
 const int *alias=owner?owner:gpu_unified_owner_aliases(ctx,n,nt,lx,ly);
 ParticleSoA *p=&ctx->h_parts;FaceCSR *f=&ctx->h_faces;
 int capacity=f->face_offset[n],count=0;GfsTraceCsr2D t={0};
 unsigned char *walls=image_reflection?calloc((size_t)capacity,1):NULL;
 GfsTraceDirected2D *edge=malloc((size_t)capacity*sizeof(*edge));if(!edge||(image_reflection&&!walls)){free(edge);free(walls);return 0;}
 if(!gfs_material_resolved_faces2d(ctx,n,nt,lx,ly,alias,image_reflection,edge,walls,&count)){free(edge);free(walls);return 0;}
 GfsTraceDirected2D *shared_fit=gfs_trace_centroid_coordinates2d(n,count,edge,walls);
 int ok=shared_fit&&gfs_trace_csr2d_build_coordinates(n,count,edge,&t,walls,shared_fit);if(!ok){fprintf(stderr,"[MATERIAL_TRACE_FAIL] phase=csr cells=%d edges=%d\n",n,count);const char *dump=getenv("GFS_MATERIAL_TRACE_FAILURE_DUMP");
  if(dump&&*dump){FILE *fp=fopen(dump,"wbx");if(fp){
   int header[2]={n,count};fwrite(header,sizeof(int),2,fp);
   fwrite(edge,sizeof(*edge),(size_t)count,fp);fclose(fp);
   /* Opt-in failed-geometry evidence; never changes the state or gate. */
   size_t pathlen=strlen(dump)+16;char *rawpath=malloc(pathlen);
   if(rawpath){snprintf(rawpath,pathlen,"%s.raw.txt",dump);FILE *raw=fopen(rawpath,"wx");
    if(raw){fprintf(raw,"%d %d %d\n",n,nt,capacity);
     for(int v=0;v<nt;v++)fprintf(raw,"%d %d %.17g %.17g %.17g\n",v,alias[v],p->x[v],p->y[v],p->weight_W[alias[v]]);
     for(int v=0;v<n;v++)for(int h=f->face_offset[v];h<f->face_offset[v+1];h++)
      fprintf(raw,"%d %d %d %.17g %.17g %.17g %.17g\n",v,h,f->neighbor_idx[h],f->c1x[h],f->c1y[h],f->c2x[h],f->c2y[h]);
     fclose(raw);
    }free(rawpath);
   }
  }}
  free(shared_fit);free(edge);free(walls);return 0;}
 GfsFvPrimitive *state=calloc((size_t)n,sizeof(*state));
 double *rhs=calloc((size_t)n*(potential_gradient?15:6),sizeof(double));
 if(!state||!rhs){free(shared_fit);free(edge);free(walls);free(state);free(rhs);gfs_trace_csr2d_free(&t);return 0;}
 GfsMaterialFaceResponse *same_rhs_cache=NULL;
 double *q=rhs+4*n,*realized=q+n,*center=potential_gradient?realized+n:NULL,*gravity_rhs=potential_gradient?center+2*n:NULL;
 double potential_power=0,potential_energy=0;
 if(potential_gradient)ok=gfs_material_volume_centroid2d(n,count,edge,p->weight_volume64,center);
 if(!ok&&potential_gradient){
  fprintf(stderr,"[MATERIAL_TRACE_FAIL] phase=native_volume_centroid\n");
  const char *prefix=getenv("GFS_MATERIAL_FAILED_GEOMETRY");
  if(prefix&&*prefix){char path[2048];snprintf(path,sizeof(path),"%s.volumes.csv",prefix);FILE*fp=fopen(path,"wx");
   if(fp){fprintf(fp,"cell,nativeV,all_normalV,kept_normalV,all_crossV,omitted_normalV,omitted_length,omitted_faces\n");
    for(int i=0;i<n;i++){double all=0,kept=0,cross=0,omit=0,olen=0;int omitted=0;
     for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
      int j=f->neighbor_idx[k];if(j<0||j>=nt)continue;
      double dx=p->x[j]-p->x[i],dy=p->y[j]-p->y[i],d=hypot(dx,dy);
      double x1=f->c1x[k],x2=f->c2x[k],y1=f->c1y[k],y2=f->c2y[k],A=hypot(x2-x1,y2-y1);
      if(!(d>0))continue;
      double contribution=.25*A*((x1+x2)*dx+(y1+y2)*dy)/d;
      all+=contribution;cross+=.5*(x1*y2-x2*y1);
      if(A<1e-12){omit+=contribution;olen+=A;omitted++;}else kept+=contribution;
     }
     fprintf(fp,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",i,p->weight_volume64[i],all,kept,cross,omit,olen,omitted);
    }fclose(fp);
   }
  }
 }
 for(int i=0;i<n;i++){
  state[i].rho=p->mass[i]/p->weight_volume64[i];
  state[i].p=(gamma-1)*p->weight_U64[i]/p->weight_volume64[i];
  state[i].u[0]=p->vx[i];state[i].u[1]=p->vy[i];
 }
 int iterations=0;double relative=0;
 double profile_trace=MPI_Wtime();
 double *alpha3=NULL;
#ifndef GFS_MATERIAL_GRADIENT_PROJECT
 const char *limit_maps=getenv("GFS_MATERIAL_LIMIT_MAPS");
 if(ok&&limit_maps&&!strcmp(limit_maps,"1")){
  alpha3=malloc((size_t)n*3*sizeof(double));
  ok=alpha3&&gfs_material_limit_maps_boundaries2d(&t,count,edge,state,alpha3,walls);
  if(!ok)fprintf(stderr,"[MATERIAL_TRACE_FAIL] phase=limit_maps\n");
 }
#endif /* Projected thermodynamics replaces the scalar limiter. */
#ifdef GFS_MATERIAL_POSITIVE_FAN
 /* Periodic candidate only: wall images need an explicit geometric adapter. */
 if(ok&&image_reflection){fprintf(stderr,"[MATERIAL_POSITIVE_FAN] walls_not_qualified=1\n");ok=0;}
 if(ok)ok=gfs_material_positive_fan2d(n,count,edge,&t);
#endif
#ifdef GFS_MATERIAL_FACE_NORMAL_LIMIT
 if(ok){ok=gfs_material_face_normal_prepare2d(&t,state);if(!ok)fprintf(stderr,"[MATERIAL_TRACE_FAIL] phase=face_normal_prepare\n");}
#endif
#ifdef GFS_MATERIAL_GRADIENT_PROJECT
 if(ok){
  const char *con=getenv("GFS_MATERIAL_CONSTRAINTS"),*damp=getenv("GFS_MATERIAL_DAMPING_CFL");
  if((con&&strcmp(con,"0"))||(damp&&strcmp(damp,"0")))ok=0;
  if(ok)ok=gfs_material_project_shared2d(&t,count,edge,state,walls,shared_fit);
  if(!ok)fprintf(stderr,"[MATERIAL_TRACE_FAIL] phase=gradient_project\n");
 }
#endif
 free(shared_fit);shared_fit=NULL;
 int point_gravity=0;
#ifdef GFS_MATERIAL_POINT_GRAVITY
 point_gravity=potential_gradient!=NULL;
#endif
 if(ok&&potential_gradient&&!point_gravity){
  ok=gfs_material_gravity_rhs2d(&t,count,edge,walls,state,center,potential_gradient,gamma,alpha3,gravity_rhs,&potential_power,&iterations);
  if(ok)for(int i=0;i<n;i++){
   for(int a=0;a<4;a++)rhs[4*i+a]=gravity_rhs[7*i+a];q[i]=gravity_rhs[7*i+4];realized[i]=gravity_rhs[7*i+5];
   potential_energy+=p->mass[i]*(potential_gradient[0]*(p->x[i]+center[2*i])+potential_gradient[1]*(p->y[i]+center[2*i+1]));
  }
 }else if(ok){
#ifdef GFS_MATERIAL_CUDA_RHS
  #ifdef GFS_MATERIAL_CUDA_RESPONSE_EXPORT
  const char *export_response=getenv("GFS_MATERIAL_CUDA_RESPONSE_EXPORT");
  if(apply&&damping_enabled&&export_response&&!strcmp(export_response,"1"))same_rhs_cache=malloc((size_t)t.faces*sizeof(*same_rhs_cache));
  ok=gfs_material_cuda_rhs_cache2d(&t,state,gamma,alpha3,rhs,same_rhs_cache)==1;
#else
  ok=gfs_material_cuda_rhs2d(&t,state,gamma,alpha3,rhs)==1;
#endif
#else
  ok=gfs_material_rhs_maps2d(&t,state,gamma,alpha3,rhs)==1;
#endif
 }

/* Diagnostic-only frozen geometry/state; never changes acceptance. */
 if(!ok||getenv("GFS_MATERIAL_TRACE_STATE")){const char *prefix=getenv("GFS_MATERIAL_TRACE_STATE");if(!prefix)prefix=getenv("GFS_MATERIAL_FAILED_GEOMETRY");if(prefix&&*prefix){
  /* Opt-in short-run stage trace; the default still saves only the first
   * immutable state. A sequence is diagnostic, never a restart checkpoint. */
  char sequence_prefix[2048];const char *sequence=getenv("GFS_MATERIAL_TRACE_SEQUENCE");
  if(sequence&&!strcmp(sequence,"1")){
   static unsigned long sequence_number=0;
   snprintf(sequence_prefix,sizeof(sequence_prefix),"%s.rhs%06lu",prefix,++sequence_number);
   prefix=sequence_prefix;
  }
  char path[2048];snprintf(path,sizeof(path),"%s.edges.csv",prefix);FILE *fp=fopen(path,"wx");
  if(fp){fprintf(fp,"i,j,dx,dy,cx,cy,area\n");for(int k=0;k<count;k++)fprintf(fp,"%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",edge[k].i,edge[k].j,edge[k].dx,edge[k].dy,edge[k].cx,edge[k].cy,edge[k].area);fclose(fp);}
  snprintf(path,sizeof(path),"%s.cells.csv",prefix);fp=fopen(path,"wx");
  if(fp){fprintf(fp,"cell,x,y,V,rho,p,vx,vy\n");for(int i=0;i<n;i++)fprintf(fp,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",i,p->x[i],p->y[i],p->weight_volume64[i],state[i].rho,state[i].p,state[i].u[0],state[i].u[1]);fclose(fp);}
  if(sequence&&!strcmp(sequence,"1")){
   snprintf(path,sizeof(path),"%s.operator.csv",prefix);fp=fopen(path,"wx");
   if(fp){fprintf(fp,"cell,id,arho,ap,av,Fx,Fy,Vdot,H\n");for(int i=0;i<n;i++)fprintf(fp,"%d,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",i,(unsigned long long)p->indx[i],alpha3?alpha3[3*i]:1.,alpha3?alpha3[3*i+1]:1.,alpha3?alpha3[3*i+2]:1.,rhs[4*i],rhs[4*i+1],rhs[4*i+2],rhs[4*i+3]);fclose(fp);}
  }
 }}
 double profile_rhs=MPI_Wtime();
 if(ok&&(!potential_gradient||point_gravity))ok=gfs_material_weight_rate_boundaries2d(n,count,edge,state,rhs,q,realized,&iterations,&relative,walls);
 if(ok&&point_gravity){
  double centroid_potential=0;
  for(int i=0;i<n&&ok;i++){
   double x[3]={p->x[i],p->y[i],0},v[3]={state[i].u[0],state[i].u[1],0};
   double grad[3]={potential_gradient[0],potential_gradient[1],0};
   GfsMaterialPointGravity source;
   ok=gfs_material_point_gravity(2,p->mass[i],x,v,grad,&source);
   if(!ok)break;
   rhs[4*i]+=source.force[0];rhs[4*i+1]+=source.force[1];
   /* No change to internal-energy rate or to geometric volume/weight rates. */
   potential_energy+=source.potential;potential_power+=source.potential_rate;
   centroid_potential+=source.potential+p->mass[i]*(grad[0]*center[2*i]+grad[1]*center[2*i+1]);
  }
  if(getenv("GFS_MATERIAL_PROFILE"))fprintf(stderr,"[MATERIAL_POINT_POTENTIAL] generator=%.17g centroid=%.17g difference=%.17g\n",potential_energy,centroid_potential,centroid_potential-potential_energy);
 }
 double profile_solve=MPI_Wtime();
 const char *profile=getenv("GFS_MATERIAL_PROFILE");
 if(profile&&!strcmp(profile,"1"))fprintf(stderr,"[MATERIAL_PROFILE] cells=%d faces=%d entries=%d trace_ms=%.9g rhs_ms=%.9g solve_ms=%.9g iterations=%d ok=%d\n",n,t.faces,t.entries,1e3*(profile_trace-profile_start),1e3*(profile_rhs-profile_trace),1e3*(profile_solve-profile_rhs),iterations,ok);
 const char *constrained=getenv("GFS_MATERIAL_CONSTRAINTS");
 if(ok&&apply&&constrained&&!strcmp(constrained,"1")){
  GfsTraceDirected2D *support=gfs_material_support_n?gfs_material_support_edges(n,count,edge,p):NULL;
  if(gfs_material_support_n&&!support)ok=0;
  GfsMaterialActiveFace2D *active=malloc((size_t)count*2*sizeof(*active));
  GfsMaterialActiveFace2D *pending=malloc((size_t)count*sizeof(*pending));
  double *projected=malloc((size_t)n*7*sizeof(double));
  int selected=active?gfs_material_predict_constraints_support2d(n,count,edge,state,p->weight_W,q,gfs_material_constraint_dt,.9,count,active,support):-1;
  if(selected<0||!projected||!pending)ok=0;
#ifdef GFS_MATERIAL_DIRECT_FV
   /* The adjoint pressure projector is a DIFFERENT operator. Never apply
    * it to direct FV. Reject/retry until a direct constraint map is derived. */
   if(ok&&selected>0){
#ifdef GFS_MATERIAL_DIRECT_PAIR
    int its=0;
    for(int pairpass=0;ok&&pairpass<=count;pairpass++){
     ok=gfs_material_direct_pair2d(&t,count,edge,state,gamma,alpha3,selected,active,.9,p->weight_W,gfs_material_constraint_dt,support,projected,&its);
     if(!ok){fprintf(stderr,"[MATERIAL_DIRECT_PAIR] active=%d pass=%d solve_failed=1\n",selected,pairpass);break;}
     for(int i=0;i<n;i++){for(int a=0;a<4;a++)rhs[4*i+a]=projected[7*i+a];q[i]=projected[7*i+4];realized[i]=projected[7*i+5];}
     int remaining=gfs_material_predict_constraints_support2d(n,count,edge,state,p->weight_W,q,gfs_material_constraint_dt,.9,count,pending,support);
     fprintf(stderr,"[MATERIAL_DIRECT_PAIR] active=%d remaining=%d graph_iterations=%d pass=%d\n",selected,remaining,its,pairpass);
     if(remaining==0)break;
     if(remaining<0){ok=0;break;}
     int added=0;
     for(int k=0;k<remaining;k++){
      int found=0;for(int a=0;a<selected;a++)if(active[a].i==pending[k].i&&active[a].j==pending[k].j&&active[a].sign==pending[k].sign&&active[a].dx==pending[k].dx&&active[a].dy==pending[k].dy){found=1;break;}
      if(!found&&selected<2*count){active[selected++]=pending[k];added++;}
     }
     if(!added||pairpass==count){ok=0;break;}
    }
    /* Never fall through to the unrelated adjoint projector. */
    selected=0;
#else
    fprintf(stderr,"[MATERIAL_DIRECT_UNSUPPORTED] active_constraints=%d\n",selected);ok=0;
#endif
   }
#endif
#ifdef GFS_MATERIAL_DIRECT_PAIR
  /* A corrected direct RHS cannot use the uncorrected free-step certificate. */
  const char *free_reuse=NULL;
#else
  const char *free_reuse=getenv("GFS_MATERIAL_FREE_RHS_REUSE");
#endif
  if(ok&&selected==0&&!potential_gradient&&free_reuse&&!strcmp(free_reuse,"1")){
   gfs_material_free_pairs=malloc((size_t)count*sizeof(*gfs_material_free_pairs));
   if(support)gfs_material_free_support=malloc((size_t)count*sizeof(*gfs_material_free_support));
   if(gfs_material_free_pairs&&(!support||gfs_material_free_support)){
    for(int k=0;k<count;k++){
     int i=edge[k].i,j=edge[k].j;if(i>=j)continue;
     if(support)gfs_material_free_support[gfs_material_free_count]=(GfsMaterialFreeSupport2D){support[k].dx,support[k].dy};
     gfs_material_free_pairs[gfs_material_free_count++]=(GfsMaterialFreePair2D){edge[k].dx,edge[k].dy,state[j].u[0]-state[i].u[0],state[j].u[1]-state[i].u[1],p->weight_W[i]-p->weight_W[j],q[i]-q[j]};
    }
    gfs_material_free_budget=gfs_material_constraint_dt;gfs_material_free_valid=1;
   }
  }
  GfsMaterialConstraintReuse2D reuse={0};
  if(ok&&selected){
#ifdef GFS_MATERIAL_ACTIVE_RELEASE
   uint64_t *seen=calloc((size_t)count+1,sizeof(*seen));int seen_count=0;
   if(!seen)ok=0;
#endif
   for(int pass=0;ok&&pass<=count;pass++){
    if(getenv("GFS_MATERIAL_ACTIVE_IDS")){
     fprintf(stderr,"[MATERIAL_ACTIVE_IDS] dt=%.17g pass=%d rows=%d pairs=",gfs_material_constraint_dt,pass,selected);
     for(int k=0;k<selected;k++)fprintf(stderr,"%s%lld:%lld:%g",k?",":"",(long long)p->indx[active[k].i],(long long)p->indx[active[k].j],active[k].sign);
     fprintf(stderr,"\n");
    }
    int its=0;
#ifdef GFS_MATERIAL_ACTIVE_RELEASE
    /* Order-independent working-set fingerprint; repeats reject, never accept. */
    uint64_t fingerprint=(uint64_t)selected;
    for(int k=0;k<selected;k++){
     uint64_t x,y;memcpy(&x,&active[k].dx,sizeof x);memcpy(&y,&active[k].dy,sizeof y);
     uint64_t h=((uint64_t)active[k].i<<32)^(uint64_t)active[k].j^x^(y*UINT64_C(0x9e3779b97f4a7c15))^(active[k].sign>0?1:2);
     h^=h>>30;h*=UINT64_C(0xbf58476d1ce4e5b9);h^=h>>27;fingerprint+=h;
    }
    for(int k=0;k<seen_count;k++)if(seen[k]==fingerprint)ok=0;
    if(!ok){fprintf(stderr,"[MATERIAL_ACTIVE_RELEASE] cycle selected=%d pass=%d\n",selected,pass);break;}
    seen[seen_count++]=fingerprint;
    int release=-1;
    if(potential_gradient){
#ifdef GFS_MATERIAL_GRAVITY_REUSE
     ok=gfs_material_gravity_active_rhs2d_reuse_status(&t,count,edge,walls,state,center,potential_gradient,gamma,alpha3,projected,&potential_power,&its,selected,active,.9,&release,&reuse);
#else
     ok=gfs_material_gravity_active_rhs2d_status(&t,count,edge,walls,state,center,potential_gradient,gamma,alpha3,projected,&potential_power,&its,selected,active,.9,&release);
#endif
    }else{
     ok=gfs_material_project_active_reuse_drive_maps2d_status(&t,count,edge,state,gamma,selected,active,.9,projected,&its,&reuse,alpha3,walls,NULL,NULL,&release);
    }
    if(!ok&&release>=0&&release<selected&&pass<count){
     fprintf(stderr,"[MATERIAL_ACTIVE_RELEASE] remove=%d selected=%d pass=%d\n",release,selected,pass);
     memmove(active+release,active+release+1,(size_t)(selected-release-1)*sizeof(*active));selected--;ok=1;continue;
    }
#else
    if(potential_gradient)ok=gfs_material_gravity_active_rhs2d(&t,count,edge,walls,state,center,potential_gradient,gamma,alpha3,projected,&potential_power,&its,selected,active,.9);
    else ok=gfs_material_project_active_reuse_boundaries_maps2d(&t,count,edge,state,gamma,selected,active,.9,projected,&its,&reuse,alpha3,walls);
#endif
    if(!ok){fprintf(stderr,"[MATERIAL_TRACE_FAIL] phase=active_projection selected=%d iterations=%d pass=%d\n",selected,its,pass);break;}
    for(int i=0;i<n;i++){for(int j=0;j<4;j++)rhs[4*i+j]=projected[7*i+j];q[i]=projected[7*i+4];realized[i]=projected[7*i+5];}
    int remaining=gfs_material_predict_constraints_support2d(n,count,edge,state,p->weight_W,q,gfs_material_constraint_dt,.9,count,pending,support);
    fprintf(stderr,"[MATERIAL_CONSTRAINT] active=%d remaining=%d iterations=%d pass=%d dt=%.17g\n",selected,remaining,its,pass,gfs_material_constraint_dt);
    if(remaining==0)break;
    if(remaining<0){ok=0;break;}
    int added=0;
    for(int k=0;k<remaining;k++){
     int found=0;for(int j=0;j<selected;j++)if(active[j].i==pending[k].i&&active[j].j==pending[k].j&&active[j].sign==pending[k].sign&&active[j].dx==pending[k].dx&&active[j].dy==pending[k].dy){found=1;break;}
     if(!found&&selected<2*count){active[selected++]=pending[k];added++;}
    }
    if(!added||pass==count){ok=0;break;}
   }
#ifdef GFS_MATERIAL_ACTIVE_RELEASE
   free(seen);
#endif
  }
#ifdef GFS_MATERIAL_REUSE_RESPONSE
  /* Same immutable RHS geometry/state/trace coefficients: reuse the face
   * responses already certified by the active constraint solve. Do not
   * retain this cache across stages or retries. */
  if(ok&&apply&&damping_enabled&&reuse.frozen&&reuse.frozen->setup_ready&&
     reuse.frozen->n==n&&reuse.frozen->nf==t.faces){
   double tick=MPI_Wtime();
   ok=gfs_material_damping_bound2d(&t,reuse.frozen->face,alpha3,p->mass,&gfs_material_last_damping_bound);
   damping_ready=ok;damping_ms=1e3*(MPI_Wtime()-tick);
  }
#endif
  gfs_material_constraint_reuse_clear2d(&reuse);
  free(support);free(active);free(pending);free(projected);
 }

 if(profile&&!strcmp(profile,"1"))fprintf(stderr,"[MATERIAL_ACTIVE_PROFILE] elapsed_ms=%.9g total_adapter_ms=%.9g ok=%d\n",1e3*(MPI_Wtime()-profile_solve),1e3*(MPI_Wtime()-profile_start),ok);
 const char *export_state=getenv("GFS_MATERIAL_EXPORT_STATE");
 if(ok&&potential_gradient&&export_state&&!strcmp(export_state,"1")){
  FILE*dump=fopen("material_frozen_state.txt","wx");if(!dump)ok=0;
  if(dump){
   fprintf(dump,"%d %d %.17g %.17g %.17g\n",n,count,gamma,potential_gradient[0],potential_gradient[1]);
   for(int i=0;i<n;i++)fprintf(dump,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",p->x[i],p->y[i],p->mass[i],p->weight_volume64[i],state[i].rho,state[i].p,state[i].u[0],state[i].u[1],center[2*i],center[2*i+1],alpha3?alpha3[3*i]:1.,alpha3?alpha3[3*i+1]:1.,alpha3?alpha3[3*i+2]:1.);
   for(int k=0;k<count;k++)fprintf(dump,"%d %d %.17g %.17g %.17g %.17g %.17g %d\n",edge[k].i,edge[k].j,edge[k].dx,edge[k].dy,edge[k].cx,edge[k].cy,edge[k].area,walls?walls[k]:0);
   if(ferror(dump))ok=0;if(fclose(dump))ok=0;
  }
 }
 double energy=0,max_volume_error=0;
 if(apply)memset(gfs_material_last_ledger,0,sizeof(gfs_material_last_ledger));
 FILE *rates=ok&&!apply?fopen("material_rates_native.bin","wbx"):NULL;
 if(!apply&&!rates)ok=0;
 for(int i=0;i<n&&ok;i++){
  double udot=-state[i].p*realized[i]+rhs[4*i+3];
  double row[7]={rhs[4*i],rhs[4*i+1],rhs[4*i+2],rhs[4*i+3],q[i],realized[i],udot};
  energy+=state[i].u[0]*row[0]+state[i].u[1]*row[1]+udot;
  max_volume_error=fmax(max_volume_error,fabs(realized[i]-row[2]));
  if(apply){
   gfs_material_last_ledger[0]+=row[3];
   gfs_material_last_ledger[1]+=-state[i].p*realized[i];
   gfs_material_last_ledger[2]+=state[i].u[0]*row[0]+state[i].u[1]*row[1];
   p->ax_out[i]=row[0]/p->mass[i];p->ay_out[i]=row[1]/p->mass[i];
   p->weight_Udot64[i]=udot;p->die_out[i]=udot;p->weight_rate[i]=q[i];p->weight_force[i]=0;
  }else ok=fwrite(row,sizeof(row),1,rates)==1;
 }
 if(rates&&fclose(rates))ok=0;
 /* Optional immutable-state comparison, once per process; no state edits. */
 if(ok&&apply&&getenv("GFS_MATERIAL_AUDIT_FIRST_RHS")){
  static int saved_first_rhs=0;
  if(!saved_first_rhs){
   FILE*f=fopen("material_first_rhs.bin","wbx");
   double header[3]={(double)n,gfs_material_constraint_dt,gamma};
   if(!f)ok=0;
   if(f){
    if(fwrite(header,sizeof(header),1,f)!=1)ok=0;
    for(int i=0;i<n&&ok;i++){
     double row[18]={p->x[i],p->y[i],state[i].u[0],state[i].u[1],p->mass[i],p->weight_volume64[i],p->weight_W[i],p->weight_U64[i],state[i].p,state[i].rho,(double)p->indx[i],rhs[4*i],rhs[4*i+1],rhs[4*i+2],rhs[4*i+3],q[i],realized[i],-state[i].p*realized[i]+rhs[4*i+3]};
     if(fwrite(row,sizeof(row),1,f)!=1)ok=0;
    }
    if(fclose(f))ok=0;
   }
   if(ok)saved_first_rhs=1;
  }
 }


/* Diagnostic/opt-in face-crossing CFL for the current material rates.
  * Uses actual Laguerre boundary normal motion, including Wdot; not the
  * legacy inertial-weight oscillator frequency. No timestep accepted here. */
 if(ok&&apply){
  double *crossing=calloc((size_t)n,sizeof(double));
  if(!crossing)ok=0;
#ifdef GFS_MATERIAL_SWEPT_CFL
  const char*ge=getenv("GFS_VARPRESS_GEOM_CFL");char*end=NULL;
  double cg=ge?strtod(ge,&end):0,*sweep=NULL;
  if(!isfinite(cg)||cg<0||cg>1||(ge&&(!*ge||*end)))ok=0;
  if(ok&&cg>0){sweep=calloc((size_t)n,sizeof(double));if(!sweep)ok=0;}
#endif
  for(int k=0;k<count&&ok;k++){
   const GfsTraceDirected2D *e=edge+k;int i=e->i,j=e->j;
   double d=hypot(e->dx,e->dy),nx=e->dx/d,ny=e->dy/d;
   double vi=nx*state[i].u[0]+ny*state[i].u[1],vj=nx*state[j].u[0]+ny*state[j].u[1];
   double wn=.5*(vi+vj)+((e->cx-.5*e->dx)*(state[i].u[0]-state[j].u[0])+(e->cy-.5*e->dy)*(state[i].u[1]-state[j].u[1]))/d+.5*(q[i]-q[j])/d;
   if(walls&&walls[k]){wn=0;vj=-vi;}
   double ci=sqrt(gamma*state[i].p/state[i].rho),cj=sqrt(gamma*state[j].p/state[j].rho);
   double signal=fmax(ci+fabs(vi-wn),cj+fabs(vj-wn));
   if(!(signal>0)||!isfinite(signal)){ok=0;break;}crossing[i]+=e->area*signal;
#ifdef GFS_MATERIAL_SWEPT_CFL
   if(sweep){
    double absolute;
    if(walls&&walls[k])absolute=e->area*fabs(vi); /* Actual stationary plane. */
    else{
     double delta[2]={e->dx,e->dy},centroid[2]={e->cx,e->cy};
     if(!gfs_material_power_segment_sweep(delta,centroid,e->area,state[i].u,state[j].u,q[i],q[j],&absolute)){ok=0;break;}
    }
    sweep[i]+=absolute;if(!isfinite(sweep[i])){ok=0;break;}
   }
#endif
  }
  double dt=HUGE_VAL;
  for(int i=0;i<n&&ok;i++){if(!(crossing[i]>0)||!isfinite(crossing[i])){ok=0;break;}dt=fmin(dt,p->weight_volume64[i]/crossing[i]);}
  if(ok&&dt>0&&isfinite(dt))gfs_material_last_face_crossing=dt;else ok=0;
#ifdef GFS_MATERIAL_SWEPT_CFL
  if(ok&&sweep){
   double bound=HUGE_VAL;int limiting=-1;
   for(int i=0;i<n;i++)if(sweep[i]>0){
    double h=cg*p->weight_volume64[i]/sweep[i];
    if(!(h>0)||!isfinite(h)){ok=0;break;}
    if(h<bound){bound=h;limiting=i;}
   }
   double time;int phase;
   if(ok&&!gfs_material_native_evaluation_time(&time,&phase))ok=0;
   if(ok){gfs_material_last_geometry_bound=bound;
    fprintf(stderr,"[MATERIAL_GEOMETRIC_CFL] t=%.17g phase=%d coefficient=%.17g dt=%.17g id=%lld\n",
     time,phase,cg,bound,limiting>=0?(long long)p->indx[limiting]:-1LL);}
  }
  free(sweep);
#endif
  free(crossing);
 }
 if(ok&&apply&&damping_enabled){
  if(!damping_ready&&same_rhs_cache){
   double tick=MPI_Wtime();ok=gfs_material_damping_bound2d(&t,same_rhs_cache,alpha3,p->mass,&gfs_material_last_damping_bound);
   damping_ready=ok;damping_ms=1e3*(MPI_Wtime()-tick);
  }
  if(ok&&!damping_ready){
   double tick=MPI_Wtime();GfsMaterialFaceResponse *cache=malloc((size_t)t.faces*sizeof(*cache));
   ok=cache&&gfs_material_constraint_cache_maps2d(&t,state,gamma,alpha3,cache)&&
       gfs_material_damping_bound2d(&t,cache,alpha3,p->mass,&gfs_material_last_damping_bound);
   free(cache);damping_ms=1e3*(MPI_Wtime()-tick);
  }
  if(profile&&!strcmp(profile,"1"))fprintf(stderr,"[MATERIAL_DAMPING_CFL] lambda_bound=%.17g dt=%.17g ms=%.9g ok=%d reused=%d\n",gfs_material_last_damping_bound,gfs_material_last_damping_bound>0?2/gfs_material_last_damping_bound:0,damping_ms,ok,damping_ready);
 }


 if(ok&&apply&&!gfs_material_native_volume_sample(n,p->indx,p->weight_volume64,realized)){
  fprintf(stderr,"[MATERIAL_VOLUME_AUDIT] invalid periodic/wall stage sampling\n");ok=0;
 }
 if(!apply||!ok)fprintf(stderr,"[MATERIAL_RHS_PROBE] success=%d iterations=%d residual=%.17g volume_error=%.17g energy_rate=%.17g\n",ok,iterations,relative,max_volume_error,energy);
 if(ok&&apply&&potential_gradient){
  gfs_material_last_potential[0]=potential_energy;
  gfs_material_last_potential[1]=potential_power;
  gfs_material_last_potential_valid=1;
 }
 if(potential_gradient)fprintf(stderr,"[MATERIAL_GRAVITY_POWER] success=%d potential=%.17g potential_rate=%.17g total_rate=%.17g\n",ok,potential_energy,potential_power,energy+potential_power);
 free(same_rhs_cache);free(alpha3);free(shared_fit);free(edge);free(walls);free(state);free(rhs);
 if(!ok||apply){if(!ok)gfs_material_free_valid=0;gfs_trace_csr2d_free(&t);return ok;}
 uint64_t h[4]={UINT64_C(0x4641434541535331),(uint64_t)n,(uint64_t)t.faces,(uint64_t)t.entries};
 FILE *out=fopen("material_trace_native.bin","wbx");if(!out){gfs_trace_csr2d_free(&t);return 0;}
 ok=fwrite(h,sizeof(h),1,out)==1;
 for(int i=0;i<n&&ok;i++){
  double V=p->weight_volume64[i],rho=p->mass[i]/V,P=(gamma-1)*p->weight_U64[i]/V;
  double row[5]={rho,P,p->vx[i],p->vy[i],V};
  if(!(V>0)||!(P>0)||!(rho>0))ok=0;
  else ok=fwrite(row,sizeof(row),1,out)==1;
 }
 if(ok)ok=fwrite(t.face,sizeof(*t.face),t.faces,out)==(size_t)t.faces&&fwrite(t.entry,sizeof(*t.entry),t.entries,out)==(size_t)t.entries;
 if(fclose(out))ok=0;
 fprintf(stderr,"[MATERIAL_TRACE_PROBE] cells=%d directed=%d faces=%d entries=%d success=%d\n",n,count,t.faces,t.entries,ok);
 gfs_trace_csr2d_free(&t);return ok;
}

static int gfs_material_trace_probe_impl(GPUContext *ctx,int n,int nt,double lx,double ly,
 double gamma,int apply,const int *owner,const unsigned char *image_reflection){
 return gfs_material_trace_probe_gravity_impl(ctx,n,nt,lx,ly,gamma,apply,owner,image_reflection,NULL);
}

/* RT is diagnostic-only until compatible gravity and its time ledger are
 * connected. Never run the periodic alias resolver on mirror sentinel IDs. */
static int gfs_material_trace_probe(GPUContext *ctx,int n,int nt,double lx,double ly,
 double gamma,int apply,double ymin,double ymax,double acceleration_x,double acceleration_y){
 if(apply)gfs_material_last_potential_valid=0;
 const char *ic=getenv("EUNHA_IC");
 if(!ic||strcmp(ic,"rt"))return gfs_material_trace_probe_gravity_impl(ctx,n,nt,lx,ly,gamma,apply,NULL,NULL,NULL);
 const char *native_rt=getenv("GFS_MATERIAL_RT_NATIVE");
 if(apply&&(!native_rt||strcmp(native_rt,"1"))){fprintf(stderr,"[MATERIAL_RT] native qualification opt-in required\n");return 0;}
 int count=0;const GfsMaterialRtOwner *pad=gfs_material_rt_owner_map(&count);
 int *alias=malloc((size_t)nt*sizeof(*alias));unsigned char *reflection=malloc((size_t)nt);
 ParticleSoA *p=&ctx->h_parts;
 int ok=alias&&reflection&&gfs_material_rt_resolve_owners(n,nt,p->indx,p->x,p->y,count,pad,lx,ymin,ymax,alias,reflection);
 double grad[2]={-acceleration_x,-acceleration_y};
 if(ok)ok=gfs_material_trace_probe_gravity_impl(ctx,n,nt,lx,ly,gamma,apply,alias,reflection,grad);
 free(alias);free(reflection);return ok;
}

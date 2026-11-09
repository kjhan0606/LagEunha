#ifdef GFS_MATERIAL_FIXED_SOLID
#include "../gfs_fixed_cylinder_contact.h"
#define CYL_MATERIAL_STREAM_MAGIC "GFSFIXS1"
static int cyl_solid_read(FILE*f){uint64_t h;return fread(&h,sizeof(h),1,f)==1&&h==gfs_fixed_cylinder_signature();}
static int cyl_solid_write(FILE*f){uint64_t h=gfs_fixed_cylinder_signature();return fwrite(&h,sizeof(h),1,f)==1;}
#else
#define CYL_MATERIAL_STREAM_MAGIC "GFSSTRM3"
static int cyl_solid_read(FILE*f){return 1;}
static int cyl_solid_write(FILE*f){return 1;}
#endif
/* RunCylinder's opt-in material boundary. The pressure/work operator stays
 * in NativeU64; this file owns only stage samples, accepted transport and IO.
 * Initial conditions/restarts must use the explicit native checkpoint codec.
 */
#include <float.h>
#include <limits.h>
#include "../gfs_material_open_stream.h"
#include "../gfs_native_checkpoint_io.h"
#include "../gfs_material_radial_image.h"
#include "../NativeU64/gfs_material_stage_time.h"
#include "material_boundary.h"
static GfsMaterialOpenStream cyl_material_stream={.dim=2};
static GfsCylinderMaterialSites cyl_material_sites;
static int64_t cyl_material_cursor;
static double cyl_material_h,cyl_material_W,cyl_material_rho,cyl_material_u,cyl_material_p;
static double cyl_material_exchange[4],cyl_material_work,cyl_material_wall,cyl_material_repair;
#ifdef GFS_MATERIAL_SEARCH_EXTENT
#include "../gfs_power_search_extent.h"
static int cyl_material_search_valid;
static double cyl_material_search_lower[2],cyl_material_search_upper[2];
int gfs_material_cylinder_search_extent(SimParameters*s,double*xmin,double*ymin,
 double*xmax,double*ymax,const double*deficit){
 if(!gfs_material_cylinder_mode()||NID(s)!=1)return -1;
 double lo[2]={*xmin,*ymin},hi[2]={*xmax,*ymax},cell=BASICCELL_CELLWIDTH(s);
 if(cyl_material_search_valid)for(int d=0;d<2;d++){
  lo[d]=fmin(lo[d],cyl_material_search_lower[d]);hi[d]=fmax(hi[d],cyl_material_search_upper[d]);
 }
 double zero[4]={0,0,0,0};
 if(gfs_power_expand_search_box(2,cell,zero,lo,hi)<0)return -1;
 int changed=0;
 if(deficit){
  changed=gfs_power_expand_search_box(2,cell,deficit,lo,hi);
  if(changed<0)return -1;
  if(changed){
   for(int d=0;d<2;d++){cyl_material_search_lower[d]=lo[d];cyl_material_search_upper[d]=hi[d];}
   cyl_material_search_valid=1;
   fprintf(stderr,"[CYLINDER_SEARCH_EXTEND] xmin=%.17g ymin=%.17g xmax=%.17g ymax=%.17g bin=%.17g\n",lo[0],lo[1],hi[0],hi[1],cell);
  }
 }
 *xmin=lo[0];*ymin=lo[1];*xmax=hi[0];*ymax=hi[1];return changed;
}
#endif
static void cyl_material_fail(SimParameters*s,const char*why){
 fprintf(stderr,"[CYLINDER_MATERIAL] %s\n",why);MPI_Abort(MPI_COMM(s),99);
}
static void cyl_material_sites_clear(void){
 free(cyl_material_sites.alias);free(cyl_material_sites.image);
 free(cyl_material_sites.owned_ids);free(cyl_material_sites.lowx);free(cyl_material_sites.lowy);
 memset(&cyl_material_sites,0,sizeof(cyl_material_sites));
}
const GfsCylinderMaterialSites*gfs_material_cylinder_sites(int n,int nt,const long long*ids){
 if(!ids||cyl_material_sites.n!=n||cyl_material_sites.nt!=nt)return NULL;
 for(int i=0;i<n;i++)if(ids[i]!=cyl_material_sites.owned_ids[i])return NULL;
 return &cyl_material_sites;
}
static void cyl_material_set(treevorostressrk4particletype*p,double x,double y,
 double vx,double vy,double m,double U,double V,double W,unsigned long long id,double gamma){
 memset(p,0,sizeof(*p));p->x=x;p->y=y;p->vx=vx;p->vy=vy;p->mass=m;p->u4if.indx=id;
 p->weight=(GfsWeightState){W,0,1,GFS_WEIGHT_STATE_TAG};
 p->weight_thermal=(GfsWeightThermal){U,0,V};p->volume=V;p->ie=U;
 p->den=m/V;p->pressure=(gamma-1)*U/V;p->csound=sqrt(gamma*p->pressure/p->den);
 p->stress.K=p->pressure/pow(p->den,gamma);
 p->stress.E_inv_xx=U;p->stress.E_inv_xy=V;p->stress.E_inv_yy=2026102603.;
 /* w2 is a legacy field. The material geometry uses weight.W exclusively. */
 p->w2=p->w2old=0;
}
void gfs_material_cylinder_padding(SimParameters*s,postype unused){
 (void)unused;double time;int phase,n=VORO_NP(s),ny=NY(s);
 if(!gfs_material_cylinder_mode()||NID(s)!=1||
    TVORORK4_DDINFO(s)[0].n_size!=sizeof(treevorostressrk4particletype)||
    !gfs_material_native_evaluation_time(&time,&phase))cyl_material_fail(s,"missing stage/payload context");
 double h=cyl_material_h,lx=SIMBOX(s).x.max,ly=SIMBOX(s).y.max,gamma=GAS_GAMMA(s);
 double center[2]={CYL_CX(s),CYL_CY(s)},R=CYL_R(s),minx=HUGE_VAL;
 treevorostressrk4particletype*real=(void*)VORORK4_TBP(s);
 size_t capacity=(size_t)3*n+(size_t)24*(ny+8)+3*cyl_material_stream.count+32;
 if(capacity>INT_MAX)cyl_material_fail(s,"padding capacity overflow");
 treevorostressrk4particletype*pad=my_malloc(capacity*sizeof(*pad));
 cyl_material_sites_clear();GfsCylinderMaterialSites*b=&cyl_material_sites;
 b->n=n;b->center[0]=center[0];b->center[1]=center[1];b->radius=R;b->lx=lx;b->ly=ly;
 b->alias=malloc((n+capacity)*sizeof(int));b->image=malloc((n+capacity)*sizeof(int));
 b->owned_ids=malloc((size_t)n*sizeof(long long));b->lowx=calloc(n+capacity,sizeof(double));b->lowy=calloc(n+capacity,sizeof(double));
 if(!pad||!b->alias||!b->image||!b->owned_ids||!b->lowx||!b->lowy)cyl_material_fail(s,"padding allocation");
 for(int i=0;i<n;i++){b->alias[i]=i;b->image[i]=-1;b->owned_ids[i]=PINDX(real+i);minx=fmin(minx,real[i].x);}
 int k=0;
 #define CYL_PAD_NEW() do{if((size_t)k>=capacity)cyl_material_fail(s,"padding count overflow");b->alias[n+k]=b->image[n+k]=-1;}while(0)
 for(int i=0;i<n;i++){
  treevorostressrk4particletype*a=real+i;
  for(int shift=-1;shift<=1;shift+=2)if((shift==1&&a->y<4*h)||(shift==-1&&a->y>ly-4*h)){
   CYL_PAD_NEW();pad[k]=*a;GfsPowerDD y=gfs_pdd_sum(a->y,shift*ly);
   pad[k].y=y.hi;b->lowy[n+k]=y.lo;b->alias[n+k]=i;k++;
  }
  double radial_distance=hypot(a->x-center[0],a->y-center[1]);
  /* Supply the complete near-side family, independently of self-face size.
   * A vanishing self face can still support a neighbour's finite wall face.
   * At 2R the image reaches the centre; the all-cell centre-site certificate
   * in the native geometry path checks that this cutoff is inactive. */
#ifdef GFS_MATERIAL_FIXED_SOLID
  int fixed_face;double fixed_normal[2];
  if(!isfinite(a->x)||!isfinite(a->y)||!(gfs_fixed_cylinder_outside_box(a->x,a->y)||gfs_fixed_cylinder_gap(a->x,a->y,&fixed_face,fixed_normal)>=0))
   cyl_material_fail(s,"invalid generator at fixed cylinder");
#else
  if(!isfinite(radial_distance)||!(radial_distance>R))
   cyl_material_fail(s,"invalid generator at cylinder");
#endif
  if(radial_distance<2*R
#ifdef GFS_MATERIAL_FIXED_SOLID
      && 0
#endif
  ){
   CYL_PAD_NEW();pad[k]=*a;double pos[2]={a->x,a->y},v[2]={a->vx,a->vy},high[2],low[2],u[2];GfsMaterialRadialImage image;
   if(!gfs_material_radial_image(2,center,R,pos,&image)||
      !gfs_material_radial_split_point(2,center,R,pos,high,low)||
      !gfs_material_radial_vector(2,&image,v,0,u))cyl_material_fail(s,"radial image");
   pad[k].x=high[0];pad[k].y=high[1];pad[k].vx=u[0];pad[k].vy=u[1];
   pad[k].u4if.indx=MAX_INDEX;SET_FLAG(pad+k,BoundaryGhostflag);
   b->lowx[n+k]=low[0];b->lowy[n+k]=low[1];b->image[n+k]=i;k++;
  }
 }
 int lo=(int)floor((-4*h-cyl_material_u*time)/h)-1,hi=(int)ceil(-cyl_material_u*time/h)+1;
 for(int col=lo;col<=hi;col++){
  double x=(col+.5)*h+cyl_material_u*time;
  if(x>=minx-.5*h||x< -4*h)continue;
  for(int j=-4;j<ny+4;j++){
   CYL_PAD_NEW();cyl_material_set(pad+k,x,(j+.5)*h,cyl_material_u,0,cyl_material_rho*h*h,
     cyl_material_p*h*h/(gamma-1),h*h,cyl_material_W,MAX_INDEX,gamma);
   SET_FLAG(pad+k,BoundaryGhostflag);k++;
  }
 }
 double period[3]={0,ly,0};
 for(size_t j=0;j<cyl_material_stream.count;j++){
  const GfsMaterialOpenSample*a=cyl_material_stream.sample+j;double pos[3];
  if(!gfs_material_open_position(2,a,time,period,pos)||pos[0]<lx||a->v[0]<=0)cyl_material_fail(s,"exterior reentry");
  for(int shift=-1;shift<=1;shift++){
   double y=pos[1]+shift*ly;if(y< -4*h||y>ly+4*h)continue;
   CYL_PAD_NEW();cyl_material_set(pad+k,pos[0],y,a->v[0],a->v[1],a->mass,a->U,a->V,a->W,MAX_INDEX,gamma);
   SET_FLAG(pad+k,BoundaryGhostflag);k++;
  }
 }
 #undef CYL_PAD_NEW
 b->nt=n+k;VORORK4_TBPP(s)=(void*)pad;VORO_NPAD(s)=k;
}
int gfs_material_cylinder_search(SimParameters*s,double end,double*xmin,double*ymin,double*xmax,double*ymax){
 if(!gfs_material_cylinder_mode()||!isfinite(end)||!(cyl_material_h>0))return 0;
 double h=cyl_material_h,lx=SIMBOX(s).x.max,ly=SIMBOX(s).y.max,period[3]={0,ly,0};
 *xmin=*ymin=-8*h;*xmax=lx+8*h;*ymax=ly+8*h;
 for(size_t k=0;k<cyl_material_stream.count;k++){
  double x[3];if(!gfs_material_open_position(2,cyl_material_stream.sample+k,end,period,x))return 0;
  *xmax=fmax(*xmax,x[0]+8*h);
 }
 BASICCELL_CELLWIDTH(s)=4*h;
#ifdef GFS_MATERIAL_SEARCH_EXTENT
 if(gfs_material_cylinder_search_extent(s,xmin,ymin,xmax,ymax,NULL)<0)return 0;
#endif
 return 1;
}
void gfs_material_cylinder_accept(double work,double wall,double repair){
 cyl_material_work+=work;cyl_material_wall+=wall;cyl_material_repair+=repair;
}
static double cyl_material_energy(SimParameters*s){
 treevorostressrk4particletype*p=(void*)VORORK4_TBP(s);long double E=0;
 for(int i=0;i<VORO_NP(s);i++)E+=p[i].weight_thermal.U+.5L*p[i].mass*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
 return E;
}
static void cyl_material_transport(const treevorostressrk4particletype*p,int sign){
 cyl_material_exchange[0]+=sign*p->mass;cyl_material_exchange[1]+=sign*p->mass*p->vx;
 cyl_material_exchange[2]+=sign*p->mass*p->vy;
 cyl_material_exchange[3]+=sign*(p->weight_thermal.U+.5*p->mass*(p->vx*p->vx+p->vy*p->vy));
}
static int cyl_material_events(SimParameters*s,double t){
 treevorostressrk4particletype*p=(void*)VORORK4_TBP(s);
 double lx=SIMBOX(s).x.max,h=cyl_material_h,gamma=GAS_GAMMA(s);int old=VORO_NP(s),k=0,removed=0,added=0;
 int64_t columns=(int64_t)floor((cyl_material_u*t+32*DBL_EPSILON*fmax(1.,cyl_material_u*t))/h+.5);
 if(columns<cyl_material_cursor)cyl_material_fail(s,"inflow cursor reversed");
 for(int i=0;i<old;i++){
  if(p[i].x>=lx){
   GfsMaterialOpenSample a={(int64_t)PINDX(p+i),t,{p[i].x,p[i].y,0},{p[i].vx,p[i].vy,0},p[i].mass,p[i].weight_thermal.U,p[i].weight_thermal.V,p[i].weight.W};
   if(a.v[0]<=0||!gfs_material_open_push(&cyl_material_stream,&a))cyl_material_fail(s,"departure state");
   cyl_material_transport(p+i,-1);removed++;
  }else {if(p[i].x<0)cyl_material_fail(s,"unexpected upstream departure");if(k!=i)p[k]=p[i];k++;}
 }
 size_t new_capacity=(size_t)k+(size_t)(columns-cyl_material_cursor)*NY(s);
 if(new_capacity>INT_MAX)cyl_material_fail(s,"owned count overflow");
 if(new_capacity>(size_t)old){
  treevorostressrk4particletype*next=my_malloc(new_capacity*sizeof(*next));
  if(!next)cyl_material_fail(s,"inflow allocation");memcpy(next,p,(size_t)k*sizeof(*p));my_free(p);p=next;
 }
 for(int64_t c=cyl_material_cursor;c<columns;c++)for(int j=0;j<NY(s);j++){
  double x=cyl_material_u*t-(c+.5)*h;if(fabs(x)<64*DBL_EPSILON)x=0;
  cyl_material_set(p+k,x,(j+.5)*h,cyl_material_u,0,cyl_material_rho*h*h,
    cyl_material_p*h*h/(gamma-1),h*h,cyl_material_W,100000000ULL+c*NY(s)+j,gamma);
  cyl_material_transport(p+k,1);added++;k++;
 }
 cyl_material_cursor=columns;VORORK4_TBP(s)=(void*)p;VORORK4_BP(s)=(void*)p;VORO_NP(s)=VORO_TNP(s)=k;
 if(added||removed)fprintf(stderr,"[CYLINDER_MATERIAL_EVENT] t=%.17g before=%d after=%d added=%d removed=%d mass_exchange=%.17g energy_exchange=%.17g\n",t,old,k,added,removed,cyl_material_exchange[0],cyl_material_exchange[3]);
 cyl_material_sites_clear();return added||removed;
}
static int cyl_material_checksum(const char*path,uint64_t*out){
 FILE*f=fopen(path,"rb");if(!f)return 0;uint64_t h=UINT64_C(14695981039346656037);unsigned char b[4096];size_t n;
 while((n=fread(b,1,sizeof(b),f)))for(size_t i=0;i<n;i++){h^=b[i];h*=UINT64_C(1099511628211);}
 int ok=!ferror(f);if(fclose(f))ok=0;*out=h;return ok;
}
/* Native checkpoint contract: geometry and inflow must not change silently
 * when a stage sequence resumes. No parameter defaults are inferred from W. */
static void cyl_material_configuration(SimParameters*s,int32_t dims[2],double v[9]){
 dims[0]=NX(s);dims[1]=NY(s);
 double a[9]={GAS_GAMMA(s),CYL_RHO(s),CYL_UINF(s),CYL_P(s),SIMBOX(s).x.max,
  SIMBOX(s).y.max,CYL_CX(s),CYL_CY(s),CYL_R(s)};memcpy(v,a,sizeof(a));
}
static void cyl_material_read_boundary(SimParameters*s,int step,double t){
 char owned[256],path[280];snprintf(owned,sizeof(owned),"cylmat.%06d.dat",step);
 snprintf(path,sizeof(path),"%s.exterior",owned);
 if(!gfs_material_open_read(path,t,&cyl_material_stream)||cyl_material_stream.dim!=2)cyl_material_fail(s,"missing exterior checkpoint");
 uint64_t actual,oh,eh;double saved,config[9],expected[9];int32_t dims[2],want[2];char magic[8];
 cyl_material_configuration(s,want,expected);snprintf(path,sizeof(path),"%s.stream",owned);FILE*f=fopen(path,"rb");
 if(!f||fread(magic,1,8,f)!=8||memcmp(magic,CYL_MATERIAL_STREAM_MAGIC,8)||fread(&saved,8,1,f)!=1||
    fread(&cyl_material_W,8,1,f)!=1||fread(&cyl_material_cursor,8,1,f)!=1||
    fread(&oh,8,1,f)!=1||fread(&eh,8,1,f)!=1||fread(dims,sizeof(dims),1,f)!=1||fread(config,sizeof(config),1,f)!=1||memcmp(dims,want,sizeof(dims))||memcmp(config,expected,sizeof(config))||!cyl_solid_read(f)||fgetc(f)!=EOF||ferror(f)||saved!=t||
    !isfinite(cyl_material_W)||cyl_material_cursor<0)cyl_material_fail(s,"stream checkpoint");
 fclose(f);if(!cyl_material_checksum(owned,&actual)||actual!=oh)cyl_material_fail(s,"owned checkpoint mismatch");
 snprintf(path,sizeof(path),"%s.exterior",owned);if(!cyl_material_checksum(path,&actual)||actual!=eh)cyl_material_fail(s,"exterior checkpoint mismatch");
}
void gfs_material_cylinder_checkpoint_boundary(SimParameters*s,const char*owned,double t){
 if(!gfs_material_cylinder_mode()||!owned||!isfinite(t))cyl_material_fail(s,"invalid boundary checkpoint context");
 int32_t dims[2];double config[9];cyl_material_configuration(s,dims,config);
 char path[280],temp[300];uint64_t oh,eh;
 if(snprintf(path,sizeof(path),"%s.exterior",owned)>=(int)sizeof(path))cyl_material_fail(s,"checkpoint path too long");
 if(!gfs_material_open_write(path,t,&cyl_material_stream)||!cyl_material_checksum(owned,&oh)||!cyl_material_checksum(path,&eh))cyl_material_fail(s,"checkpoint write");
 snprintf(path,sizeof(path),"%s.stream",owned);snprintf(temp,sizeof(temp),"%s.tmp",path);FILE*f=fopen(temp,"wb");
 if(!f||fwrite(CYL_MATERIAL_STREAM_MAGIC,1,8,f)!=8||fwrite(&t,8,1,f)!=1||fwrite(&cyl_material_W,8,1,f)!=1||fwrite(&cyl_material_cursor,8,1,f)!=1||fwrite(&oh,8,1,f)!=1||fwrite(&eh,8,1,f)!=1||fwrite(dims,sizeof(dims),1,f)!=1||fwrite(config,sizeof(config),1,f)!=1||!cyl_solid_write(f)||fclose(f)||rename(temp,path))cyl_material_fail(s,"checkpoint manifest write");
}
static void cyl_material_write(SimParameters*s,int step,double t,double dt){
 gfs_native_checkpoint_outdata(s,step,t,dt,"cylmat");
 char owned[256];snprintf(owned,sizeof(owned),"cylmat.%06d.dat",step);
 gfs_material_cylinder_checkpoint_boundary(s,owned,t);
}
extern void gfs_rk4_set_step_context(int,double,double);
static int cyl_material_run(SimParameters*s,int step){
#ifdef GFS_MATERIAL_SEARCH_EXTENT
 cyl_material_search_valid=0; /* computational cache is not a physical checkpoint variable */
#endif
 if(NID(s)!=1||step<1||GAS_AVMODE(s)<1||GAS_EVOLMETHOD(s)!=3||!GAS_GPU_ENABLED(s)||
    NX(s)<16||NY(s)<4||!(CYL_UINF(s)>0)||GAS_Kappa(s)!=0)
  cyl_material_fail(s,"requires one-rank GPU material mode and explicit versioned cylmat restart");
#ifdef GFS_MATERIAL_FIXED_SOLID
 if(CYL_CX(s)!=2.5||CYL_CY(s)!=2.||CYL_R(s)!=.5)cyl_material_fail(s,"fixed solid coordinates mismatch");
#endif
 double lx=SIMBOX(s).x.max,ly=SIMBOX(s).y.max;
 Xmin_HydroExam(s)=Ymin_HydroExam(s)=0;Xmax_HydroExam(s)=lx;Ymax_HydroExam(s)=ly;
 cyl_material_h=lx/NX(s);cyl_material_rho=CYL_RHO(s);cyl_material_u=CYL_UINF(s);cyl_material_p=CYL_P(s);
 if(fabs(cyl_material_h-ly/NY(s))>32*DBL_EPSILON*cyl_material_h)cyl_material_fail(s,"square initial spacing required");
 postype t,dt;gfs_native_checkpoint_readdata(s,&t,&dt,step,"cylmat");
 if(Xmax_HydroExam(s)!=lx||Ymax_HydroExam(s)!=ly)cyl_material_fail(s,"checkpoint domain mismatch");
 cyl_material_read_boundary(s,step,t);VORORK4_TBPP(s)=NULL;VORO_NPAD(s)=0;
 double initial_energy=cyl_material_energy(s),start=t;
 const char*stop=getenv("HYDRO_TSTOP"),*cap_env=getenv("GFS_MATERIAL_DT");
 double end=stop?strtod(stop,NULL):AMAX(s),cap=cap_env?strtod(cap_env,NULL):0;
 if(!(end>t)||!(cap>0)||!isfinite(end)||!isfinite(cap))cyl_material_fail(s,"invalid end/dt");
 double interval=ASTEP(s)>0?ASTEP(s):.1,next_output=t+interval;
 if(!gfs_material_cylinder_refresh(s,t))cyl_material_fail(s,"restart geometry");
#ifdef GFS_MATERIAL_FIXED_SOLID
 if(getenv("GFS_MATERIAL_FIXED_IC")){
  if(t!=0||cyl_material_cursor!=0)cyl_material_fail(s,"fixed-solid IC requires t=0");
  treevorostressrk4particletype*ic=(void*)VORORK4_TBP(s);
  for(int i=0;i<VORO_NP(s);i++){
   double V=ic[i].weight_thermal.V;
   if(!(V>0)||!isfinite(V))cyl_material_fail(s,"fixed IC volume");
   cyl_material_set(ic+i,ic[i].x,ic[i].y,ic[i].vx,ic[i].vy,cyl_material_rho*V,cyl_material_p*V/(GAS_GAMMA(s)-1),V,ic[i].weight.W,PINDX(ic+i),GAS_GAMMA(s));
  }
  cyl_material_write(s,step+1,t,0);
  fprintf(stderr,"[CYLINDER_FIXED_IC_COMPLETE] cells=%td solid=%llu\n",VORO_NP(s),(unsigned long long)gfs_fixed_cylinder_signature());
  return 1;
 }
#endif
 while(t<end){
  double next_entry=(cyl_material_cursor+.5)*cyl_material_h/cyl_material_u;
  double ceiling=fmin(cap,fmin(end-t,next_entry-t));
  if(!(ceiling>0)||t+ceiling==t)cyl_material_fail(s,"nonadvancing ownership schedule");
  gfs_rk4_set_step_context(step,t,ceiling);
  dt=exam2d_vph_kdk_int_blend(s,gfs_material_cylinder_padding,cyl_w2Measure2D,
     searchCellRk4Neighbors2D,findCellRk4BP2D,cyl_evolBP,mkLinkedList2D_oExam,cyl_postStage);
  if(!(dt>0)||dt>ceiling*(1+8*DBL_EPSILON))cyl_material_fail(s,"invalid accepted duration");
  t+=dt;step++;
  if(cyl_material_events(s,t)&&!gfs_material_cylinder_refresh(s,t))cyl_material_fail(s,"handoff geometry");
  double e=cyl_material_energy(s)-initial_energy-cyl_material_exchange[3]-cyl_material_work;
  fprintf(stderr,"[CYLINDER_MATERIAL_ACCEPT] t=%.17g dt=%.17g cells=%td energy_error=%.17g image_wall_work=%.17g contact_ke=%.17g ledger_residual=%.17g\n",t,(double)dt,VORO_NP(s),e,cyl_material_wall,cyl_material_repair,e-cyl_material_repair);
  if(fabs(e-cyl_material_repair)>1e-9*fmax(1.,fabs(initial_energy)))cyl_material_fail(s,"accepted energy ledger");
  if(t>=next_output||t>=end){cyl_material_write(s,step,t,dt);next_output=t+interval;}
 }
 fprintf(stderr,"[CYLINDER_MATERIAL_COMPLETE] start=%.17g end=%.17g cells=%td physical_relative_energy_error=%.17g\n",start,(double)t,VORO_NP(s),(cyl_material_energy(s)-initial_energy-cyl_material_exchange[3]-cyl_material_work+cyl_material_wall)/initial_energy);
 gfs_material_open_free(&cyl_material_stream);cyl_material_sites_clear();return 1;
}

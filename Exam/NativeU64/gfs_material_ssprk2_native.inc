#include "../gfs_material_ssprk2.h"
/* Optional saved-state wall diagnostic. Never changes geometry or state. */
static void gfs_material_contact_cell_diagnostic(SimParameters*s,long long id,
 int attempt,int stage,const char*phase,double h){
 treevorostressrk4particletype*p=(void*)VORORK4_TBP(s);
 for(ptrdiff_t i=0;i<VORO_NP(s);i++)if((long long)PINDX(p+i)==id){
  fprintf(stderr,"[MATERIAL_CONTACT_CELL] t=%.17g attempt=%d stage=%d phase=%s dt=%.17g id=%lld x=%.17g y=%.17g vx=%.17g vy=%.17g mass=%.17g U=%.17g V=%.17g W=%.17g Fx=%.17g Fy=%.17g Udot=%.17g Wdot=%.17g\n",
   gfs_rk4_context_time,attempt,stage,phase,h,id,(double)p[i].x,(double)p[i].y,
   (double)p[i].vx,(double)p[i].vy,(double)p[i].mass,p[i].weight_thermal.U,
   p[i].weight_thermal.V,p[i].weight.W,p[i].mass*p[i].ax,p[i].mass*p[i].ay,
   p[i].weight_thermal.Udot,p[i].weight_rk4.k4W);
  return;
 }
 fprintf(stderr,"[MATERIAL_CONTACT_CELL] missing_id=%lld attempt=%d stage=%d phase=%s\n",id,attempt,stage,phase);
}
/* Optional pair diagnostic only; never edits state or acceptance. */
static void gfs_material_pair_diagnostic(SimParameters*s,int attempt,int stage,
 const char*phase,double h,const double period[3]){
 const char*env=getenv("GFS_MATERIAL_PAIR_DIAGNOSTIC");if(!env)return;
 long long id0,id1;if(sscanf(env,"%lld,%lld",&id0,&id1)!=2)return;
 treevorostressrk4particletype*p=(void*)VORORK4_TBP(s),*a=NULL,*b=NULL;
 for(ptrdiff_t i=0;i<VORO_NP(s);i++){if((long long)PINDX(p+i)==id0)a=p+i;if((long long)PINDX(p+i)==id1)b=p+i;}
 if(!a||!b)return;
 double dx=b->x-a->x,dy=b->y-a->y;
 if(period[0]>0)dx-=period[0]*nearbyint(dx/period[0]);
 if(period[1]>0)dy-=period[1]*nearbyint(dy/period[1]);
 double dw=a->weight.W-b->weight.W,d2=dx*dx+dy*dy;
 fprintf(stderr,"[MATERIAL_PAIR_DIAG] attempt=%d stage=%d phase=%s h=%.17g dx=%.17g dy=%.17g dvx=%.17g dvy=%.17g dw=%.17g dq=%.17g slack=%.17g ratio=%.17g\n",attempt,stage,phase,h,dx,dy,b->vx-a->vx,b->vy-a->vy,dw,a->weight_rk4.k4W-b->weight_rk4.k4W,.9*d2-fabs(dw),fabs(dw)/d2);
}

#ifdef GFS_MATERIAL_VOLUME_PRIMARY
/* Rebuild only geometry during inverse iterations. Successful final padding
 * and device faces must remain alive for the one ordinary RHS evaluation. */
static int gfs_material_ssprk2_geometry(SimParameters*s,
 double xmin,double ymin,double xmax,double ymax,double gamma,
 void (*padding)(SimParameters*,postype),Voro2D_point*(*neighbors)(SimParameters*,int,int,int*),
 treevorork4particletype*(*cells)(SimParameters*,int,int,int*),
 void (*linked)(SimParameters*,postype,postype,postype,postype,postype,void(*)(SimParameters*,postype)),
 int target_stage,double h){
 for(int iteration=0;;iteration++){
  updateDenW2Pressure2DBlend(s,xmin,ymin,xmax,ymax,gamma,padding,neighbors,cells,linked,1);
  if(gpu_geometry_trial_failed())return 0;
  if(!target_stage)return 1;
  int status=gfs_material_native_volume_correct(s,target_stage,h,iteration);
  if(status==1)return 1;
  if(status!=2){
   fprintf(stderr,"[MATERIAL_VOLUME_COORDINATE_REJECT] target_stage=%d iteration=%d\n",target_stage,iteration);
   return 0;
  }
  /* updateDen already freed BASICCELL. Only its retained padding is owned
   * here; the next call reconstructs aliases/images with the new weights. */
  if(VORORK4_TBPP(s)){my_free(VORORK4_TBPP(s));VORORK4_TBPP(s)=NULL;}
  VORO_NPAD(s)=0;gpu_get_context()->tess_faces_on_device=0;
 }
}
#endif

/* Global material stage operator, using existing geometry/rollback.
 * The MPI adapter supplies a consistent global view and handles ownership. */
static double gfs_material_ssprk2_native(SimParameters *s,
 void (*padding)(SimParameters *,postype),Voro2D_point *(*neighbors)(SimParameters *,int,int,int *),
 treevorork4particletype *(*cells)(SimParameters *,int,int,int *),
 void (*linked)(SimParameters *,postype,postype,postype,postype,postype,void (*)(SimParameters *,postype)),
 void (*post)(SimParameters *)){
 const char *ic=getenv("EUNHA_IC"),*hs=getenv("GFS_MATERIAL_DT");
 const char *energy_cell=getenv("GFS_MATERIAL_ENERGY_CELL");double energy_tx=0,energy_ty=0;
 const int energy_selected=energy_cell&&sscanf(energy_cell,"%lf,%lf",&energy_tx,&energy_ty)==2;
 double cap=hs?strtod(hs,NULL):0;
 const char *rt_gate=getenv("GFS_MATERIAL_RT_NATIVE");
 int material_rt=ic&&!strcmp(ic,"rt")&&rt_gate&&!strcmp(rt_gate,"1");
 int material_noh=gfs_material_noh_native_mode();
 int material_cyl=gfs_material_cylinder_mode();
 const char *contact_env=getenv("GFS_MATERIAL_CONTACT_AUDIT_ID");long long contact_id=0;
 int contact_audit=material_cyl&&contact_env&&sscanf(contact_env,"%lld",&contact_id)==1;
 const char *volume_audit_env=getenv("GFS_MATERIAL_VOLUME_AUDIT");
 int volume_audit=volume_audit_env&&!strcmp(volume_audit_env,"1");
#ifdef GFS_MATERIAL_VOLUME_PRIMARY
 const char*volume_coordinate=getenv("GFS_MATERIAL_VOLUME_COORDINATE");
 int volume_mode=0;
 if(volume_coordinate&&strcmp(volume_coordinate,"0")){
  if(!strcmp(volume_coordinate,"endpoint"))volume_mode=1;
  else if(!strcmp(volume_coordinate,"stages"))volume_mode=2;
  else{fprintf(stderr,"[MATERIAL_VOLUME_COORDINATE] unknown mode\n");MPI_Abort(MPI_COMM(s),99);}
  volume_audit=1;
 }
#endif
 gfs_material_native_volume_audit_end();
 const char *ek=getenv("GFS_MATERIAL_EXTERNAL_KICK");
 int external=ek&&!strcmp(ek,"1");
 const char *predict_env=getenv("GFS_MATERIAL_CFL_PREDICT");
 int predict_cfl=predict_env&&!strcmp(predict_env,"1")&&!external;
 const char *consistent_env=getenv("GFS_MATERIAL_CONSISTENT_STAGE_DT");
 int consistent_dt=consistent_env&&!strcmp(consistent_env,"1");
 const char *warm_env=getenv("GFS_MATERIAL_CFL_WARM_START");
 int warm_start=consistent_dt&&!external&&warm_env&&!strcmp(warm_env,"1");
 /* Only the preceding accepted step supplies a CFL slope. This is a
  * proposal, never an acceptance rule; stage-two checks and rollback remain.
  * Exact time continuity prevents reuse after a restart/time reset. */
 static double prior_cfl_end_time=-HUGE_VAL,prior_cfl_slope=0,older_cfl_slope=0,prior_cfl_dt=0,older_cfl_dt=0,prior_cfl_bound=0;
 double predicted_slope=predict_cfl&&prior_cfl_end_time==gfs_rk4_context_time?prior_cfl_slope:0;
 double gravity[3]={GAS_ACCX(s),GAS_ACCY(s),0};
#ifdef GFS_MATERIAL_TEST_GRAVITY_Y
 /* Isolated qualification build; production takes acceleration from s. */
 gravity[1]=GFS_MATERIAL_TEST_GRAVITY_Y;
#endif
 if(NID(s)!=1||post!=(material_cyl?cyl_postStage:periodic_postStage_blend)||GAS_EVOLMETHOD(s)!=3||!gfs_native_periodic_weight_kdk_mode()||
    !(cap>0)||!isfinite(cap)||!ic||(strcmp(ic,"gresho")&&strcmp(ic,"kh")&&strcmp(ic,"sod_periodic")&&strcmp(ic,"sedov2d")&&strcmp(ic,"shear_periodic")&&!material_rt&&!material_noh&&!material_cyl)||(!external&&!material_rt&&(GAS_ACCX(s)!=0||GAS_ACCY(s)!=0))||
    ((material_noh||material_cyl)&&(external||!GAS_GPU_ENABLED(s)||GAS_Kappa(s)!=0||!getenv("GFS_MATERIAL_FACE_CFL")||strcmp(getenv("GFS_MATERIAL_FACE_CFL"),"1")))||
    (material_rt&&(external||GAS_ACCX(s)!=0||Xmin_HydroExam(s)!=0||Ymin_HydroExam(s)!=0))||
    GAS_VISCOSITY(s)!=0||GAS_AlphaVis(s)!=0||GAS_BetaVis(s)!=0){
  fprintf(stderr,"[MATERIAL_SSP] unsupported qualification settings\n");MPI_Abort(MPI_COMM(s),99);
 }
 ptrdiff_t n=VORO_NP(s),total=VORO_TNP(s);size_t bytes=n*sizeof(treevorostressrk4particletype);
 treevorostressrk4particletype *original=malloc(bytes),*base=malloc(bytes);GfsWeightBaseIndex *map=malloc(n*sizeof(*map));
 double *first_force=malloc((size_t)n*2*sizeof(double));
 double *contact_snapshot=contact_audit?malloc((size_t)n*4*sizeof(double)):NULL;
 if(contact_audit&&!contact_snapshot)MPI_Abort(MPI_COMM(s),99);
 if(!original||!base||!map||!first_force)MPI_Abort(MPI_COMM(s),99);memcpy(base,VORORK4_TBP(s),bytes);memcpy(original,base,bytes);
 for(ptrdiff_t i=0;i<n;i++){map[i]=(GfsWeightBaseIndex){PINDX(base+i),i};if(base[i].weight.pi!=0)MPI_Abort(MPI_COMM(s),99);}
 qsort(map,n,sizeof(*map),gfs_weight_base_order);
 for(ptrdiff_t i=1;i<n;i++)if(map[i-1].id==map[i].id)MPI_Abort(MPI_COMM(s),99);
 if(material_noh)for(ptrdiff_t i=0;i<n;i++)if(map[i].id!=(unsigned long long)i){
  fprintf(stderr,"[NOH_NATIVE_BOUNDARY] owned IDs must preserve the IC grid labels\n");MPI_Abort(MPI_COMM(s),99);
 }
 /* Noh's search domain below spans 2.25 L. Use N spatial-index bins across
  * it, rather than inheriting the periodic 4 L/N bin width. This changes
  * neither HydroGridSize nor the physical exterior samples; the GPU's
  * complete adaptive power-cell search still certifies all candidates. */
 if(material_noh)BASICCELL_CELLWIDTH(s)=2.25*Xmax_HydroExam(s)/NX(s);
#ifdef GFS_MATERIAL_SEARCH_EXTENT
 /* Scope runtime capacity to this driver/domain. A fresh process or a time
  * rewind resets it; SSP rollback inside this call keeps valid capacity. */
 if(!gfs_material_search_begin(s,(material_cyl||material_noh)?0:(material_rt?2:1),gfs_rk4_context_time))MPI_Abort(MPI_COMM(s),99);
#endif
 double gamma=GAS_GAMMA(s),cell=BASICCELL_CELLWIDTH(s),halo=gfs_kepler_halo_width(cell);
 double xmin=Xmin_HydroExam(s)-halo,ymin=Ymin_HydroExam(s)-halo,xmax=Xmax_HydroExam(s)+halo,ymax=Ymax_HydroExam(s)+halo;
 double period[3]={Xmax_HydroExam(s)-Xmin_HydroExam(s),Ymax_HydroExam(s)-Ymin_HydroExam(s),0};
 if(material_rt)period[1]=0; /* Reflecting y walls are not periodic images. */
 if(material_noh){
  /* Same finite exterior sample domain as the boundary qualification driver.
   * Search bins cover it; physical material cells are never wrapped. */
  double L=Xmax_HydroExam(s);xmin=ymin=-.625*L;xmax=ymax=1.625*L;
  period[0]=period[1]=0;linked=mkLinkedList2D_material_noh;
 }
 if(material_cyl){
  period[0]=0;
  if(!gfs_material_cylinder_search(s,gfs_rk4_context_time+cap,&xmin,&ymin,&xmax,&ymax))MPI_Abort(MPI_COMM(s),99);
 }
 const char *support_env=getenv("GFS_MATERIAL_STEP_SUPPORT");
 if(support_env&&!strcmp(support_env,"1")){
  long long*ids=malloc((size_t)n*sizeof(*ids));double*xy=malloc((size_t)2*n*sizeof(double));
  if(!ids||!xy)MPI_Abort(MPI_COMM(s),99);
  for(ptrdiff_t i=0;i<n;i++){ids[i]=(long long)PINDX(original+i);xy[2*i]=original[i].x;xy[2*i+1]=original[i].y;}
  int good=gfs_material_step_support_begin(n,ids,xy,period[0],period[1]);free(ids);free(xy);if(!good)MPI_Abort(MPI_COMM(s),99);
 }else gfs_material_step_support_end();
 float dtold=GAS_dtold(s);double floors=gfs_stage_floor_l;int nfloors=gfs_stage_floor_n_l;long long pairs=gfs_pair_face_ends;
 double ceiling=fmin(cap,gfs_rk4_step_ceiling>0?gfs_rk4_step_ceiling:cap),h=0,begin=MPI_Wtime();
 if(warm_start&&prior_cfl_end_time==gfs_rk4_context_time&&prior_cfl_bound>0&&isfinite(prior_cfl_bound)){
  ceiling=fmin(ceiling,prior_cfl_bound);
  fprintf(stderr,"[MATERIAL_CFL_WARM_START] proposal=%.17g prior_bound=%.17g\n",ceiling,prior_cfl_bound);
 }
 /* Cylinder f(t,y) is independent of the trial duration. Keep one successful
  * step-start host evaluation, not GPU geometry or rejected predictor states.
  * This cache exists only inside this call; every accepted step frees it. */
 treevorostressrk4particletype *first_rhs=NULL;
 double first_bound=0,first_ledger[3]={0},first_boundary=0,first_wall=0;
 double first_box[4]={0},first_floors=0;int first_nfloors=0,first_ready=0;
 long long first_pairs=0;
 int used_stage_cfl_retry=0,prediction_ready=0;
 int attempt_limit=12;
#ifdef GFS_MATERIAL_TEST_STOP_TIME
 /* Fault-injection build only: exercise complete accepted-state recovery
  * after real ownership events without changing any accepted update. */
 if(gfs_rk4_context_time>=GFS_MATERIAL_TEST_STOP_TIME)attempt_limit=0;
#endif
 for(int attempt=0;attempt<attempt_limit;attempt++){
  double current_box[4]={xmin,ymin,xmax,ymax};
#ifdef GFS_MATERIAL_SEARCH_EXTENT
  if(material_cyl&&gfs_material_cylinder_search_extent(s,&current_box[0],&current_box[1],
     &current_box[2],&current_box[3],NULL)<0)MPI_Abort(MPI_COMM(s),99);
#endif
  if(first_ready&&memcmp(current_box,first_box,sizeof(first_box))){
   first_ready=0;
   fprintf(stderr,"[MATERIAL_FIRST_RHS_INVALIDATE] reason=search_extent attempt=%d\n",attempt);
  }
  if(volume_audit){
   int ok=first_ready?gfs_material_native_volume_audit_rewind_first((int)n):
                      gfs_material_native_volume_audit_begin((int)n);
   if(!ok)MPI_Abort(MPI_COMM(s),99);
  }
  h=0; /* Stage-one failure must halve this trial ceiling, not a stale h. */
  double ledger[2][3]={{0}},mixing=0,thermal_time_increment=0,retry_bound=0,stage_bound[2]={0};
  double potential_initial=0,potential_final=0,potential_rates[2]={0,0};
  double boundary_power[2]={0,0},wall_power[2]={0,0},wall_repair[2]={0,0};
  int reject=0,retry_from_stage2=0;double external_work=0;gpu_geometry_trial_begin();
  if(external){
   memcpy(base,original,bytes);
   treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
   for(ptrdiff_t i=0;i<n;i++){
    ptrdiff_t j=gfs_weight_base_find(map,n,PINDX(p+i));if(j<0){reject=1;break;}
    GfsMaterialStepState a={p[i].mass,{p[i].x,p[i].y,0},{p[i].vx,p[i].vy,0},p[i].weight_thermal.U,p[i].weight.W},q;double work;
    if(!gfs_material_external_kick(2,.5*ceiling,gravity,&a,&q,&work)){reject=1;break;}
    p[i].vx=q.v[0];p[i].vy=q.v[1];base[j].vx=q.v[0];base[j].vy=q.v[1];external_work+=work;
   }
  }
  for(int stage=1;stage<=2&&!reject;stage++){
   if(!gfs_material_evaluation_set(gfs_rk4_context_time,stage==1?ceiling:h,stage)){reject=1;break;}
#ifdef GFS_MATERIAL_TEST_REJECT_ONCE
   /* Isolated test build only: reject after the predictor changed state. */
   static int material_test_rejected=0;
   if(stage==2&&!material_test_rejected){material_test_rejected=1;reject=1;break;}
#endif
   gfs_material_constraint_step_budget(stage==1?ceiling:h);
   double bound;
   if(stage==1&&first_ready){
    /* Rollback already discarded device caches. Restore the complete host
     * payload/order produced by the first successful evaluation. Reapply the
     * predictor and contact at the NEW h; stage2 rebuilds all device geometry. */
    if(VORO_NP(s)!=n||VORO_TNP(s)!=total||VORORK4_TBPP(s)||gpu_get_context()->initialized)
     MPI_Abort(MPI_COMM(s),99);
    memcpy(VORORK4_TBP(s),first_rhs,bytes);
    gfs_stage_floor_l=first_floors;gfs_stage_floor_n_l=first_nfloors;gfs_pair_face_ends=first_pairs;
    bound=first_bound;memcpy(ledger[0],first_ledger,sizeof(first_ledger));
    boundary_power[0]=first_boundary;wall_power[0]=first_wall;
    gfs_material_evaluation_clear();
    fprintf(stderr,"[MATERIAL_FIRST_RHS_REUSE] t=%.17g attempt=%d bound=%.17g\n",gfs_rk4_context_time,attempt,bound);
   }else{
    double material_stage_profile_start=MPI_Wtime();
#ifdef GFS_MATERIAL_VOLUME_PRIMARY
    if(!gfs_material_ssprk2_geometry(s,xmin,ymin,xmax,ymax,gamma,padding,neighbors,cells,linked,
       volume_mode==2&&stage==2?1:0,h)){reject=1;break;}
#else
    updateDenW2Pressure2DBlend(s,xmin,ymin,xmax,ymax,gamma,padding,neighbors,cells,linked,1);
    if(gpu_geometry_trial_failed()){reject=1;break;}
#endif
    double material_stage_profile_geometry=MPI_Wtime();
#ifdef GFS_MATERIAL_SEARCH_EXTENT
    /* Geometry may have expanded its value-argument box. Carry that same
     * box into force reconstruction, including its host fallback path. */
    if(gfs_material_search_extent(s,&xmin,&ymin,&xmax,&ymax,NULL)<0){reject=1;break;}
#endif
    bound=getAccVoro2DBlend_GPU(s,xmin,ymin,xmax,ymax,VoroAccuracyOrder(s),GAS_COURANT(s),gamma,padding,neighbors,cells,linked);
    gfs_material_evaluation_clear();
    if(getenv("GFS_MATERIAL_PROFILE"))fprintf(stderr,"[MATERIAL_STAGE_PROFILE] stage=%d density_geometry_ms=%.9g acceleration_ms=%.9g\n",stage,1e3*(material_stage_profile_geometry-material_stage_profile_start),1e3*(MPI_Wtime()-material_stage_profile_geometry));
    if(!(bound>0)||!isfinite(bound)){reject=1;break;}
    gfs_material_rhs_ledger(ledger[stage-1]);
    if((material_noh||material_cyl)&&!gfs_material_native_boundary_power(&boundary_power[stage-1])){reject=1;break;}
    if(material_cyl&&!gfs_material_native_wall_power(&wall_power[stage-1])){reject=1;break;}
    if(stage==1&&material_cyl&&!external){
     if(!first_rhs)first_rhs=malloc(bytes);
     if(!first_rhs)MPI_Abort(MPI_COMM(s),99);
     memcpy(first_rhs,VORORK4_TBP(s),bytes);
     first_bound=bound;memcpy(first_ledger,ledger[0],sizeof(first_ledger));
     first_boundary=boundary_power[0];first_wall=wall_power[0];
     first_floors=gfs_stage_floor_l;first_nfloors=gfs_stage_floor_n_l;first_pairs=gfs_pair_face_ends;
     double box[4]={xmin,ymin,xmax,ymax};
#ifdef GFS_MATERIAL_SEARCH_EXTENT
     if(gfs_material_cylinder_search_extent(s,&box[0],&box[1],&box[2],&box[3],NULL)<0)
      MPI_Abort(MPI_COMM(s),99);
#endif
     memcpy(first_box,box,sizeof(first_box));first_ready=1;
    }
   }
   stage_bound[stage-1]=bound;
   if(material_rt){double potential;
    if(!gfs_material_potential_ledger(&potential,&potential_rates[stage-1])){reject=1;break;}
    if(stage==1)potential_initial=potential;
   }
   if(stage==1){h=fmin(ceiling,bound);
    /* Slopes describe step-midpoint derivatives. Extrapolate to the
     * next midpoint on nonuniform time intervals, with no safety tuning.
     * If history is absent, retain the preceding slope alone. */
    if(!prediction_ready&&predict_cfl&&prior_cfl_end_time==gfs_rk4_context_time&&older_cfl_dt>0)
     predicted_slope=prior_cfl_slope+(prior_cfl_slope-older_cfl_slope)*(prior_cfl_dt+h)/(prior_cfl_dt+older_cfl_dt);
    /* Freeze this history-only proposal during rollback retries. Re-solving
     * its midpoint extrapolation after every shrink creates an unnecessary
     * fixed-point iteration. Actual stage CFL bounds are still reevaluated. */
    prediction_ready=1;
    if(predicted_slope<0&&isfinite(predicted_slope)){
     double proposal=bound/(1-predicted_slope);
     if(proposal>0&&isfinite(proposal))h=fmin(h,proposal);
    }
    /* Kick duration must equal the accepted hydro duration; discard the
     * entire trial if its first CFL estimate shortens the interval. */
    if((external||consistent_dt)&&h<ceiling){
     /* Cylinder evaluates f(t,y) without a trial-duration budget.
       * State, phase-one ledgers and volume rates are unchanged here.
       * Reuse only this just-computed RHS; stage-two/endpoint checks remain. */
     if(!external&&(material_cyl||gfs_material_free_rhs_valid_for_dt(h))){
      fprintf(stderr,"[MATERIAL_FREE_RHS_REUSE] evaluated=%.17g accepted_proposal=%.17g\n",ceiling,h);
     }else{
      if(consistent_dt)fprintf(stderr,"[MATERIAL_STAGE_BUDGET_RETRY] evaluated=%.17g proposed=%.17g\n",ceiling,h);
      retry_bound=h;h=ceiling;reject=1;break;
     }
    }
   }else if(h>bound*(1+2e-7)){retry_bound=bound;retry_from_stage2=1;reject=1;break;}
   gfs_material_pair_diagnostic(s,attempt,stage,"before",h,period);
   if(contact_audit)gfs_material_contact_cell_diagnostic(s,contact_id,attempt,stage,"rhs",h);
   treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
   for(ptrdiff_t i=0;i<n;i++){
    ptrdiff_t j=gfs_weight_base_find(map,n,PINDX(p+i));if(j<0){reject=1;break;}
    GfsMaterialStepState b={base[j].mass,{base[j].x,base[j].y,0},{base[j].vx,base[j].vy,0},base[j].weight_thermal.U,base[j].weight.W};
    GfsMaterialStepState a={p[i].mass,{p[i].x,p[i].y,0},{p[i].vx,p[i].vy,0},p[i].weight_thermal.U,p[i].weight.W},out;
    GfsMaterialStepRate r={{p[i].mass*p[i].ax,p[i].mass*p[i].ay,0},p[i].weight_thermal.Udot,p[i].weight_rk4.k4W};
    if(stage==1){first_force[2*j]=r.force[0];first_force[2*j+1]=r.force[1];}
    else{
     double dx=r.force[0]-first_force[2*j],dy=r.force[1]-first_force[2*j+1];
     double mix=h*h*(dx*dx+dy*dy)/(8*p[i].mass);mixing+=mix;
     if(material_cyl){
      /* Convex combination of kinetic energies, including the actual
       * predictor wall impulse. Diagnostic only; never changes state. */
      double dvx=a.v[0]-b.v[0]+h*r.force[0]/a.mass;
      double dvy=a.v[1]-b.v[1]+h*r.force[1]/a.mass;
      thermal_time_increment+=.125*a.mass*(dvx*dvx+dvy*dvy)-.25*h*h/a.mass*
       (first_force[2*j]*first_force[2*j]+first_force[2*j+1]*first_force[2*j+1]+
        r.force[0]*r.force[0]+r.force[1]*r.force[1]);
     }else thermal_time_increment-=mix;
    }
    if(!gfs_material_ssprk2_stage(2,stage,h,period,&b,&a,&r,&out)){reject=1;break;}
    /* Diagnostic only: isolate conservative Euler-stage kinetic subtraction
     * for one spatially selected material cell; no state/gate changes. */
    if(energy_selected&&stage==1){
     if(fabs(b.x[0]-energy_tx)<.001&&fabs(b.x[1]-energy_ty)<.001){
      double loss=h*h*(r.force[0]*r.force[0]+r.force[1]*r.force[1])/(2*a.mass);
      fprintf(stderr,"[MATERIAL_CELL_ENERGY] t=%.17g attempt=%d id=%lld h=%.17g U=%.17g Udot=%.17g linear_increment=%.17g kinetic_subtraction=%.17g trialU=%.17g V=%.17g vx=%.17g vy=%.17g Fx=%.17g Fy=%.17g\n",gfs_rk4_context_time,attempt,(long long)PINDX(p+i),h,a.U,r.Udot,h*r.Udot,loss,out.U,p[i].weight_thermal.V,a.v[0],a.v[1],r.force[0],r.force[1]);
     }
    }

    p[i].x=out.x[0];p[i].y=out.x[1];p[i].vx=out.v[0];p[i].vy=out.v[1];p[i].weight.W=out.W;p[i].weight.pi=0;
    p[i].weight_thermal.U=out.U;p[i].ie=(float)out.U;
    double V=p[i].weight_thermal.V,rho=p[i].mass/V,P=(gamma-1)*out.U/V;
    p[i].den=rho;p[i].pressure=P;p[i].csound=sqrt(gamma*P/rho);p[i].stress.K=P/pow(rho,gamma);
   }
   gfs_material_pair_diagnostic(s,attempt,stage,"after",h,period);
   if(reject)break;
   if(material_noh){
    /* One rank: keep immutable material labels and physical positions. No
     * periodic poststage or reservoir particles in the migration list. */
   }else if(material_rt){
    /* Reject crossing trials; no wrapping, clamping or wall-energy repair. */
    for(ptrdiff_t i=0;i<n;i++)if(!(p[i].y>Ymin_HydroExam(s)&&p[i].y<Ymax_HydroExam(s))){reject=1;break;}
    if(reject)break;
    migrateTreeVorork4Particles(s);
   }else if(material_cyl){
    /* Preserve the original obstacle repair. Its KE change is explicit,
     * never added to U; SSPRK2 weights the predictor repair by one half. */
    long double before=0,after=0;
    for(ptrdiff_t i=0;i<n;i++)before+=.5L*p[i].mass*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
    if(contact_audit){
     gfs_material_contact_cell_diagnostic(s,contact_id,attempt,stage,"pre_contact",h);
     for(ptrdiff_t i=0;i<n;i++){
      ptrdiff_t j=gfs_weight_base_find(map,n,PINDX(p+i));if(j<0)MPI_Abort(MPI_COMM(s),99);
      contact_snapshot[4*j]=p[i].x;contact_snapshot[4*j+1]=p[i].y;
      contact_snapshot[4*j+2]=p[i].vx;contact_snapshot[4*j+3]=p[i].vy;
     }
    }
    post(s);p=(void*)VORORK4_TBP(s);
    if(contact_audit){
     gfs_material_contact_cell_diagnostic(s,contact_id,attempt,stage,"post_contact",h);
     for(ptrdiff_t i=0;i<VORO_NP(s);i++){
      ptrdiff_t j=gfs_weight_base_find(map,n,PINDX(p+i));if(j<0)MPI_Abort(MPI_COMM(s),99);
      double *q=contact_snapshot+4*j;
      double dy=p[i].y-q[1];if(period[1]>0)dy-=period[1]*nearbyint(dy/period[1]);
      /* Suppress the harmless fmod rounding of periodic y. Actual radial
       * contact changes x or velocity as well, except at the vertical pole. */
      if(p[i].x!=q[0]||p[i].vx!=q[2]||p[i].vy!=q[3]||fabs(dy)>8*DBL_EPSILON*period[1])
       fprintf(stderr,"[MATERIAL_CONTACT_MOVE] t=%.17g attempt=%d stage=%d dt=%.17g id=%lld oldx=%.17g oldy=%.17g x=%.17g y=%.17g dx=%.17g dy=%.17g dvx=%.17g dvy=%.17g\n",
        gfs_rk4_context_time,attempt,stage,h,(long long)PINDX(p+i),q[0],q[1],(double)p[i].x,(double)p[i].y,
        p[i].x-q[0],dy,p[i].vx-q[2],p[i].vy-q[3]);
     }
    }
    for(ptrdiff_t i=0;i<VORO_NP(s);i++)after+=.5L*p[i].mass*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
    wall_repair[stage-1]=(double)(after-before);
   }else post(s);
   if(VORO_NP(s)!=n){reject=1;break;}
  }
  if(!reject){
   if(!gfs_material_evaluation_set(gfs_rk4_context_time,h,3))reject=1;
  }
  if(!reject){
#ifdef GFS_MATERIAL_VOLUME_PRIMARY
   reject=!gfs_material_ssprk2_geometry(s,xmin,ymin,xmax,ymax,gamma,padding,neighbors,cells,linked,
       volume_mode?2:0,h);
#else
   updateDenW2Pressure2DBlend(s,xmin,ymin,xmax,ymax,gamma,padding,neighbors,cells,linked,1);
   reject=gpu_geometry_trial_failed();
#endif
  }
#ifdef GFS_MATERIAL_SEARCH_EXTENT
  if(!reject&&gfs_material_search_extent(s,&xmin,&ymin,&xmax,&ymax,NULL)<0)reject=1;
#endif
  if(!reject&&material_rt){
   /* Qualification ledger: evaluate the same endpoint potential used by
    * the RHS (centroid normally; generator only in the explicitly compiled
    * point-gravity experiment). Never repair the thermal state. */
   gfs_material_constraint_step_budget(h);
   double endpoint_bound=getAccVoro2DBlend_GPU(s,xmin,ymin,xmax,ymax,VoroAccuracyOrder(s),GAS_COURANT(s),gamma,padding,neighbors,cells,linked);
   double endpoint_power;
   if(!(endpoint_bound>0)||!isfinite(endpoint_bound)||!gfs_material_potential_ledger(&potential_final,&endpoint_power))reject=1;
  }
  if(!reject&&external){
   treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
   for(ptrdiff_t i=0;i<n;i++){
    GfsMaterialStepState a={p[i].mass,{p[i].x,p[i].y,0},{p[i].vx,p[i].vy,0},p[i].weight_thermal.U,p[i].weight.W},q;double work;
    if(!gfs_material_external_kick(2,.5*h,gravity,&a,&q,&work)){reject=1;break;}
    p[i].vx=q.v[0];p[i].vy=q.v[1];external_work+=work;
   }
  }
  gfs_material_evaluation_clear();
  if(contact_audit&&!reject)gfs_material_contact_cell_diagnostic(s,contact_id,attempt,3,"endpoint",h);
  if((material_noh||material_cyl)&&VORORK4_TBPP(s)){
   /* Endpoint geometry keeps its padding alive for an optional RHS. The
    * material endpoint needs no further force; release this trial's samples. */
   my_free(VORORK4_TBPP(s));VORORK4_TBPP(s)=NULL;VORO_NPAD(s)=0;
  }
  long double E=0,E0=0,U=0,U0=0,volume_defect_sum=0;
  double minU=HUGE_VAL,minV=HUGE_VAL,volume_relative_max=0,volume_defect_max=0;
  long long volume_worst_id=-1;
  if(!reject){
   treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
   for(ptrdiff_t i=0;i<n;i++){
    ptrdiff_t j=gfs_weight_base_find(map,n,PINDX(p+i));
    if(j<0||p[i].mass!=base[j].mass||!gfs_weight_particle_valid(p+i,GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_THERMAL64)){reject=1;break;}
    E+=p[i].weight_thermal.U+.5L*p[i].mass*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
    E0+=original[i].weight_thermal.U+.5L*original[i].mass*(original[i].vx*original[i].vx+original[i].vy*original[i].vy);
    U+=p[i].weight_thermal.U;U0+=base[i].weight_thermal.U;
    minU=fmin(minU,p[i].weight_thermal.U);minV=fmin(minV,p[i].weight_thermal.V);
    if(volume_audit){double defect,relative;
     if(!gfs_material_native_volume_audit_cell((long long)PINDX(p+i),p[i].weight_thermal.V,h,&defect,&relative)){
      fprintf(stderr,"[MATERIAL_VOLUME_AUDIT] incomplete stage/ID sampling at endpoint\n");MPI_Abort(MPI_COMM(s),99);
     }
     volume_defect_sum+=defect;volume_defect_max=fmax(volume_defect_max,fabs(defect));
     if(relative>volume_relative_max){volume_relative_max=relative;volume_worst_id=(long long)PINDX(p+i);}
    }
   }
  }
  gpu_geometry_trial_end();
  if(!reject){
   if(volume_audit)fprintf(stderr,"[MATERIAL_VOLUME_STEP] t=%.17g dt=%.17g max_relative=%.17g max_absolute=%.17g signed_sum=%.17Lg worst_id=%lld\n",gfs_rk4_context_time+h,h,volume_relative_max,volume_defect_max,volume_defect_sum,volume_worst_id);
   gfs_material_native_volume_audit_end();
   if(prior_cfl_end_time==gfs_rk4_context_time){older_cfl_slope=prior_cfl_slope;older_cfl_dt=prior_cfl_dt;}
   else{older_cfl_slope=0;older_cfl_dt=0;}
   prior_cfl_dt=h;prior_cfl_bound=stage_bound[1];
   prior_cfl_end_time=gfs_rk4_context_time+h;
   prior_cfl_slope=(stage_bound[1]-stage_bound[0])/h;
   if(predict_cfl)fprintf(stderr,"[MATERIAL_CFL_PREDICT] start_bound=%.17g end_bound=%.17g prior_slope=%.17g next_slope=%.17g\n",stage_bound[0],stage_bound[1],predicted_slope,prior_cfl_slope);
   if(material_rt)fprintf(stderr,"[MATERIAL_RT_ENERGY] t=%.17g potential0=%.17g potential1=%.17g integrated_potential_rate=%.17g total_energy_residual=%.17Lg\n",gfs_rk4_context_time+h,potential_initial,potential_final,.5*h*(potential_rates[0]+potential_rates[1]),E-E0+(potential_final-potential_initial));
   if(external)fprintf(stderr,"[MATERIAL_EXTERNAL_WORK] work=%.17g residual=%.17Lg\n",external_work,E-E0-external_work);
   if(material_noh||material_cyl)fprintf(stderr,"[MATERIAL_BOUNDARY_WORK] t=%.17g work=%.17g residual=%.17Lg\n",gfs_rk4_context_time+h,.5*h*(boundary_power[0]+boundary_power[1]),E-E0-.5*h*(boundary_power[0]+boundary_power[1]));
   if(material_cyl)gfs_material_cylinder_accept(.5*h*(boundary_power[0]+boundary_power[1]),
     .5*h*(wall_power[0]+wall_power[1]),.5*wall_repair[0]+wall_repair[1]);
   fprintf(stderr,"[MATERIAL_SSP_BUDGET] t=%.17g heat=%.17g pressure_work=%.17g kinetic_work=%.17g mixing=%.17g time_increment=%.17g dU=%.17Lg closure=%.17Lg\n",gfs_rk4_context_time+h,.5*h*(ledger[0][0]+ledger[1][0]),.5*h*(ledger[0][1]+ledger[1][1]),.5*h*(ledger[0][2]+ledger[1][2]),mixing,thermal_time_increment,U-U0,U-U0-.5*h*(ledger[0][0]+ledger[1][0]+ledger[0][1]+ledger[1][1])-thermal_time_increment);
#ifdef GFS_MATERIAL_SEARCH_EXTENT
   /* Retain capacity for the next material step but do not expose it to
    * unrelated drivers between calls in the same process. */
   gfs_material_search_capacity.active=0;
#endif
   fprintf(stderr,"[MATERIAL_SSP_ACCEPT] t=%.17g dt=%.17g E=%.17Lg dE=%.17Lg Umin=%.17g Vmin=%.17g retries=%d wall_s=%.6f\n",gfs_rk4_context_time+h,h,E,E-E0,minU,minV,attempt,MPI_Wtime()-begin);gfs_material_step_support_end();free(first_rhs);free(contact_snapshot);free(first_force);free(original);free(base);free(map);return h;}
  gfs_ssp_restore(s,original,n,total,dtold,floors,nfloors,pairs);
  /* A valid second-stage CFL bound is a better first retry than halving.
   * Recompute both stages and retain every acceptance check. If it fails
   * again (or failure was not a valid CFL bound), use the existing halving.
   * Only a certified step-start evaluation may survive; rejected predictor
   * state, stage2 samples and device geometry never survive the restore. */
  /* Stage-one budget reevaluation must not consume the single directed
   * stage-two CFL retry. All rollback, stage checks and fallback halving
   * remain unchanged; never iterate directed CFL retries indefinitely. */
  int use_bound=(attempt==0||(retry_from_stage2&&!used_stage_cfl_retry))&&retry_bound>0&&retry_bound<h;
  double next_ceiling=use_bound?nextafter(retry_bound,0.0):.5*(h>0?h:ceiling);
  if(use_bound&&retry_from_stage2)used_stage_cfl_retry=1;
  /* A rejected duration must advance the retry to a different representable
   * endpoint. Otherwise an asymptotic CFL predictor can spend all attempts
   * reducing h by sub-time-ULP amounts. Use the existing conservative halving;
   * do not accept a mismatched stage budget or relax a physical constraint. */
  if(next_ceiling>0&&next_ceiling<ceiling&&
     gfs_rk4_context_time+next_ceiling==gfs_rk4_context_time+ceiling){
   next_ceiling*=.5;
   fprintf(stderr,"[MATERIAL_RETRY_TIME_PROGRESS] old=%.17g next=%.17g\n",ceiling,next_ceiling);
  }
  ceiling=next_ceiling;
  fprintf(stderr,"[MATERIAL_SSP_RETRY] ceiling=%.17g attempt=%d trial_dt=%.17g stage1_bound=%.17g stage2_bound=%.17g retry_bound=%.17g from_stage2=%d directed_bound=%d\n",ceiling,attempt,h,stage_bound[0],stage_bound[1],retry_bound,retry_from_stage2,use_bound);
 }
 gfs_rk4_failure_snapshot(s,original,n,"material_ssp_stop",0,dtold);gfs_material_native_volume_audit_end();gfs_material_step_support_end();free(first_rhs);free(contact_snapshot);free(first_force);free(original);free(base);free(map);MPI_Abort(MPI_COMM(s),99);return 0;
}

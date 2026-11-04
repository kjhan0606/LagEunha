/* Coupled prototype using original native GPU geometry and pressure,
 * shared OpenMP implicit viscosity, and the original central-force KDK
 * subcycle or free periodic drift. Explicit one-rank W0 boundaries are
 * supported. Multi-rank evolution is not qualified by the periodic tests.
 */
#include "../gfs_variational_core.h"
#include "gfs_unified_planar_audit.inc"
#include "gfs_unified_flow_audit.inc"
typedef struct {unsigned long long id;double V,U;} GfsUnifiedDriftState;
static int gfs_unified_saved_compare(const void *a,const void *b){
    unsigned long long x=((const GfsUnifiedDriftState*)a)->id,y=((const GfsUnifiedDriftState*)b)->id;
    return x<y?-1:x>y;
}
static double gfs_unified_kdk(
        SimParameters *simpar,
        void (*paddingAllTreeParticles)(SimParameters *,postype),
        Voro2D_point *(*find2DNeighborBP)(SimParameters *,int,int,int *),
        treevorork4particletype *(*find2DCellBP)(SimParameters *,int,int,int *),
        void mkLinkedList2D(SimParameters *,postype,postype,postype,postype,postype,
                           void (*)(SimParameters *,postype)),
        void (*postStage)(SimParameters *))
{
    const char *implicit=getenv("GFS_UNIFIED_IMPLICIT"),*ic=getenv("EUNHA_IC");
    const char *raw=getenv("GFS_GPU_RAW_FACES"),*finalv=getenv("GFS_GPU_FINAL_VOLUME");
    const char *wall=getenv("GFS_KEPLER_KDK_FIXEDPOLY_WALL_EVENTS");
    int periodic=gfs_unified_periodic_mode();
    int boundary_valid=periodic
        ? postStage==periodic_postStage_blend&&ic&&strcmp(ic,"kepler")
        : ic&&!strcmp(ic,"kepler")&&wall&&!strcmp(wall,"1")&&postStage==fixed_polygon_postStage_blend;
    if(NID(simpar)!=1||!implicit||strcmp(implicit,"1")||!boundary_valid||
       !raw||strcmp(raw,"1")||!finalv||strcmp(finalv,"1")||
       !gfs_native_u64_w0_kdk_mode()||GAS_EVOLMETHOD(simpar)!=3||
       !GAS_GPU_ENABLED(simpar)||GAS_Kappa(simpar)!=0||gfs_weight_native_active()||
       GAS_ENTROPY_MODE(simpar)!=0||GAS_FCENTROID(simpar)!=0||GAS_XSPHEPS(simpar)!=0||
       GAS_VISCOSITY(simpar)!=0||GAS_AlphaVis(simpar)!=0||GAS_BetaVis(simpar)!=0||
       GAS_ACCX(simpar)!=0||GAS_ACCY(simpar)!=0){
        fprintf(stderr,"[UNIFIED_KDK] requires one-rank W0, implicit common viscosity, raw final geometry, explicit supported boundary and no extra forces/stress/steering\n");
        MPI_Abort(MPI_COMM(simpar),99);
    }
    double begin=MPI_Wtime(),geometry_seconds=0,pressure_seconds=0,viscosity_seconds=0,drift_seconds=0;
    double gamma=GAS_GAMMA(simpar),courant=GAS_COURANT(simpar),OoA=VoroAccuracyOrder(simpar);
    double cell=BASICCELL_CELLWIDTH(simpar),xmin=Xmin_HydroExam(simpar)-cell,ymin=Ymin_HydroExam(simpar)-cell;
    double xmax=Xmax_HydroExam(simpar)+cell,ymax=Ymax_HydroExam(simpar)+cell;
    double eta=.01,cx=3,cy=3,eps=.01175;const char *env;
    if((env=getenv("GFS_KDK_KEPLER_SUBCYCLE_ETA")))eta=atof(env);
    if((env=getenv("EUNHA_KEPLER_CX")))cx=atof(env);
    if((env=getenv("EUNHA_KEPLER_CY")))cy=atof(env);
    if((env=getenv("EUNHA_KEPLER_EPS")))eps=atof(env);
    if(!(eta>0)||!isfinite(eta))MPI_Abort(MPI_COMM(simpar),99);
    postype warm_ax=0,warm_ay=0;kepler_accel(0,0,&warm_ax,&warm_ay);if(!periodic)gfs_kdk_fixedpoly_wall_init();
    static int cache_ready,last_step;
    static double last_end,geometry_key[6];
    static SimParameters *cached_simpar;
    static void (*cached_boundary)(SimParameters *);
    double key[6]={xmin,ymin,xmax,ymax,cell,OoA};
    const char *reuse_env=getenv("GFS_UNIFIED_REUSE_GEOMETRY");
    int reuse=(!reuse_env||strcmp(reuse_env,"0"))&&cache_ready&&cached_simpar==simpar&&
        cached_boundary==postStage&&last_step+1==gfs_rk4_context_step&&last_end==gfs_rk4_context_time&&
        !memcmp(key,geometry_key,sizeof(key))&&gpu_unified_geometry_reusable(simpar);
    cache_ready=0;
    double mark=MPI_Wtime(),dt;
    if(reuse){
        /* The last kick changed v. Recompute CFL, fit, and pressure from
         * current states while retaining only the unchanged geometry. */
        dt=gpu_unified_pressure_refresh(simpar);
    }else{
        updateDenW2Pressure2DBlend(simpar,xmin,ymin,xmax,ymax,gamma,
            paddingAllTreeParticles,find2DNeighborBP,find2DCellBP,mkLinkedList2D,1);
        geometry_seconds+=MPI_Wtime()-mark;mark=MPI_Wtime();
        dt=getAccVoro2DBlend_GPU(simpar,xmin,ymin,xmax,ymax,OoA,courant,gamma,
            paddingAllTreeParticles,find2DNeighborBP,find2DCellBP,mkLinkedList2D);
    }
    pressure_seconds+=MPI_Wtime()-mark;
    const char *factor_env=getenv("GFS_UNIFIED_DT_FACTOR");
    double dt_factor=factor_env?atof(factor_env):1;
    if(!(dt_factor>0)||dt_factor>1||!isfinite(dt_factor))MPI_Abort(MPI_COMM(simpar),99);
    dt*=dt_factor;
    if(gfs_rk4_step_ceiling>0)dt=fmin(dt,gfs_rk4_step_ceiling);
    if(!(dt>0)||!isfinite(dt)||(!periodic&&VORO_NPAD(simpar)!=0))MPI_Abort(MPI_COMM(simpar),99);
    int n=VORO_NP(simpar),bad=0;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    gfs_unified_planar_audit(simpar,"initial_pressure",gfs_rk4_context_time,gamma);
    gfs_unified_flow_audit(simpar,gfs_rk4_context_time,gamma);
    static int initial_audit_done=0;
    const char *initial_audit=getenv("GFS_UNIFIED_INITIAL_AUDIT");
    if(initial_audit&&!strcmp(initial_audit,"1")&&!initial_audit_done){
        if(gpu_unified_viscosity_rate_audit(simpar,gamma))MPI_Abort(MPI_COMM(simpar),99);
        initial_audit_done=1;
        fprintf(stderr,"[UNIFIED_INITIAL_AUDIT] t=%.17g accepted_steps_this_process=0 pressure_only_acceleration=1\n",gfs_rk4_context_time);
        const char *only=getenv("GFS_UNIFIED_INITIAL_ONLY");
        if(only&&!strcmp(only,"1")){MPI_Finalize();exit(0);}
    }
    static double E_reference=0;static int reference_ready=0;
    if(!reference_ready){
        long double E=0;for(int i=0;i<n;i++)E+=p[i].weight_thermal.U+.5L*p[i].mass*(p[i].vx*p[i].vx+p[i].vy*p[i].vy)+p[i].mass*kepler_phi(p[i].x,p[i].y);
        E_reference=(double)E;reference_ready=1;
    }
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){p[i].vx+=.5*dt*p[i].ax;p[i].vy+=.5*dt*p[i].ay;}
    gfs_unified_planar_audit(simpar,"first_pressure_kick",gfs_rk4_context_time,gamma);
    GfsViscosityStats D0,D1;mark=MPI_Wtime();
    int rc=gpu_unified_viscosity_half_step(simpar,gamma,.5*dt,&D0);
    if(rc){fprintf(stderr,"[UNIFIED_KDK] first viscosity failed rc=%d\n",rc);MPI_Abort(MPI_COMM(simpar),99);}
    viscosity_seconds+=MPI_Wtime()-mark;
    gfs_unified_planar_audit(simpar,"first_viscosity",gfs_rk4_context_time,gamma);
    static GfsUnifiedDriftState *saved;static int capacity=0;
    if(n>capacity){void *q=realloc(saved,(size_t)n*sizeof(*saved));if(!q)MPI_Abort(MPI_COMM(simpar),99);saved=q;capacity=n;}
    for(int i=0;i<n;i++)saved[i]=(GfsUnifiedDriftState){(unsigned long long)PINDX(p+i),p[i].weight_thermal.V,p[i].weight_thermal.U};
    qsort(saved,n,sizeof(*saved),gfs_unified_saved_compare);
    for(int i=1;i<n;i++)if(saved[i].id==saved[i-1].id)MPI_Abort(MPI_COMM(simpar),99);
    mark=MPI_Wtime();double wall_Jx=0,wall_Jy=0,wall_Jz=0,wall_roundoff=0;unsigned long long impacts=0;
#pragma omp parallel for schedule(static) reduction(max:bad) reduction(+:wall_Jx,wall_Jy,wall_Jz,wall_roundoff,impacts)
    for(int i=0;i<n;i++){
        if(periodic){p[i].x+=dt*p[i].vx;p[i].y+=dt*p[i].vy;continue;}
        GfsKdkFixedPolyWallStats w;double x,y,vx,vy;
        int sub=gfs_kdk_lts_kepler_drift(p[i].x,p[i].y,p[i].vx,p[i].vy,dt,p[i].mass,1,eta,cx,cy,eps,0,0,&x,&y,&vx,&vy,&w);
        if(sub<0){bad=1;continue;}p[i].x=x;p[i].y=y;p[i].vx=vx;p[i].vy=vy;
        impacts+=w.impacts;wall_Jx+=w.wall_jx;wall_Jy+=w.wall_jy;wall_Jz+=w.wall_torque;wall_roundoff+=w.kinetic_roundoff;
    }
    if(bad){fprintf(stderr,"[UNIFIED_KDK] drift failure\n");MPI_Abort(MPI_COMM(simpar),99);}
    postStage(simpar);drift_seconds+=MPI_Wtime()-mark;
    if(VORO_NP(simpar)!=n)MPI_Abort(MPI_COMM(simpar),99);
    mark=MPI_Wtime();
    updateDenW2Pressure2DBlend(simpar,xmin,ymin,xmax,ymax,gamma,
        paddingAllTreeParticles,find2DNeighborBP,find2DCellBP,mkLinkedList2D,dt);
    geometry_seconds+=MPI_Wtime()-mark;p=(void*)VORORK4_TBP(simpar);
    double max_adiabatic_identity_error=0;
#pragma omp parallel for schedule(static) reduction(max:bad,max_adiabatic_identity_error)
    for(int i=0;i<n;i++){
        GfsUnifiedDriftState key={(unsigned long long)PINDX(p+i),0,0};
        const GfsUnifiedDriftState *old=bsearch(&key,saved,n,sizeof(*saved),gfs_unified_saved_compare);double U1;
        if(!old||!gfs_variational_adiabatic(old->U,old->V,p[i].weight_thermal.V,gamma,&U1)){bad=1;continue;}
        double identity=U1/old->U*pow(p[i].weight_thermal.V/old->V,gamma-1)-1;
        max_adiabatic_identity_error=fmax(max_adiabatic_identity_error,fabs(identity));
        p[i].weight_thermal.U=U1;p[i].ie=(float)U1;
        double V=p[i].weight_thermal.V,rho=p[i].mass/V,P=(gamma-1)*U1/V;
        p[i].volume=(float)V;p[i].den=(float)rho;p[i].pressure=(float)P;p[i].csound=(float)sqrt(gamma*P/rho);
    }
    if(bad)MPI_Abort(MPI_COMM(simpar),99);
    gfs_unified_planar_audit(simpar,"drift_adiabat",gfs_rk4_context_time+dt,gamma);
    /* Consume the newly built GPU geometry through the normal native wrapper.
     * Pressure is refreshed again after dissipation changes U. */
    mark=MPI_Wtime();
    (void)getAccVoro2DBlend_GPU(simpar,xmin,ymin,xmax,ymax,OoA,courant,gamma,
        paddingAllTreeParticles,find2DNeighborBP,find2DCellBP,mkLinkedList2D);
    pressure_seconds+=MPI_Wtime()-mark;mark=MPI_Wtime();
    gfs_unified_planar_audit(simpar,"final_pressure_before_heat",gfs_rk4_context_time+dt,gamma);
    rc=gpu_unified_viscosity_half_step(simpar,gamma,.5*dt,&D1);
    if(rc){fprintf(stderr,"[UNIFIED_KDK] second viscosity failed rc=%d\n",rc);MPI_Abort(MPI_COMM(simpar),99);}
    viscosity_seconds+=MPI_Wtime()-mark;mark=MPI_Wtime();
    gfs_unified_planar_audit(simpar,"second_viscosity",gfs_rk4_context_time+dt,gamma);
    (void)gpu_unified_pressure_refresh(simpar);pressure_seconds+=MPI_Wtime()-mark;
    p=(void*)VORORK4_TBP(simpar);
    gfs_unified_planar_audit(simpar,"final_pressure_after_heat",gfs_rk4_context_time+dt,gamma);
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){p[i].vx+=.5*dt*p[i].ax;p[i].vy+=.5*dt*p[i].ay;}
    gfs_unified_planar_audit(simpar,"final_pressure_kick",gfs_rk4_context_time+dt,gamma);
    long double E=0;double Umin=HUGE_VAL,Vmin=HUGE_VAL;
    for(int i=0;i<n;i++){
        double U=p[i].weight_thermal.U,V=p[i].weight_thermal.V;
        if(!(U>0)||!(V>0)||!isfinite(p[i].vx)||!isfinite(p[i].vy)||p[i].weight.W!=0||p[i].weight.pi!=0)bad=1;
        Umin=fmin(Umin,U);Vmin=fmin(Vmin,V);
        E+=U+.5L*p[i].mass*(p[i].vx*p[i].vx+p[i].vy*p[i].vy)+p[i].mass*kepler_phi(p[i].x,p[i].y);
    }
    if(bad)MPI_Abort(MPI_COMM(simpar),99);
    last_step=gfs_rk4_context_step;last_end=gfs_rk4_context_time+dt;
    gfs_unified_flow_audit(simpar,last_end,gamma);
    memcpy(geometry_key,key,sizeof(key));cached_simpar=simpar;cached_boundary=postStage;cache_ready=1;
    fprintf(stderr,"[UNIFIED_GEOMETRY_REUSE] reused=%d new_geometry_builds=%d\n",reuse,reuse?1:2);
    fprintf(stderr,"[UNIFIED_KDK] step=%d t=%.17g dt=%.17g N=%d Etot=%.17Lg dErel=%.9e Umin=%.9e Vmin=%.9e heat=%.9e viscosity_energy_defect=%.9e CG=%d/%d active=%d/%d wall_hits=%llu wall_J=(%.9e,%.9e,%.9e) wall_roundoff=%.9e geometry_s=%.6f pressure_s=%.6f viscosity_s=%.6f drift_s=%.6f total_s=%.6f\n",
        gfs_rk4_context_step,gfs_rk4_context_time+dt,dt,n,E,((double)E-E_reference)/fabs(E_reference),Umin,Vmin,D0.heat+D1.heat,D0.energy_defect+D1.energy_defect,D0.iterations,D1.iterations,D0.active_iterations,D1.active_iterations,impacts,wall_Jx,wall_Jy,wall_Jz,wall_roundoff,geometry_seconds,pressure_seconds,viscosity_seconds,drift_seconds,MPI_Wtime()-begin);
    fprintf(stderr,"[UNIFIED_THERMAL_AUDIT] gamma=%.17g max_adiabat_identity_error=%.17g viscosity_min_heat=(%.17g,%.17g) viscosity_equation_residual=(%.17g,%.17g)\n",gamma,max_adiabatic_identity_error,D0.min_heat,D1.min_heat,D0.relative_equation_residual,D1.relative_equation_residual);
    return dt;
}

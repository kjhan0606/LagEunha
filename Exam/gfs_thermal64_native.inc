/* Included after all legacy helpers and before the native RK4 driver. */
static void gfs_t64_dispatch_guard(SimParameters *simpar){
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    int enabled=gfs_thermal64_requested(),tagged=0;
    for(ptrdiff_t k=0;k<VORO_NP(simpar);k++)if(GFS_T64_TAGGED(p[k]))tagged++;
    if(!enabled){
        if(tagged){fprintf(stderr,"[THERMAL64] tagged state cannot enter legacy RK4\n");MPI_Abort(MPI_COMM(simpar),99);}
        return;
    }
    /* Moment initialization rewrites the legacy thermal mirror and exits
     * before thermal64 initialization. Until a joint energy/state contract
     * exists, allowing both flags can serialize a stale authoritative U64. */
    const char *moment_init=getenv("GFS_MOMENT_INIT_ONLY");
    if(moment_init && atoi(moment_init)>0){
        fprintf(stderr,"[THERMAL64] moment initialization cannot share the thermal64 payload contract\n");
        MPI_Abort(MPI_COMM(simpar),99);
    }
    gfs_es_config();
#ifdef USE_CUDA
    const char *direct=getenv("GFS_RK4_DIRECT_U"),*ssp=getenv("GFS_SSPRK3");
    const char *local=getenv("GFS_GPU_LOCAL_PERIODIC_SYNC");
    const char *tx=getenv("GFS_RK4_TRANSACTION"),*txdiag=getenv("GFS_RK4_TX_ENTROPY_DIAGNOSTIC");
    const char *rot=getenv("GFS_LAGUERRE_ROTATION"),*charge=getenv("GFS_FACE_CHARGE_LIMIT");
    const char *geomface=getenv("GFS_GEOM_FACE_VEL"),*initonly=getenv("GFS_THERMAL64_INIT_ONLY");
    int init_gate=initonly && atoi(initonly)>0;
    if(NID(simpar)==1 && GAS_GPU_ENABLED(simpar) && gpu_varpress_w0_enabled() &&
       (init_gate || (tx && atoi(tx)>0 && (!txdiag || atoi(txdiag)==0))) &&
       (!rot || atoi(rot)==0) && (!charge || atoi(charge)==0) && (!geomface || atoi(geomface)==0) &&
       GAS_AVMODE(simpar)==5 && GAS_Kappa(simpar)==0 && GAS_ENTROPY_MODE(simpar)==0 &&
       GAS_XSPHEPS(simpar)==0 && GAS_FCENTROID(simpar)==0 && GAS_VISCOSITY(simpar)==0 &&
       GAS_AlphaVis(simpar)==0 && GAS_BetaVis(simpar)==0 && GAS_CDAMAX(simpar)==0 &&
       GAS_HYPERVALPHA(simpar)==0 && sedov_phase1_on() && !gfs_dual_energy_eta() &&
       !gfs_es_on && !gfs_half_on && !gfs_dynw2_on() && !gfs_w2lin_on() && !sedov_lagvol_on() && !gfs_vaud_on() &&
       (!direct || atoi(direct)==0) && (!ssp || atoi(ssp)==0) && (!local || atoi(local)==0) &&
       getenv("EUNHA_IC") && !strcmp(getenv("EUNHA_IC"),"kepler") &&
       getenv("GFS_GPU_RAW_FACES") && atoi(getenv("GFS_GPU_RAW_FACES"))==1 &&
       getenv("GFS_GPU_FINAL_VOLUME") && atoi(getenv("GFS_GPU_FINAL_VOLUME"))==1 &&
       TVORORK4_DDINFO(simpar)[0].n_size==sizeof(*p) && (tagged==0 || tagged==VORO_NP(simpar)))return;
#endif
    fprintf(stderr,"[THERMAL64] unsupported or mixed payload contract; smooth periodic one-rank native W0 only\n");
    MPI_Abort(MPI_COMM(simpar),99);
}
static void gfs_t64_eos(SimParameters *simpar){
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
        double U=GFS_T64_U(p[k]),V=GFS_T64_V(p[k]);
        double P=(GAS_GAMMA(simpar)-1.)*U/V,rho=p[k].mass/V;
        if(!GFS_T64_TAGGED(p[k]) || !(U>0) || !(V>0) || !(P>0) || !isfinite(P) || !(rho>0) || !isfinite(rho)){
            fprintf(stderr,"[THERMAL64] invalid actual thermal state before EOS\n");MPI_Abort(MPI_COMM(simpar),99);
        }
        p[k].ie=(float)U;p[k].volume=(float)V;p[k].den=(float)rho;
        p[k].pressure=(float)P;p[k].csound=(float)sqrt(GAS_GAMMA(simpar)*P/rho);
    }
}
static int gfs_t64_initialize(SimParameters *simpar,
    void (*padding)(SimParameters*,postype),
    Voro2D_point *(*neighbors)(SimParameters*,int,int,int*),
    treevorork4particletype *(*cell)(SimParameters*,int,int,int*),
    void (*linked)(SimParameters*,postype,postype,postype,postype,postype,void(*)(SimParameters*,postype))){
    if(!gfs_thermal64_requested())return 0;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    if(VORO_NP(simpar)>0 && !GFS_T64_TAGGED(p[0])){
        for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
            GFS_T64_U(p[k])=(double)p[k].ie;GFS_T64_V(p[k])=(double)p[k].volume;
            GFS_T64_RHS(p[k])=0;p[k].stress.E_inv_yy=GFS_T64_TAG;
        }
        /* Exactly promote accepted heat, then obtain actual double geometry
         * at those same positions. No entropy-based heat replacement. */
        exam2dUpdateVol(simpar,padding,neighbors,cell,linked);
        VORORK4_TBPP(simpar)=NULL;VORO_BASICCELL(simpar)=NULL;
        gfs_t64_eos(simpar);
        fprintf(stderr,"[THERMAL64_INIT] accepted U promoted exactly; actual double V, unchanged particle stride=%zu\n",sizeof(*p));
    }else gfs_t64_eos(simpar);
    p=(void*)VORORK4_TBP(simpar);
    for(ptrdiff_t k=0;k<VORO_NP(simpar);k++)GFS_T64_ANCHOR(p[k])=GFS_T64_U(p[k]);
    static int init_snapshot_done=0;
    const char *only=getenv("GFS_THERMAL64_INIT_ONLY");
    if(only && atoi(only)>0 && !init_snapshot_done){
        init_snapshot_done=1;
        const char *step_s=getenv("KEPLER_RESTART_STEP");
        const char *time_s=getenv("GFS_THERMAL64_INIT_TIME");
        char *step_end=NULL,*time_end=NULL;
        long init_step=step_s?strtol(step_s,&step_end,10):-1;
        double init_time=time_s?strtod(time_s,&time_end):NAN;
        if(!step_s || !time_s || step_end==step_s || *step_end || init_step<0 || init_step>2147483647L ||
           time_end==time_s || *time_end || !isfinite(init_time)){
            fprintf(stderr,"[THERMAL64_INIT_ONLY] exact restart step and time are required\n");
            MPI_Abort(MPI_COMM(simpar),99);
        }
        gfs_rk4_set_step_context((int)init_step,init_time,0.0);
        gfs_rk4_failure_snapshot(simpar,p,(int)VORO_NP(simpar),"thermal64_initial",GAS_dtold(simpar),GAS_dtold(simpar));
        fprintf(stderr,"[THERMAL64_INIT_ONLY] zero physical time advanced; step=%ld time=%.17g\n",init_step,init_time);
        MPI_Finalize();exit(0);
    }
    return 1;
}

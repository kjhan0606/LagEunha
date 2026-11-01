/* Experimental native SSPRK(3,3), fixed masses and W=0. Included after the
 * shared geometry, gravity and diagnostic helpers in exam.c. Defaults OFF.
 * Rejected attempts never reach physical cumulative ledgers or output clocks.
 */
#ifdef USE_CUDA
static long gfs_ssp_rhs_count=0,gfs_ssp_accepted_count=0;

static long gfs_ssp_limit(const char *name,long fallback){
    const char *s=getenv(name); if(!s || !s[0]) return fallback;
    char *end=NULL; long v=strtol(s,&end,10);
    if(*end || v<1 || v>1000000000L){ fprintf(stderr,"[SSPRK3] invalid %s\n",name); MPI_Abort(MPI_COMM_WORLD,99); }
    return v;
}

static int gfs_ssp_valid(SimParameters *simpar,int rhs){
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    int bad=0,global=0;
    for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
        if(!(p[k].ie>0) || !(p[k].mass>0) || !(p[k].volume>0) ||
           !isfinite(p[k].ie) || !isfinite(p[k].mass) || !isfinite(p[k].volume) ||
           !isfinite(p[k].x) || !isfinite(p[k].y) || !isfinite(p[k].vx) || !isfinite(p[k].vy) ||
           p[k].w2!=0 || p[k].mass!=p[k].rk4.k1K) bad=1;
        if(rhs && (!isfinite(p[k].ax) || !isfinite(p[k].ay) || !isfinite(p[k].die))) bad=1;
    }
    MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,MPI_COMM(simpar));
    return !global;
}

static int gfs_ssp_eos(SimParameters *simpar){
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    int bad=0,global=0;
    for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
        double pressure=(GAS_GAMMA(simpar)-1.)*(double)p[k].ie/p[k].volume;
        p[k].den=p[k].mass/p[k].volume;p[k].pressure=pressure;
        p[k].csound=sqrt(GAS_GAMMA(simpar)*pressure/p[k].den);
        if(!(p[k].den>0) || !(p[k].pressure>0) || !(p[k].csound>0) ||
           !isfinite(p[k].den) || !isfinite(p[k].pressure) || !isfinite(p[k].csound)) bad=1;
    }
    MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,MPI_COMM(simpar));
    return !global;
}

static void gfs_ssp_restore(SimParameters *simpar,const void *backup,ptrdiff_t np,
    ptrdiff_t total,float dtold,double floors,int nfloors,long long pairs){
    /* pmigrate reads decomposition pivots; it changes particle ownership and
     * allocations, not the pivots. Each rank restores its own original full
     * payload/count. No allgather of global particle state is needed. */
    if(VORORK4_TBPP(simpar)) my_free(VORORK4_TBPP(simpar));
    VORORK4_TBPP(simpar)=NULL; VORO_NPAD(simpar)=0;
    my_free(VORORK4_TBP(simpar));
    VORORK4_TBP(simpar)=my_malloc((size_t)(np?np:1)*sizeof(treevorostressrk4particletype));
    memcpy(VORORK4_TBP(simpar),backup,(size_t)np*sizeof(treevorostressrk4particletype));
    VORO_NP(simpar)=np; VORO_TNP(simpar)=total; NPSUM(simpar)=total;
    GAS_dtold(simpar)=dtold;
    gfs_stage_floor_l=floors; gfs_stage_floor_n_l=nfloors; gfs_pair_face_ends=pairs;
    /* All reusable device geometry, gradients, density seeds and CSR indices
     * refer to the discarded attempt. Rebuild them from restored host state. */
    gpu_free(gpu_get_context());
    VORO_BASICCELL(simpar)=NULL;
}

static void gfs_ssp_ledger(SimParameters *simpar,double z[7]){
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    double local[7]={0,0,0,0,0,0,0};
    postype cx=SIMBOX(simpar).x.max*.5,cy=SIMBOX(simpar).y.max*.5;
    const char *ic=getenv("EUNHA_IC"); if(ic && !strcmp(ic,"kepler")) kepler_center(&cx,&cy);
    for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
        double m=p[k].mass,vx=p[k].vx,vy=p[k].vy;
        local[0]+=p[k].ie+.5*m*(vx*vx+vy*vy);
        local[1]+=m*kepler_phi(p[k].x,p[k].y); local[2]+=m;
        local[3]+=m*vx;local[4]+=m*vy;
        local[5]+=m*((p[k].x-cx)*vy-(p[k].y-cy)*vx);local[6]+=p[k].ie;
    }
    MPI_Allreduce(local,z,7,MPI_DOUBLE,MPI_SUM,MPI_COMM(simpar));
}

static double gfs_ssprk3_step(SimParameters *simpar,
    void (*padding)(SimParameters*,postype),
    double (*measureW2)(SimParameters*,postype,postype,postype),
    Voro2D_point *(*neighbors)(SimParameters*,int,int,int*),
    treevorork4particletype *(*cell)(SimParameters*,int,int,int*),
    int (*target)(treevorork4particletype*,postype,postype),
    void (*linked)(SimParameters*,postype,postype,postype,postype,postype,void(*)(SimParameters*,postype)),
    void (*post)(SimParameters*)){
    (void)measureW2;
    gfs_es_config();
    const char *raw=getenv("GFS_GPU_RAW_FACES"),*fv=getenv("GFS_GPU_FINAL_VOLUME");
    const char *ic=getenv("EUNHA_IC");
    if(!GAS_GPU_ENABLED(simpar) || GAS_AVMODE(simpar)!=5 || !GAS_USEMUSCL(simpar) ||
       !sedov_phase1_on() || GAS_Kappa(simpar)!=0 || GAS_ENTROPY_MODE(simpar)!=0 ||
       GAS_XSPHEPS(simpar)!=0 || GAS_FCENTROID(simpar)!=0 || GAS_VISCOSITY(simpar)!=0 ||
       GAS_AlphaVis(simpar)!=0 || GAS_BetaVis(simpar)!=0 || GAS_CDAMAX(simpar)!=0 ||
       GAS_ACCX(simpar)!=0 || GAS_ACCY(simpar)!=0 ||
       GAS_HYPERVALPHA(simpar)!=0 || gfs_dual_energy_eta()>0 || gfs_es_on || gfs_half_on ||
       gfs_dynw2_on() || gfs_w2lin_on() || sedov_lagvol_on() || gfs_vaud_on() ||
       !raw || atoi(raw)<=0 || !fv || atoi(fv)<=0 || post!=periodic_postStage_blend ||
       !ic || (strcmp(ic,"kepler") && strcmp(ic,"shock_probe") && strcmp(ic,"affine_probe")) ||
       TVORORK4_DDINFO(simpar)[0].n_size!=sizeof(treevorostressrk4particletype)){
        fprintf(stderr,"[SSPRK3] unsupported contract: requires periodic native GPU phase1 av5/MUSCL, fixed masses/W0, raw complete geometry, no steering/viscosity/entropy repair/dynamic weights\n");
        MPI_Abort(MPI_COMM(simpar),99);
    }
    const long max_rhs=gfs_ssp_limit("GFS_SSP_MAX_RHS",256),max_steps=gfs_ssp_limit("GFS_SSP_MAX_STEPS",64);
    const ptrdiff_t np0=VORO_NP(simpar),ntotal=VORO_TNP(simpar);
    size_t bytes=(size_t)np0*sizeof(treevorostressrk4particletype);
    void *backup=malloc(bytes?bytes:1);
    if(!backup) MPI_Abort(MPI_COMM(simpar),99);
    memcpy(backup,VORORK4_TBP(simpar),bytes);
    const float dtold=GAS_dtold(simpar);
    const double floor0=gfs_stage_floor_l;
    const int nfloor0=gfs_stage_floor_n_l;
    const long long pair0=gfs_pair_face_ends;
    const double t0=gfs_rk4_context_time;
    double e0[7],e1[7];gfs_ssp_ledger(simpar,e0);
    if(gfs_ssp_accepted_count>=max_steps || gfs_ssp_rhs_count>=max_rhs) goto budget;
    double ceiling=gfs_rk4_step_ceiling>0?gfs_rk4_step_ceiling:1.e100;
    static int forced_done=0;
    for(int attempt=0;attempt<=8;attempt++){
        double h=ceiling,retry_bound=1.e100,hyd=0,wg=0,hab=0;
        int reject=0;
        treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
        for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
            /* Base state and unwrapped displacement migrate with each cell. */
            p[k].rk4.k1x=p[k].x;p[k].rk4.k1y=p[k].y;
            p[k].rk4.k1vx=p[k].vx;p[k].rk4.k1vy=p[k].vy;
            p[k].rk4.k1ie=p[k].ie;p[k].rk4.k1K=p[k].mass;
            p[k].rk4.k2x=0;p[k].rk4.k2y=0;
        }
        gpu_geometry_trial_begin();
        for(int stage=0;stage<3;stage++){
            if(gfs_ssp_rhs_count>=max_rhs){
                gpu_geometry_trial_end();
                gfs_ssp_restore(simpar,backup,np0,ntotal,dtold,floor0,nfloor0,pair0);
                goto budget;
            }
            if(!gfs_ssp_valid(simpar,0)){reject=1;break;}
            /* Prepare an actual stage state: complete geometry, then EOS on
             * all real cells, then the ordinary halo/gradient/force path.
             * This explicitly avoids using the previous stage's rho/P cache. */
            exam2dUpdateVol(simpar,padding,neighbors,cell,linked);
            VORORK4_TBPP(simpar)=NULL;VORO_BASICCELL(simpar)=NULL;
            if(gpu_geometry_trial_failed() || !gfs_ssp_valid(simpar,0)){reject=1;break;}
            if(!gfs_ssp_eos(simpar)){reject=1;break;}
            const postype width=gfs_kepler_halo_width(BASICCELL_CELLWIDTH(simpar));
            const postype xmin=Xmin_HydroExam(simpar)-width,ymin=Ymin_HydroExam(simpar)-width;
            const postype xmax=Xmax_HydroExam(simpar)+width,ymax=Ymax_HydroExam(simpar)+width;
            updateDenW2Pressure2DBlend(simpar,xmin,ymin,xmax,ymax,GAS_GAMMA(simpar),padding,neighbors,cell,linked,stage?h:1);
            VORO_BASICCELL(simpar)=NULL;
            if(gpu_geometry_trial_failed()){reject=1;break;}
            double cfl=getAccVoro2DBlend_GPU(simpar,xmin,ymin,xmax,ymax,VoroAccuracyOrder(simpar),
                GAS_COURANT(simpar),GAS_GAMMA(simpar),padding,neighbors,cell,linked);
            VORORK4_TBPP(simpar)=NULL;
            gfs_ssp_rhs_count++;
            if(!gfs_ssp_valid(simpar,1)){reject=1;break;}
            const char *shadow_only=getenv("GFS_ALE_SHADOW_ONLY");
            if(stage==0 && shadow_only && atoi(shadow_only)>0){
                const char *shadow_call=getenv("GFS_ALE_SHADOW_CALL");
                if(!shadow_call || atoi(shadow_call)!=1) MPI_Abort(MPI_COMM(simpar),99);
                if(MYID(simpar)==0) fprintf(stderr,"[ALE_SHADOW_ONLY] current-stage RHS captured; no evolution/accepted step/output state\n");
                gpu_geometry_trial_end();free(backup);
                MPI_Barrier(MPI_COMM(simpar));MPI_Finalize();exit(0);
            }
            p=(void*)VORORK4_TBP(simpar);
            double hu_local=1.e100,hu=1.e100;
            for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
                double ud=p[k].die-p[k].mass*(p[k].vx*p[k].ax+p[k].vy*p[k].ay);
                if(ud<0 && -p[k].ie/ud<hu_local) hu_local=-p[k].ie/ud;
            }
            MPI_Allreduce(&hu_local,&hu,1,MPI_DOUBLE,MPI_MIN,MPI_COMM(simpar));
            double orbit=kepler_dt_limit(simpar,p,sizeof(*p),(int)VORO_NP(simpar));
            double bound=fmin(cfl,fmin(orbit,.9*hu));
            if(stage==0) h=fmin(h,bound);
            retry_bound=fmin(retry_bound,bound);
            if(!(h>0) || !isfinite(h) || !(t0+h>t0)) {reject=1;break;}
            if(h>bound*(1.+16.*2.2204460492503131e-16)){
                reject=1;
                if(MYID(simpar)==0) fprintf(stderr,"[SSP_REJECT] step=%d attempt=%d stage=%d h=%.17g bound=%.17g reason=FE_bound\n",gfs_rk4_context_step,attempt,stage+1,h,bound);
                break;
            }
            char label[80];snprintf(label,sizeof(label),"ssp_a%d_s%d_rhs",attempt,stage+1);
            gfs_rk4_context_time=t0+(stage==0?0:stage==1?h:.5*h);
            gfs_rk4_stage_diagnostic(simpar,label,h,1);gfs_rk4_context_time=t0;
            double rates_local[3]={0,0,0},rates[3];
            for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
                postype gx=0,gy=0;kepler_accel(p[k].x,p[k].y,&gx,&gy);
                rates_local[0]+=p[k].die;rates_local[1]+=p[k].mass*(p[k].vx*gx+p[k].vy*gy);
                rates_local[2]+=fabs(p[k].die);
            }
            MPI_Allreduce(rates_local,rates,3,MPI_DOUBLE,MPI_SUM,MPI_COMM(simpar));
            const double weight=stage<2?1./6.:2./3.;
            hyd+=weight*h*rates[0];wg+=weight*h*rates[1];hab+=weight*h*rates[2];
            /* Forward-Euler endpoint, before any convex combination or EOS. */
            for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
                postype gx=0,gy=0;kepler_accel(p[k].x,p[k].y,&gx,&gy);
                double ud=p[k].die-p[k].mass*(p[k].vx*p[k].ax+p[k].vy*p[k].ay);
                p[k].rk4.k2x+=h*p[k].vx;p[k].rk4.k2y+=h*p[k].vy;
                p[k].x=p[k].rk4.k1x+p[k].rk4.k2x;p[k].y=p[k].rk4.k1y+p[k].rk4.k2y;
                p[k].vx+=h*(p[k].ax+GAS_ACCX(simpar)+gx);
                p[k].vy+=h*(p[k].ay+GAS_ACCY(simpar)+gy);
                p[k].ie=(postype)p[k].ie+h*ud;
            }
            if(!gfs_ssp_valid(simpar,0)){reject=1;break;}
            post(simpar);
            /* Stage1 FE endpoint is the next RHS input and is checked there.
             * Stages2/3 have separate unmixed FE geometries to validate. */
            if(stage>0){
                exam2dUpdateVol(simpar,padding,neighbors,cell,linked);
                VORORK4_TBPP(simpar)=NULL;VORO_BASICCELL(simpar)=NULL;
                if(gpu_geometry_trial_failed() || !gfs_ssp_valid(simpar,0)){reject=1;break;}
                p=(void*)VORORK4_TBP(simpar);
                double b=stage==1?.25:2./3.,a=1.-b;
                for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
                    p[k].rk4.k2x*=b;p[k].rk4.k2y*=b;
                    p[k].x=p[k].rk4.k1x+p[k].rk4.k2x;p[k].y=p[k].rk4.k1y+p[k].rk4.k2y;
                    p[k].vx=a*p[k].rk4.k1vx+b*p[k].vx;
                    p[k].vy=a*p[k].rk4.k1vy+b*p[k].vy;
                    p[k].ie=a*p[k].rk4.k1ie+b*p[k].ie;
                }
                if(!gfs_ssp_valid(simpar,0)){reject=1;break;}
                post(simpar);
            }
            const char *forced=getenv("GFS_SSP_FORCE_REJECT_ONCE");
            if(stage==0 && !forced_done && forced && atoi(forced)>0){
                forced_done=1;reject=1;
                if(MYID(simpar)==0) fprintf(stderr,"[SSP_REJECT] step=%d attempt=%d stage=1 reason=verification_injection\n",gfs_rk4_context_step,attempt);
                break;
            }
        }
        if(!reject){
            exam2dUpdateVol(simpar,padding,neighbors,cell,linked);
            VORORK4_TBPP(simpar)=NULL;VORO_BASICCELL(simpar)=NULL;
            if(gpu_geometry_trial_failed() || !gfs_ssp_valid(simpar,0)) reject=1;
            if(!reject && !gfs_ssp_eos(simpar)) reject=1;
        }
        if(!reject){
            p=(void*)VORORK4_TBP(simpar);
            double min_local=1.e100,minU;
            for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
                if(p[k].ie<min_local) min_local=p[k].ie;
            }
            MPI_Allreduce(&min_local,&minU,1,MPI_DOUBLE,MPI_MIN,MPI_COMM(simpar));
            gfs_ssp_ledger(simpar,e1);gfs_ssp_accepted_count++;
            if(MYID(simpar)==0){
                static FILE *ledger=NULL;
                if(!ledger){ledger=fopen("ssprk3_steps.tsv","wx");if(!ledger) MPI_Abort(MPI_COMM(simpar),99);
                    fprintf(ledger,"step\tt_start\th\tattempts\trhs_total\tmin_U\tEhyd\tEpot\tmass\tpx\tpy\tLz\tUint\thydro_integral\tgravity_work\thydro_abs_integral\tfinite_step_hydro_residual\tfinite_step_potential_residual\tfloor_injection\n");}
                fprintf(ledger,"%d\t%.17g\t%.17g\t%d\t%ld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t0\n",
                    gfs_rk4_context_step,t0,h,attempt+1,gfs_ssp_rhs_count,minU,e1[0],e1[1],e1[2],e1[3],e1[4],e1[5],e1[6],hyd,wg,hab,e1[0]-e0[0]-hyd-wg,e1[1]-e0[1]+wg);
                fflush(ledger);
                fprintf(stderr,"[SSP_ACCEPT] step=%d t=%.17g h=%.17g attempts=%d rhs=%ld minU=%.17g floors=0\n",gfs_rk4_context_step,t0,h,attempt+1,gfs_ssp_rhs_count,minU);
            }
            gpu_geometry_trial_end();free(backup);return h;
        }
        gpu_geometry_trial_end();
        gfs_ssp_restore(simpar,backup,np0,ntotal,dtold,floor0,nfloor0,pair0);
        ceiling=fmin(.5*h,retry_bound);
        if(MYID(simpar)==0) fprintf(stderr,"[SSP_ROLLBACK] step=%d attempt=%d restored_full_payload=1 dtold=%.17g next_ceiling=%.17g\n",gfs_rk4_context_step,attempt,(double)dtold,ceiling);
        if(!(ceiling>0) || !(t0+ceiling>t0)) break;
    }
budget:
    gpu_geometry_trial_end();gfs_rk4_context_time=t0;
    if(NID(simpar)==1)
        gfs_rk4_failure_snapshot(simpar,backup,(int)np0,"ssp_rejected_start",dtold,dtold);
    free(backup);
    fprintf(stderr,"[SSPRK3_STOP] budget/retry/time-progress limit rhs=%ld accepted=%ld; no invalid trial committed\n",gfs_ssp_rhs_count,gfs_ssp_accepted_count);
    MPI_Abort(MPI_COMM(simpar),99);return 0;
}
#endif

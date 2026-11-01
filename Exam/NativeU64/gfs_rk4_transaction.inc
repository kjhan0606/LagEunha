/* Default-off, smooth W0 RK4 trial acceptance. No energy/state repair.
 * Initially restricted to one MPI rank; the old MPI RK4 path is unchanged.
 * Restore uses the pre-existing SSP transaction's full-payload/cache reset.
 */
#ifdef USE_CUDA
static long gfs_rk4_tx_rhs_total=0;
typedef struct { size_t id; double U,V,m; int seen; double dlog,cast; } GfsRK4Base;
typedef struct {
    int on,attempt,diag; ptrdiff_t np,total; size_t bytes;
    void *backup; GfsRK4Base *base;
    float dtold; double floors; int nfloors; long long pairs;
    double ceiling,external,budget,span,factor,defect,quant,h;
    const char *reason;
} GfsRK4Transaction;
static int gfs_rk4_base_cmp(const void *a,const void *b){
    size_t x=((const GfsRK4Base*)a)->id,y=((const GfsRK4Base*)b)->id;
    return (x>y)-(x<y);
}
static double gfs_rk4_tx_parameter(const char *name,double fallback){
    const char *s=getenv(name);if(!s || !s[0])return fallback;
    char *end=NULL;double v=strtod(s,&end);
    if(*end || !isfinite(v)){fprintf(stderr,"[RK4_TX] invalid %s\n",name);MPI_Abort(MPI_COMM_WORLD,99);}
    return v;
}
static int gfs_rk4_tx_valid(SimParameters *simpar,GfsRK4Transaction *tx,int rhs){
    if(!tx->on)return 1;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    if(gpu_geometry_trial_failed()){tx->reason="geometry";return 0;}
    if(VORO_NP(simpar)!=tx->np){tx->reason="ownership_count";return 0;}
    for(ptrdiff_t k=0;k<tx->np;k++){
        double U=gfs_t64_U(p+k),V=gfs_t64_V(p+k);
        int tag_ok=!gfs_thermal64_requested() || GFS_T64_TAGGED(p[k]);
        if(!tag_ok || !(U>0) || !(V>0) || !(p[k].ie>0) || !(p[k].volume>0) || !(p[k].mass>0) || p[k].w2!=0 ||
           !isfinite(U) || !isfinite(V) || !isfinite(p[k].ie) || !isfinite(p[k].volume) || !isfinite(p[k].mass) ||
           !isfinite(p[k].x) || !isfinite(p[k].y) || !isfinite(p[k].vx) || !isfinite(p[k].vy) ||
           (rhs && (!isfinite(p[k].ax) || !isfinite(p[k].ay) || !isfinite(p[k].die) ||
             (gfs_thermal64_requested() && !isfinite(GFS_T64_RHS(p[k])))))){
            tx->reason="nonpositive_or_nonfinite";return 0;
        }
    }
    return 1;
}
static void gfs_rk4_tx_begin(SimParameters *simpar,GfsRK4Transaction *tx){
    memset(tx,0,sizeof(*tx));
    const char *s=getenv("GFS_RK4_TRANSACTION");tx->on=s && atoi(s)>0;
    if(!tx->on)return;
    const char *diag=getenv("GFS_RK4_TX_ENTROPY_DIAGNOSTIC");tx->diag=diag && atoi(diag)>0;
    gfs_es_config();
    if(NID(simpar)!=1 || !GAS_GPU_ENABLED(simpar) || !gpu_varpress_w0_enabled() ||
       GAS_AVMODE(simpar)!=5 || GAS_Kappa(simpar)!=0 || GAS_ENTROPY_MODE(simpar)!=0 ||
       !sedov_phase1_on() || GAS_XSPHEPS(simpar)!=0 || GAS_FCENTROID(simpar)!=0 ||
       GAS_VISCOSITY(simpar)!=0 || GAS_AlphaVis(simpar)!=0 || GAS_BetaVis(simpar)!=0 ||
       GAS_CDAMAX(simpar)!=0 || GAS_HYPERVALPHA(simpar)!=0 || gfs_dual_energy_eta()>0 || gfs_es_on || gfs_half_on ||
       gfs_dynw2_on() || gfs_w2lin_on() || sedov_lagvol_on() || gfs_vaud_on() ||
       !getenv("GFS_GPU_RAW_FACES") || atoi(getenv("GFS_GPU_RAW_FACES"))!=1 ||
       !getenv("GFS_GPU_FINAL_VOLUME") || atoi(getenv("GFS_GPU_FINAL_VOLUME"))!=1 ||
       !getenv("EUNHA_IC") || strcmp(getenv("EUNHA_IC"),"kepler") ||
       TVORORK4_DDINFO(simpar)[0].n_size!=sizeof(treevorostressrk4particletype)){
        fprintf(stderr,"[RK4_TX] unsupported contract: one-rank native smooth W0 only\n");MPI_Abort(MPI_COMM(simpar),99);
    }
    double c=gfs_rk4_tx_parameter("GFS_VARPRESS_GEOM_CFL",0);
    tx->budget=gfs_rk4_tx_parameter("GFS_RK4_TX_ENTROPY_BUDGET",.001);
    double ref=gfs_rk4_tx_parameter("GFS_RK4_TX_REF_TIME",-1);
    double end=gfs_rk4_tx_parameter("GFS_RK4_TX_BUDGET_END",10);
    tx->factor=gfs_rk4_tx_parameter("GFS_RK4_TX_H_FACTOR",1);
    tx->span=end-ref;
    if(c!=.2 || ref<0 || ref>gfs_rk4_context_time || !(end>gfs_rk4_context_time) ||
       !(tx->budget>0) || !(tx->factor>0 && tx->factor<=1)){
        fprintf(stderr,"[RK4_TX] invalid fixed C=.2/entropy budget/reference/step factor\n");MPI_Abort(MPI_COMM(simpar),99);
    }
    tx->np=VORO_NP(simpar);tx->total=VORO_TNP(simpar);
    tx->bytes=(size_t)tx->np*sizeof(treevorostressrk4particletype);
    tx->backup=malloc(tx->bytes?tx->bytes:1);
    tx->base=malloc((size_t)(tx->np?tx->np:1)*sizeof(*tx->base));
    if(!tx->backup || !tx->base)MPI_Abort(MPI_COMM(simpar),99);
    memcpy(tx->backup,VORORK4_TBP(simpar),tx->bytes);
    treevorostressrk4particletype *p=tx->backup;
    for(ptrdiff_t k=0;k<tx->np;k++)tx->base[k]=(GfsRK4Base){
        (size_t)PINDX((treevorork4particletype*)(p+k)),gfs_t64_U(p+k),gfs_t64_V(p+k),p[k].mass};
    qsort(tx->base,(size_t)tx->np,sizeof(*tx->base),gfs_rk4_base_cmp);
    for(ptrdiff_t k=1;k<tx->np;k++)if(tx->base[k-1].id==tx->base[k].id)MPI_Abort(MPI_COMM(simpar),99);
    tx->dtold=GAS_dtold(simpar);tx->floors=gfs_stage_floor_l;tx->nfloors=gfs_stage_floor_n_l;tx->pairs=gfs_pair_face_ends;
    tx->external=gfs_rk4_step_ceiling;tx->ceiling=tx->external>0?tx->external:1.e100;
    if(!gfs_rk4_tx_valid(simpar,tx,0)){
        fprintf(stderr,"[RK4_TX] invalid accepted start\n");MPI_Abort(MPI_COMM(simpar),99);
    }
}
static int gfs_rk4_tx_bound(SimParameters *simpar,GfsRK4Transaction *tx,double h,double bound){
    if(!tx->on)return 1;
    gfs_rk4_tx_rhs_total++;
    if(!gfs_rk4_tx_valid(simpar,tx,1))return 0;
    double orbit=kepler_dt_limit(simpar,VORORK4_TBP(simpar),sizeof(treevorostressrk4particletype),(int)VORO_NP(simpar));
    if(!(bound>0) || !isfinite(bound) || h>fmin(bound,orbit)*(1.+2.e-7)){
        tx->reason="later_stage_dt";return 0;
    }
    return 1;
}
static int gfs_rk4_tx_final(SimParameters *simpar,GfsRK4Transaction *tx,double h){
    if(!tx->on)return 1;
    if(!gfs_rk4_tx_valid(simpar,tx,0))return 0;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    tx->defect=tx->quant=0;
    for(ptrdiff_t k=0;k<tx->np;k++)tx->base[k].seen=0;
    for(ptrdiff_t k=0;k<tx->np;k++){
        GfsRK4Base key={.id=(size_t)PINDX((treevorork4particletype*)(p+k))};
        GfsRK4Base *b=bsearch(&key,tx->base,(size_t)tx->np,sizeof(*b),gfs_rk4_base_cmp);
        if(!b || b->seen || b->m!=p[k].mass){tx->reason="id_or_mass";return 0;}
        b->seen=1;
        double U=gfs_t64_U(p+k),V=gfs_t64_V(p+k),gm1=(double)GAS_GAMMA(simpar)-1;
        double d=log1p((U-b->U)/b->U)+gm1*log1p((V-b->V)/b->V);
        double rho=p[k].mass/V,pressure=gm1*U/V;
        if(!isfinite(d) || !(rho>0) || !(pressure>0) || !isfinite(rho) || !isfinite(pressure)){
            tx->reason="final_entropy_or_EOS_invalid";return 0;
        }
        b->dlog=d;
        b->cast=tx->diag?(gfs_thermal64_requested()?((double)p[k].ie-U)/U:
            (U-p[k].stress.E_inv_xy)/p[k].stress.E_inv_xy):0;
        if(!isfinite(b->cast)){tx->reason="nonfinite_final_cast";return 0;}
        tx->defect=fmax(tx->defect,fabs(d));
        /* Storage resolution only; never added to the acceptance tolerance. */
        double q=((double)nextafterf(p[k].ie,INFINITY)-U)/U+
            gm1*((double)nextafterf(p[k].volume,INFINITY)-V)/V;
        tx->quant=fmax(tx->quant,q);
    }
    static int forced_done=0;
    const char *force=getenv("GFS_RK4_TX_FORCE_REJECT_ONCE");
    if(force && atoi(force)>0 && !forced_done){forced_done=1;tx->reason="verification_injection";return 0;}
    if(!tx->diag && tx->defect>tx->budget*h/tx->span){tx->reason="final_entropy_defect";return 0;}
    tx->h=h;
    fprintf(stderr,"[RK4_TX_ACCEPT] step=%d t=%.17g h=%.17g attempts=%d rhs_total=%ld max_log_entropy_defect=%.9e allowed=%.9e storage_ulp=%.9e\n",
        gfs_rk4_context_step,gfs_rk4_context_time,h,tx->attempt+1,gfs_rk4_tx_rhs_total,tx->defect,tx->budget*h/tx->span,tx->quant);
    return 1;
}
static int gfs_rk4_tx_reject(SimParameters *simpar,GfsRK4Transaction *tx,double h){
    gpu_geometry_trial_end();
    gfs_ssp_restore(simpar,tx->backup,tx->np,tx->total,tx->dtold,tx->floors,tx->nfloors,tx->pairs);
    int exact=memcmp(VORORK4_TBP(simpar),tx->backup,tx->bytes)==0 &&
        VORO_NP(simpar)==tx->np && VORO_TNP(simpar)==tx->total && GAS_dtold(simpar)==tx->dtold &&
        gfs_stage_floor_l==tx->floors && gfs_stage_floor_n_l==tx->nfloors && gfs_pair_face_ends==tx->pairs;
    fprintf(stderr,"[RK4_TX_ROLLBACK] step=%d attempt=%d reason=%s h=%.17g payload_count_dtold_ledgers_exact=%d max_defect=%.9e storage_ulp=%.9e\n",
        gfs_rk4_context_step,tx->attempt+1,tx->reason?tx->reason:"invalid",h,exact,tx->defect,tx->quant);
    if(!exact)MPI_Abort(MPI_COMM(simpar),99);
    tx->attempt++;tx->ceiling=.5*h;
    if(tx->attempt>=8 || !(tx->ceiling>0) || !(gfs_rk4_context_time+tx->ceiling>gfs_rk4_context_time)){
        gfs_rk4_failure_snapshot(simpar,tx->backup,(int)tx->np,"tx_rejected_start",tx->dtold,tx->dtold);
        gfs_rk4_step_ceiling=tx->external;
        free(tx->backup);free(tx->base);
        fprintf(stderr,"[RK4_TX_STOP] eight-rejection/time-progress limit; accepted start preserved, no heat repair\n");
        MPI_Abort(MPI_COMM(simpar),99);return 0;
    }
    return 1;
}
typedef struct {size_t id;double signed_sum,abs_sum,max_cast;} GfsRK4EntropyLedger;
static void gfs_rk4_tx_entropy_commit(const GfsRK4Transaction *tx){
    if(!tx->diag)return;
    static GfsRK4EntropyLedger *ledger=NULL;static ptrdiff_t count=0;
    static double start_time=0;static long steps=0,rejections=0;
    if(!ledger){
        count=tx->np;ledger=calloc((size_t)(count?count:1),sizeof(*ledger));
        if(!ledger)MPI_Abort(MPI_COMM_WORLD,99);
        for(ptrdiff_t k=0;k<count;k++)ledger[k].id=tx->base[k].id;
        start_time=gfs_rk4_context_time;
    }
    if(count!=tx->np)MPI_Abort(MPI_COMM_WORLD,99);
    double signed_max=0,abs_max=0,cast_max=0;
    for(ptrdiff_t k=0;k<count;k++){
        if(ledger[k].id!=tx->base[k].id)MPI_Abort(MPI_COMM_WORLD,99);
        ledger[k].signed_sum+=tx->base[k].dlog;
        ledger[k].abs_sum+=fabs(tx->base[k].dlog);
        ledger[k].max_cast=fmax(ledger[k].max_cast,fabs(tx->base[k].cast));
        signed_max=fmax(signed_max,fabs(ledger[k].signed_sum));
        abs_max=fmax(abs_max,ledger[k].abs_sum);cast_max=fmax(cast_max,ledger[k].max_cast);
    }
    steps++;rejections+=tx->attempt;
    fprintf(stderr,"[RK4_TX_ENTROPY_DIAGNOSTIC] step=%d since=%.17g steps=%ld rejections=%ld max_abs_signed=%.9e max_sum_abs=%.9e max_final_cast_relative=%.9e budget=%.9e diagnostic_only=1\n",
        gfs_rk4_context_step,start_time,steps,rejections,signed_max,abs_max,cast_max,tx->budget);
    const char *target=getenv("HYDRO_TSTOP");
    if(target && fabs(gfs_rk4_context_time+tx->h-atof(target))<1.e-12){
        FILE *f=fopen("rk4_tx_entropy_cells.tsv","wx");if(!f)MPI_Abort(MPI_COMM_WORLD,99);
        fprintf(f,"id\tsum_signed_log_entropy\tsum_abs_log_entropy\tmax_final_cast_relative\n");
        for(ptrdiff_t k=0;k<count;k++)fprintf(f,"%zu\t%.17g\t%.17g\t%.17g\n",ledger[k].id,ledger[k].signed_sum,ledger[k].abs_sum,ledger[k].max_cast);
        if(fclose(f))MPI_Abort(MPI_COMM_WORLD,99);
    }
}
static void gfs_rk4_tx_commit(GfsRK4Transaction *tx){
    if(!tx->on)return;
    gfs_rk4_tx_entropy_commit(tx);
    gpu_geometry_trial_end();gfs_rk4_step_ceiling=tx->external;
    free(tx->backup);free(tx->base);
}
#else
typedef struct {int on;} GfsRK4Transaction;
#endif

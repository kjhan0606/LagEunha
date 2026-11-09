/* Bounded periodic smooth canonical (x,p,W,pi) KDK qualification.
 * No dissipation, mass transfer, clipping, or production default changes.
 * U(V) is the same adiabatic potential used by the qualified W0 KDK.
 */
typedef struct {unsigned long long id;ptrdiff_t index;} GfsWeightBaseIndex;
static int gfs_weight_base_order(const void *a,const void *b){
    unsigned long long x=((const GfsWeightBaseIndex*)a)->id,y=((const GfsWeightBaseIndex*)b)->id;
    return (x>y)-(x<y);
}
static ptrdiff_t gfs_weight_base_find(const GfsWeightBaseIndex *map,ptrdiff_t n,unsigned long long id){
    GfsWeightBaseIndex key={id,0};
    const GfsWeightBaseIndex *p=bsearch(&key,map,(size_t)n,sizeof(*map),gfs_weight_base_order);
    return p?p->index:-1;
}
static double gfs_native_weight_kdk(
    SimParameters *s,void (*padding)(SimParameters *,postype),
    Voro2D_point *(*neighbors)(SimParameters *,int,int,int *),
    treevorork4particletype *(*cells)(SimParameters *,int,int,int *),
    void (*linked)(SimParameters *,postype,postype,postype,postype,postype,
                   void (*)(SimParameters *,postype)),
    void (*post)(SimParameters *))
{
    const char *blocked[]={"GFS_UNIFIED_HYDRO","GFS_KEPLER_FIXED_POLYGON",
        "GFS_KEPLER_WALL_EVENTS","GFS_WEIGHT_NATIVE_RHS_ONLY","GFS_WALL_RHS_PROBE_ONLY",
        "GFS_RK4_TRANSACTION","GFS_THERMAL64","GFS_DYN_W2","GFS_W2_LINEAR"};
    int bad=0;
    for(size_t k=0;k<sizeof(blocked)/sizeof(blocked[0]);k++){
        const char *v=getenv(blocked[k]);if(v&&strcmp(v,"0"))bad=1;
    }
    const char *raw=getenv("GFS_GPU_RAW_FACES"),*fv=getenv("GFS_GPU_FINAL_VOLUME");
    const char *ic=getenv("EUNHA_IC");
    if(bad||NID(s)!=1||post!=periodic_postStage_blend||GAS_EVOLMETHOD(s)!=3||
       !GAS_GPU_ENABLED(s)||!gpu_varpress_w0_enabled()||
       !sedov_phase1_on()||GAS_AVMODE(s)!=5||GAS_Kappa(s)!=0||
       GAS_ENTROPY_MODE(s)!=0||GAS_VISCOSITY(s)!=0||GAS_AlphaVis(s)!=0||
       GAS_BetaVis(s)!=0||GAS_CDAMAX(s)!=0||GAS_HYPERVALPHA(s)!=0||
       GAS_FCENTROID(s)!=0||GAS_XSPHEPS(s)!=0||GAS_ACCX(s)!=0||GAS_ACCY(s)!=0||
       !raw||strcmp(raw,"1")||!fv||strcmp(fv,"1")||!ic||!strcmp(ic,"kepler")){
        fprintf(stderr,"[WEIGHT_KDK] unsupported smooth periodic validation contract\n");
        MPI_Abort(MPI_COMM(s),99);
    }
    const int active_weight=gfs_weight_native_active();
    const ptrdiff_t n=VORO_NP(s),total=VORO_TNP(s);
    size_t bytes=(size_t)n*sizeof(treevorostressrk4particletype);
    treevorostressrk4particletype *base=malloc(bytes);
    double *startV=malloc((size_t)n*sizeof(double));
    if(!base||!startV||n<1)MPI_Abort(MPI_COMM(s),99);
    memcpy(base,VORORK4_TBP(s),bytes);
    GfsWeightBaseIndex *map=malloc((size_t)n*sizeof(*map));
    if(!map)MPI_Abort(MPI_COMM(s),99);
    for(ptrdiff_t i=0;i<n;i++){
        map[i]=(GfsWeightBaseIndex){(unsigned long long)PINDX(base+i),i};
        if(!gfs_weight_particle_valid(base+i,
            GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_THERMAL64))
            MPI_Abort(MPI_COMM(s),99);
    }
    qsort(map,(size_t)n,sizeof(*map),gfs_weight_base_order);
    for(ptrdiff_t i=1;i<n;i++)if(map[i].id==map[i-1].id)MPI_Abort(MPI_COMM(s),99);
    const float dtold=GAS_dtold(s);const double floors=gfs_stage_floor_l;
    const int nfloors=gfs_stage_floor_n_l;const long long pairs=gfs_pair_face_ends;
    double gamma=GAS_GAMMA(s),cell=BASICCELL_CELLWIDTH(s),halo=gfs_kepler_halo_width(cell);
    double xmin=Xmin_HydroExam(s)-halo,ymin=Ymin_HydroExam(s)-halo;
    double xmax=Xmax_HydroExam(s)+halo,ymax=Ymax_HydroExam(s)+halo;
    double ceiling=gfs_rk4_step_ceiling>0?gfs_rk4_step_ceiling:.02;
    /* Use the first admissible CFL as a dyadic root for this continuous run.
     * Quantization only shortens steps; the later-stage bound and rollback
     * remain mandatory. This avoids rediscovering an almost-identical rejected
     * upper step after every accepted half step. It is global, not local LTS. */
    static SimParameters *bin_sim=NULL;
    static double bin_root=0,bin_last_time=-1;
    if(bin_sim!=s || bin_last_time!=gfs_rk4_context_time)bin_root=0;
    double h=0,begin=MPI_Wtime();int attempt;
    for(attempt=0;attempt<12;attempt++){
        int reject=0;
        gpu_geometry_trial_begin();
        updateDenW2Pressure2DBlend(s,xmin,ymin,xmax,ymax,gamma,padding,neighbors,cells,linked,1);
        if(gpu_geometry_trial_failed()){reject=1;goto retry;}
        h=getAccVoro2DBlend_GPU(s,xmin,ymin,xmax,ymax,VoroAccuracyOrder(s),
            GAS_COURANT(s),gamma,padding,neighbors,cells,linked);
        if(!(h>0)||!isfinite(h)){reject=1;goto retry;}
        if(!(bin_root>0))bin_root=h;
        double binned=bin_root;
        while(binned>h)binned*=.5;
        h=fmin(binned,ceiling);
        if(!(h>0)||!isfinite(h)){reject=1;goto retry;}
        {
            treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
            for(ptrdiff_t i=0;i<n;i++){
                ptrdiff_t j=gfs_weight_base_find(map,n,(unsigned long long)PINDX(p+i));
                if(j<0){reject=1;break;}
                /* Retain the geometric start volume, not a rounded IC proxy. */
                startV[j]=p[i].weight_thermal.V;
                p[i].vx+=.5*h*p[i].ax;p[i].vy+=.5*h*p[i].ay;
                if(active_weight)p[i].weight.pi+=.5*h*p[i].weight_rk4.k4pi;
                p[i].x+=h*p[i].vx;p[i].y+=h*p[i].vy;
                if(active_weight)p[i].weight.W+=h*p[i].weight.pi/p[i].weight.inertia;
                if(!isfinite(p[i].x)||!isfinite(p[i].y)||!isfinite(p[i].weight.W))reject=1;
            }
        }
        if(reject)goto retry;
        post(s);
        if(VORO_NP(s)!=n){reject=1;goto retry;}
        updateDenW2Pressure2DBlend(s,xmin,ymin,xmax,ymax,gamma,padding,neighbors,cells,linked,h);
        if(gpu_geometry_trial_failed()){reject=1;goto retry;}
        {
            treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
            for(ptrdiff_t i=0;i<n;i++){
                double U;
                ptrdiff_t j=gfs_weight_base_find(map,n,(unsigned long long)PINDX(p+i));
                if(j<0||p[i].mass!=base[j].mass||
                   !gfs_variational_adiabatic(base[j].weight_thermal.U,startV[j],
                                              p[i].weight_thermal.V,gamma,&U)){reject=1;break;}
                p[i].weight_thermal.U=U;p[i].ie=(float)U;
                double V=p[i].weight_thermal.V,rho=p[i].mass/V,P=(gamma-1)*U/V;
                p[i].volume=(float)V;p[i].den=(float)rho;p[i].pressure=(float)P;
                p[i].csound=(float)sqrt(gamma*P/rho);
            }
        }
        if(reject)goto retry;
        {
            double bound=getAccVoro2DBlend_GPU(s,xmin,ymin,xmax,ymax,VoroAccuracyOrder(s),
                GAS_COURANT(s),gamma,padding,neighbors,cells,linked);
            if(!isfinite(bound)||!(bound>0)||h>bound*(1+2e-7)){reject=1;goto retry;}
            treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
            long double E=0,E0=0,M=0;double minV=HUGE_VAL,minU=HUGE_VAL;
            for(ptrdiff_t i=0;i<n;i++){
                p[i].vx+=.5*h*p[i].ax;p[i].vy+=.5*h*p[i].ay;
                if(active_weight)p[i].weight.pi+=.5*h*p[i].weight_rk4.k4pi;
                ptrdiff_t j=gfs_weight_base_find(map,n,(unsigned long long)PINDX(p+i));
                if(j<0||p[i].mass!=base[j].mass||
                   !gfs_weight_particle_valid(p+i,GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_THERMAL64)){
                    reject=1;break;
                }
                E+=p[i].weight_thermal.U+.5L*p[i].mass*(p[i].vx*p[i].vx+p[i].vy*p[i].vy)+gfs_weight_energy(&p[i].weight);
                E0+=base[i].weight_thermal.U+.5L*base[i].mass*(base[i].vx*base[i].vx+base[i].vy*base[i].vy)+gfs_weight_energy(&base[i].weight);
                M+=p[i].mass;minV=fmin(minV,p[i].weight_thermal.V);minU=fmin(minU,p[i].weight_thermal.U);
            }
            if(reject)goto retry;
            fprintf(stderr,"[WEIGHT_KDK_ACCEPT] step=%d t=%.17g dt=%.17g E=%.17Lg step_dE=%.17Lg M=%.17Lg Vmin=%.17g Umin=%.17g retries=%d wall_s=%.6f\n",
                gfs_rk4_context_step,gfs_rk4_context_time+h,h,E,E-E0,M,minV,minU,attempt,MPI_Wtime()-begin);
        }
        bin_sim=s;bin_last_time=gfs_rk4_context_time+h;
        gpu_geometry_trial_end();free(map);free(startV);free(base);return h;
retry:
        gpu_geometry_trial_end();
        gfs_ssp_restore(s,base,n,total,dtold,floors,nfloors,pairs);
        ceiling=.5*(h>0?h:ceiling);
        fprintf(stderr,"[WEIGHT_KDK_RETRY] attempt=%d ceiling=%.17g\n",attempt+1,ceiling);
    }
    /* Rollback above restores the last accepted step-start state. Preserve it
     * for a diagnostic replay; never serialize the rejected trial as accepted. */
    gfs_rk4_failure_snapshot(s,base,n,"weight_kdk_stop",0,dtold);
    fprintf(stderr,"[WEIGHT_KDK_STOP] retry budget exhausted; start state restored and diagnostic snapshot saved\n");
    free(map);free(startV);free(base);MPI_Abort(MPI_COMM(s),99);return 0;
}

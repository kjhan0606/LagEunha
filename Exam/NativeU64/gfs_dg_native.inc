/* Opt-in fresh Euler projection only; no polynomial timestep yet. */
#ifdef USE_CUDA
static void gfs_dg_gauss_rule(int n,double *node,double *weight){
    for(int i=0;i<(n+1)/2;i++){
        long double z=cosl(acosl(-1.L)*(i+.75L)/(n+.5L)),dp=0;
        for(int iteration=0;iteration<100;iteration++){
            long double p0=1,p1=z;
            for(int k=2;k<=n;k++){long double q=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=q;}
            dp=n*(z*p1-p0)/(z*z-1);
            long double next=z-p1/dp;
            if(fabsl(next-z)<4*LDBL_EPSILON){z=next;break;}z=next;
        }
        long double p0=1,p1=z;
        for(int k=2;k<=n;k++){long double q=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=q;}
        dp=n*(z*p1-p0)/(z*z-1);
        double w=1/((1-z*z)*dp*dp);
        node[i]=(1-z)*.5L;node[n-1-i]=(1+z)*.5L;weight[i]=weight[n-1-i]=w;
    }
}

#include "gfs_enriched_native.inc"

static void gfs_dg_projection_gate(SimParameters *simpar,
    void (*padding)(SimParameters*,postype),
    Voro2D_point *(*neighbors)(SimParameters*,int,int,int*),
    treevorork4particletype *(*cell)(SimParameters*,int,int,int*),
    void (*linked)(SimParameters*,postype,postype,postype,postype,postype,void(*)(SimParameters*,postype)),
    void (*post)(SimParameters*)){
    const char *raw=getenv("GFS_GPU_RAW_FACES"),*fv=getenv("GFS_GPU_FINAL_VOLUME"),*profile=getenv("EUNHA_KEPLER_PROFILE");
    const char *icname=getenv("EUNHA_IC"),*aff=getenv("GFS_DG_AFFINE");
    if(!GAS_GPU_ENABLED(simpar) || GAS_AVMODE(simpar)!=5 || GAS_Kappa(simpar)!=0 ||
       GAS_FCENTROID(simpar)!=0 || GAS_XSPHEPS(simpar)!=0 || GAS_VISCOSITY(simpar)!=0 ||
       GAS_AlphaVis(simpar)!=0 || GAS_BetaVis(simpar)!=0 || GAS_CDAMAX(simpar)!=0 ||
       GAS_HYPERVALPHA(simpar)!=0 || GAS_ENTROPY_MODE(simpar)!=0 || gfs_rk4_context_time!=0 ||
       !raw || atoi(raw)<=0 || !fv || atoi(fv)<=0 || post!=periodic_postStage_blend ||
       (profile && profile[0]) || !icname || (strcmp(icname,"kepler") && strcmp(icname,"affine_probe")) ||
       TVORORK4_DDINFO(simpar)[0].n_size!=sizeof(treevorostressrk4particletype)){
        fprintf(stderr,"[DG_PROJECT] unsupported fresh native geometry/IC contract\n");MPI_Abort(MPI_COMM(simpar),99);
    }
    GfsDgConfig cfg={0};cfg.order=(int)gfs_moment_env("GFS_DG_ORDER",2);
    if(cfg.order!=2 && cfg.order!=4 && cfg.order!=6)MPI_Abort(MPI_COMM(simpar),99);
    cfg.ic.gamma=GAS_GAMMA(simpar);cfg.ic.pressure=1.e-6;cfg.ic.kind=aff && aff[0]?1:2;
    if(cfg.ic.kind==1){
        double v[10];const char *s=aff;
        for(int k=0;k<10;k++){
            char *end=NULL;v[k]=strtod(s,&end);
            if(end==s || !isfinite(v[k]) || (k<9?*end!=',':*end!=0))MPI_Abort(MPI_COMM(simpar),99);
            s=end+(k<9);
        }
        cfg.ic.rho=v[0];cfg.ic.pressure=v[1];cfg.ic.cx=v[2];cfg.ic.cy=v[3];cfg.ic.vx=v[4];cfg.ic.vy=v[5];
        for(int k=0;k<4;k++)cfg.ic.B[k]=v[6+k];
        const char *slope=getenv("GFS_DG_DENSITY_SLOPE");
        if(slope && slope[0]){
            char *end=NULL;cfg.density_slope[0]=strtod(slope,&end);
            if(end==slope || *end!=',')MPI_Abort(MPI_COMM(simpar),99);
            const char *second=end+1;cfg.density_slope[1]=strtod(second,&end);
            if(end==second || *end || !isfinite(cfg.density_slope[0]) || !isfinite(cfg.density_slope[1]))MPI_Abort(MPI_COMM(simpar),99);
        }
    }else{
        cfg.ic.cx=gfs_moment_env("EUNHA_KEPLER_CX",.5*SIMBOX(simpar).x.max);
        cfg.ic.cy=gfs_moment_env("EUNHA_KEPLER_CY",.5*SIMBOX(simpar).y.max);
        cfg.ic.eps=gfs_moment_env("EUNHA_KEPLER_EPS",0);
        cfg.ic.rin=gfs_moment_env("EUNHA_KEPLER_RIN",2);cfg.ic.rout=gfs_moment_env("EUNHA_KEPLER_ROUT",8);
        if(!(cfg.ic.eps>0))MPI_Abort(MPI_COMM(simpar),99);
    }
    cfg.ng[0]=cfg.order+2;cfg.ng[1]=cfg.order+4;cfg.ng[2]=cfg.order+3;
    for(int k=0;k<3;k++)gfs_dg_gauss_rule(cfg.ng[k],cfg.node[k],cfg.weight[k]);
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    for(ptrdiff_t k=0;k<VORO_NP(simpar);k++)if(p[k].w2!=0 || gfs_moment_tagged(p+k))MPI_Abort(MPI_COMM(simpar),99);
    exam2dUpdateVol(simpar,padding,neighbors,cell,linked);
    VORORK4_TBPP(simpar)=NULL;VORO_BASICCELL(simpar)=NULL;
    int np=VORO_NP(simpar),bad=0,any=0;
    const char *enriched=getenv("GFS_DG_ENRICHED");
    if(enriched && atoi(enriched)>0){
        if(cfg.order!=2)MPI_Abort(MPI_COMM(simpar),99);
        gfs_enriched_projection_output(simpar,np,&cfg);
    }
    GfsDgProjectionAudit *audit=calloc((size_t)(np?np:1),sizeof(*audit));if(!audit)MPI_Abort(MPI_COMM(simpar),99);
    bad=gpu_dg_project_initial(gpu_get_context(),np,MYID(simpar),&cfg,audit);
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_SUM,MPI_COMM(simpar));
    int n=gfs_dg_nb(cfg.order);unsigned long long record_size=288+40*n;
    char path[160];snprintf(path,sizeof(path),"dg_project_r%d.bin",MYID(simpar));
    FILE *out=fopen(path,"wx");if(!out)MPI_Abort(MPI_COMM(simpar),99);
    unsigned long long header[8]={0x4746534447503031ULL,1,(unsigned long long)np,(unsigned long long)cfg.order,(unsigned long long)n,DG_NMETRICS,record_size,
        (unsigned long long)(cfg.ng[0]+16*cfg.ng[1]+256*cfg.ng[2])};
    fwrite(header,sizeof(header),1,out);
    int nonpositive[2]={0,0},total_nonpositive[2];
    for(int i=0;i<np;i++){
        unsigned long long id=(size_t)PINDX((treevorork4particletype*)(p+i));double xy[2]={p[i].x,p[i].y};
        fwrite(&id,sizeof(id),1,out);fwrite(xy,sizeof(xy),1,out);
        fwrite(&audit[i].basis.area,sizeof(double),6,out);fwrite(audit[i].basis.mean,sizeof(double),n,out);
        for(int c=0;c<4;c++)fwrite(audit[i].polynomial.c[c],sizeof(double),n,out);
        fwrite(audit[i].integral,sizeof(double),4,out);fwrite(audit[i].coarse_integral,sizeof(double),4,out);
        fwrite(audit[i].metrics,sizeof(double),DG_NMETRICS,out);
        int status[4]={audit[i].status,audit[i].nonpositive_rho,audit[i].nonpositive_p,audit[i].worst_kind};fwrite(status,sizeof(status),1,out);
        if(audit[i].nonpositive_rho)nonpositive[0]++;if(audit[i].nonpositive_p)nonpositive[1]++;
    }
    if(ferror(out) || fclose(out))MPI_Abort(MPI_COMM(simpar),99);
    free(audit);MPI_Allreduce(nonpositive,total_nonpositive,2,MPI_INT,MPI_SUM,MPI_COMM(simpar));
    if(any){
        fprintf(stderr,"[DG_PROJECT] projection/geometry failed cells=%d\n",any);
        /* The binary audit already preserves every per-cell failure code. */
        MPI_Abort(MPI_COMM(simpar),99);
    }
    if(MYID(simpar)==0)fprintf(stderr,"[DG_PROJECT_ONLY] order=%d completed; nonpositive_rho_cells=%d nonpositive_pressure_cells=%d; no evolution, clipping or physics pass implied\n",
        cfg.order,total_nonpositive[0],total_nonpositive[1]);
    MPI_Barrier(MPI_COMM(simpar));MPI_Finalize();exit(0);
}
#endif

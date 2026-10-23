/* Included by exam.c, after geometry and RK diagnostics. Defaults OFF. */
#include "gfs_moment_initial.h"
#ifdef USE_CUDA
typedef char gfs_moment_payload_fits[(sizeof(GfsMomentPayload)==24*sizeof(double) &&
    offsetof(RK4,w2backup)==24*sizeof(double) && sizeof(PosType)==sizeof(double))?1:-1];
#endif

static int gfs_moment_tagged(const treevorostressrk4particletype *p){
    if(sizeof(PosType)!=sizeof(double) || offsetof(RK4,w2backup)!=24*sizeof(double))return 0;
    GfsMomentPayload q;memcpy(&q,&p->rk4,sizeof(q));return q.tag==GFS_MOMENT_TAG;
}

#ifdef USE_CUDA
static double gfs_moment_env(const char *name,double value){
    const char *s=getenv(name);if(!s || !s[0])return value;
    char *end=NULL;double v=strtod(s,&end);
    if(end==s || *end || !isfinite(v)){fprintf(stderr,"[MOMENT_INIT] invalid %s\n",name);MPI_Abort(MPI_COMM_WORLD,99);}
    return v;
}

static void gfs_moment_initial_gate(SimParameters *simpar,
    void (*padding)(SimParameters*,postype),
    Voro2D_point *(*neighbors)(SimParameters*,int,int,int*),
    treevorork4particletype *(*cell)(SimParameters*,int,int,int*),
    void (*linked)(SimParameters*,postype,postype,postype,postype,postype,void(*)(SimParameters*,postype)),
    void (*post)(SimParameters*)){
    const char *raw=getenv("GFS_GPU_RAW_FACES"),*fv=getenv("GFS_GPU_FINAL_VOLUME");
    const char *profile=getenv("EUNHA_KEPLER_PROFILE"),*icname=getenv("EUNHA_IC");
    const char *aff=getenv("GFS_MOMENT_AFFINE");
    if(!GAS_GPU_ENABLED(simpar) || GAS_AVMODE(simpar)!=5 || GAS_Kappa(simpar)!=0 ||
       GAS_FCENTROID(simpar)!=0 || GAS_XSPHEPS(simpar)!=0 || GAS_VISCOSITY(simpar)!=0 ||
       GAS_AlphaVis(simpar)!=0 || GAS_BetaVis(simpar)!=0 || GAS_CDAMAX(simpar)!=0 ||
       GAS_HYPERVALPHA(simpar)!=0 || GAS_ENTROPY_MODE(simpar)!=0 ||
       !raw || atoi(raw)<=0 || !fv || atoi(fv)<=0 || post!=periodic_postStage_blend ||
       (profile && profile[0]) || !icname ||
       (strcmp(icname,"kepler") && strcmp(icname,"affine_probe")) ||
       TVORORK4_DDINFO(simpar)[0].n_size!=sizeof(treevorostressrk4particletype)){
        fprintf(stderr,"[MOMENT_INIT] unsupported native geometry/IC contract\n");MPI_Abort(MPI_COMM(simpar),99);
    }
    GfsMomentIC ic={0};ic.gamma=GAS_GAMMA(simpar);ic.pressure=1.e-6;
    ic.kind=aff && aff[0]?1:2;
    if(ic.kind==1){
        double values[10];const char *s=aff;
        for(int k=0;k<10;k++){
            char *end=NULL;values[k]=strtod(s,&end);
            if(end==s || !isfinite(values[k]) || (k<9?*end!=',':*end!=0))MPI_Abort(MPI_COMM(simpar),99);
            s=end+(k<9);
        }
        ic.rho=values[0];ic.pressure=values[1];ic.cx=values[2];ic.cy=values[3];
        ic.vx=values[4];ic.vy=values[5];for(int k=0;k<4;k++)ic.B[k]=values[6+k];
    }else{
        ic.cx=gfs_moment_env("EUNHA_KEPLER_CX",.5*SIMBOX(simpar).x.max);
        ic.cy=gfs_moment_env("EUNHA_KEPLER_CY",.5*SIMBOX(simpar).y.max);
        ic.eps=gfs_moment_env("EUNHA_KEPLER_EPS",0);
        ic.rin=gfs_moment_env("EUNHA_KEPLER_RIN",2);ic.rout=gfs_moment_env("EUNHA_KEPLER_ROUT",8);
        if(!(ic.eps>0) || !(ic.rin>0) || !(ic.rout>ic.rin))MPI_Abort(MPI_COMM(simpar),99);
    }
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    int bad=0,any=0,tagged=0;for(ptrdiff_t k=0;k<VORO_NP(simpar);k++){
        if(p[k].w2!=0)bad++;
        if(gfs_moment_tagged(p+k))tagged++;
        else if(gfs_rk4_context_time!=0)bad++;
    }
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_SUM,MPI_COMM(simpar));if(any)MPI_Abort(MPI_COMM(simpar),99);
    int total_tagged=0;
    MPI_Allreduce(&tagged,&total_tagged,1,MPI_INT,MPI_SUM,MPI_COMM(simpar));
    if(total_tagged!=0 && total_tagged!=VORO_TNP(simpar))MPI_Abort(MPI_COMM(simpar),99);
    /* A fresh IC is integrated only at t=0. A moment checkpoint is read as
     * state, never reinitialized from the analytic Kepler profile. */
    exam2dUpdateVol(simpar,padding,neighbors,cell,linked);
    VORORK4_TBPP(simpar)=NULL;VORO_BASICCELL(simpar)=NULL;
    int np=VORO_NP(simpar);
    GfsMomentPayload *input=calloc((size_t)(np?np:1),sizeof(*input));
    GfsMomentInitialAudit *audit=calloc((size_t)(np?np:1),sizeof(*audit));
    if(!input || !audit)MPI_Abort(MPI_COMM(simpar),99);
    p=(void*)VORORK4_TBP(simpar);
    for(int k=0;k<np;k++)if(gfs_moment_tagged(p+k))memcpy(input+k,&p[k].rk4,sizeof(*input));
    bad=gpu_moment_initialize(gpu_get_context(),np,&ic,input,audit);
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_SUM,MPI_COMM(simpar));
    char path[160];snprintf(path,sizeof(path),"moment_initial_r%d.bin",MYID(simpar));
    FILE *out=fopen(path,"wx");if(!out)MPI_Abort(MPI_COMM(simpar),99);
    unsigned long long header[6]={0x4746534d4f4d3031ULL,1,(unsigned long long)np,sizeof(*audit),
        offsetof(treevorostressrk4particletype,rk4),sizeof(GfsMomentPayload)};
    fwrite(header,sizeof(header),1,out);
    for(int k=0;k<np;k++){
        unsigned long long id=(size_t)PINDX((treevorork4particletype*)(p+k));
        double xy[2]={p[k].x,p[k].y};
        fwrite(&id,sizeof(id),1,out);fwrite(xy,sizeof(xy),1,out);fwrite(audit+k,sizeof(*audit),1,out);
    }
    if(ferror(out) || fclose(out))MPI_Abort(MPI_COMM(simpar),99);
    if(any){
        int shown=0;for(int k=0;k<np;k++)if(audit[k].valid!=1 && shown++ <8)
            fprintf(stderr,"[MOMENT_INIT_INVALID] rank=%d id=%zu status=%d residual=%.17g iterations=%d\n",
                MYID(simpar),(size_t)PINDX((treevorork4particletype*)(p+k)),audit[k].valid,audit[k].moment_residual,audit[k].iterations);
        MPI_Abort(MPI_COMM(simpar),99);
    }
    for(int k=0;k<np;k++){
        memcpy(&p[k].rk4,&audit[k].payload,sizeof(*input));
        GfsMomentPrimitive r;gfs_moment_primitive(&audit[k].payload.s,&r);
        p[k].mass=audit[k].payload.s.a[GFM_M];p[k].ie=r.U;
        p[k].vx=r.v[0]-r.B[0]*r.c[0]-r.B[1]*r.c[1];
        p[k].vy=r.v[1]-r.B[2]*r.c[0]-r.B[3]*r.c[1];
        p[k].den=p[k].mass/p[k].volume;
        p[k].pressure=(ic.gamma-1)*r.U/p[k].volume;
        p[k].csound=sqrt(ic.gamma*p[k].pressure/p[k].den);
    }
    /* Unique caller-owned directory; step+1 does not overwrite the input
     * restart. The existing native KH writer preserves the entire payload. */
    extern void kh_outdata(SimParameters*,int,postype,postype);
    if(MYID(simpar)==0){
        snprintf(path,sizeof(path),"khout.%06d.dat",gfs_rk4_context_step+1);
        FILE *existing=fopen(path,"rb");
        if(existing){fclose(existing);fprintf(stderr,"[MOMENT_INIT] checkpoint already exists: %s\n",path);MPI_Abort(MPI_COMM(simpar),99);}
    }
    kh_outdata(simpar,gfs_rk4_context_step+1,gfs_rk4_context_time,0);
    free(input);free(audit);
    if(MYID(simpar)==0)fprintf(stderr,"[MOMENT_INIT_ONLY] native CUDA quadrature/density fit completed; full-stride moment checkpoint exported; no evolution, weight solve or t6-10 validation\n");
    MPI_Barrier(MPI_COMM(simpar));MPI_Finalize();exit(0);
}
#endif

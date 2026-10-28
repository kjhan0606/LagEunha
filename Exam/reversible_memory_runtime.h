/* Included in exam.c after retry/configuration helpers. Experimental MPI1
 * periodic host implementation; row kernel is shared with the CUDA candidate.
 * Dedicated ID state and acceleration keep thermal work unmodified. */
#include "../Voro/Laguerre/reversible_memory.h"
#include "../Voro/Laguerre/reversible_memory_checkpoint.h"
typedef struct {int id;double dx,dy;} RmemCandidate;
static rm_owner_state *rmem_state;
static treevorostressrk4particletype **rmem_particle;
static double *rmem_acc,*rmem_rate,*rmem_columns;
static int *rmem_head,*rmem_next;
static long rmem_n;
static double rmem_power,rmem_trial_work,rmem_aux_before;
static int rmem_on(void){const char *s=getenv("GFS_REVERSIBLE_MEMORY");return s && atoi(s)>0;}
static int rmem_corotate(void){const char *s=getenv("GFS_RMEM_COROTATE");return s && atoi(s)>0;}
static void rmem_fail(SimParameters *s,const char *why){fprintf(stderr,"[RMEM] %s\n",why);MPI_Abort(MPI_COMM(s),2);}
static void rmem_map(SimParameters *s){
    if(NID(s)!=1 || TVORORK4_DDINFO(s)[0].n_size!=sizeof(treevorostressrk4particletype))rmem_fail(s,"requires MPI1 stress-sized particles");
    long n=VORO_NP(s);
    if(n<=0 || n>10000000)rmem_fail(s,"invalid particle count");
    if(!rmem_n){
        rmem_n=n;rmem_state=calloc(n,sizeof(*rmem_state));rmem_particle=calloc(n,sizeof(*rmem_particle));
        rmem_acc=calloc(2*n,sizeof(double));rmem_rate=calloc(2*n,sizeof(double));
        rmem_columns=calloc(n,sizeof(double));rmem_next=calloc(n,sizeof(int));
        if(!rmem_state || !rmem_particle || !rmem_acc || !rmem_rate || !rmem_columns || !rmem_next)rmem_fail(s,"allocation failed");
        for(long i=0;i<n;i++)rms_init(rmem_state+i,i,0,0);
    }
    if(n!=rmem_n)rmem_fail(s,"particle count changed");
    memset(rmem_particle,0,n*sizeof(*rmem_particle));
    treevorostressrk4particletype *p=(treevorostressrk4particletype*)VORORK4_TBP(s);
    for(long j=0;j<n;j++){
        long id=(long)PINDX(p+j);
        if(id<0 || id>=n || rmem_particle[id])rmem_fail(s,"requires unique dense persistent IDs");
        rmem_particle[id]=p+j;
    }
}
static uint64_t rmem_hydro_hash(void){
    uint64_t h=UINT64_C(14695981039346656037);
    for(long i=0;i<rmem_n;i++){
        treevorostressrk4particletype *p=rmem_particle[i];
        uint64_t id=i;double a[9]={p->x,p->y,p->vx,p->vy,p->mass,p->ie,p->w2,p->volume,p->stress.K};
        h=rms_hash_bytes(h,&id,sizeof(id));h=rms_hash_bytes(h,a,sizeof(a));
    }
    return h;
}
void gfs_rmem_snapshot_read(SimParameters *s,const char *path,double time){
    char side[4096];if(snprintf(side,sizeof(side),"%s.rmem",path)>=(int)sizeof(side))rmem_fail(s,"snapshot path too long");
    FILE *f=fopen(side,"rb");int exists=f!=NULL;if(f)fclose(f);
    if(!rmem_on()){if(exists)rmem_fail(s,"auxiliary snapshot requires GFS_REVERSIBLE_MEMORY=1");return;}
    rmem_map(s);
    if(exists){if(rms_checkpoint_read(side,rmem_state,rmem_n,rmem_hydro_hash(),time))rmem_fail(s,"sidecar validation failed");}
    else {
        const char *init=getenv("GFS_RMEM_INIT_ZERO");
        if(time!=0 && !(init && atoi(init)>0))rmem_fail(s,"missing sidecar: explicit zero-state branch required");
        fprintf(stderr,"[RMEM_INIT] zero auxiliary state from %s time=%.17g; new model branch\n",path,time);
    }
}
void gfs_rmem_snapshot_write(SimParameters *s,const char *path,double time){
    if(!rmem_on())return;
    if(!rmem_n && time!=0)rmem_fail(s,"uninitialized auxiliary snapshot");
    rmem_map(s);char side[4096];
    if(snprintf(side,sizeof(side),"%s.rmem",path)>=(int)sizeof(side))rmem_fail(s,"snapshot path too long");
    if(rms_checkpoint_write(side,rmem_state,rmem_n,rmem_hydro_hash(),time))rmem_fail(s,"sidecar write failed");
}
static double rmem_energy(void){double e=0;for(long i=0;i<rmem_n;i++)e+=rms_energy(rmem_state+i);return e;}
static int rmem_cmp(const void *a,const void *b){int i=((const RmemCandidate*)a)->id,j=((const RmemCandidate*)b)->id;return (i>j)-(i<j);}
static int rmem_mod(int a,int n){int r=a%n;return r<0?r+n:r;}
static double rmem_rhs(SimParameters *s){
    if(!rmem_on())return 1e300;
    rmem_map(s);
    double L=KH_SIMBOX(s).x.max-KH_SIMBOX(s).x.min;
    double Ly=KH_SIMBOX(s).y.max-KH_SIMBOX(s).y.min;
    if(!(L>0) || L!=Ly || KH_SIMBOX(s).x.min!=0 || KH_SIMBOX(s).y.min!=0)rmem_fail(s,"requires zero-origin square periodic box");
    int bins=(int)(sqrt((double)rmem_n)/4);if(bins<4)bins=4;
    if(!rmem_head){rmem_head=malloc(bins*bins*sizeof(int));if(!rmem_head)rmem_fail(s,"grid allocation failed");}
    for(int b=0;b<bins*bins;b++)rmem_head[b]=-1;
    double width=L/bins;
    for(long i=0;i<rmem_n;i++){
        treevorostressrk4particletype *p=rmem_particle[i];
        if(!isfinite(p->x) || !isfinite(p->y) || !(p->mass>0) || !(p->volume>0) || !(p->ie>0))rmem_fail(s,"invalid stage state");
        int ix=rmem_mod((int)floor(p->x/width),bins),iy=rmem_mod((int)floor(p->y/width),bins),b=ix+bins*iy;
        rmem_next[i]=rmem_head[b];rmem_head[b]=i;
    }
    memset(rmem_acc,0,2*rmem_n*sizeof(double));memset(rmem_rate,0,2*rmem_n*sizeof(double));
    memset(rmem_columns,0,rmem_n*sizeof(double));rmem_power=0;
    double max_row=0,max_spin=0;int row_errors=0;
    int corotate=rmem_corotate();
    const char *parallel_env=getenv("GFS_RMEM_PARALLEL");
    int parallel=parallel_env && atoi(parallel_env)>0;
    /* Initialize the frozen-boundary helper's static configuration on the
     * calling thread. Rows are independent; ordered scatter preserves the
     * original owner order, including the timestep column sums. */
    (void)gfs_kepler_frozen(rmem_particle[0]->x,rmem_particle[0]->y);
    #pragma omp parallel for if(parallel) ordered schedule(dynamic,1)
    for(long owner=0;owner<rmem_n;owner++){
        RmemCandidate candidates[4096];double dx[4096],dy[4096],row[4096],work[7*4096],meta[2];
        double gradient_x[4096],gradient_y[4096];
        int row_failure=0;
        treevorostressrk4particletype *p=rmem_particle[owner];
        if(gfs_kepler_frozen(p->x,p->y))continue;
        double reach=fmin(2*width,.49*L);int count=0,slot=-1,code=2;
        int ix=rmem_mod((int)floor(p->x/width),bins),iy=rmem_mod((int)floor(p->y/width),bins);
        while(code==2){
            if(!(reach>0 && reach<L/2)){row_failure=1;goto row_ready;}
            int span=(int)ceil(reach/width),full=2*span+1>=bins;count=0;
            int lo=full?0:-span,hi=full?bins-1:span;
            for(int ky=lo;ky<=hi;ky++)for(int kx=lo;kx<=hi;kx++){
                int bx=full?kx:rmem_mod(ix+kx,bins),by=full?ky:rmem_mod(iy+ky,bins);
                for(int id=rmem_head[bx+bins*by];id>=0;id=rmem_next[id]){
                    double xx=rmem_particle[id]->x-p->x,yy=rmem_particle[id]->y-p->y;
                    xx-=L*nearbyint(xx/L);yy-=L*nearbyint(yy/L);
                    if(xx*xx+yy*yy>=reach*reach)continue;
                    if(count>=4096){row_failure=2;goto row_ready;}
                    candidates[count++]=(RmemCandidate){id,xx,yy};
                }
            }
            qsort(candidates,count,sizeof(*candidates),rmem_cmp);
            for(int j=0;j<count;j++){dx[j]=candidates[j].dx;dy[j]=candidates[j].dy;if(candidates[j].id==owner)slot=j;}
            if(count<=12){reach*=1.5;continue;}
            if(corotate)code=rm_quadratic_row_gradient(count,dx,dy,slot,reach,12,row,work,meta,gradient_x,gradient_y);
            else code=rm_quadratic_row(count,dx,dy,slot,reach,12,row,work,meta);
            if(code==2){reach*=1.5;continue;}
            if(code){row_failure=4;goto row_ready;}
        }
row_ready:
        #pragma omp ordered
        {
        if(row_failure)row_errors|=row_failure;
        else {
        double gam=GAS_GAMMA(s),scale=sqrt((double)p->mass)*sqrt(gam*(gam-1)*(double)p->ie/(double)p->mass)/sqrt((double)p->volume);
        double row_norm=0;
        for(int j=0;j<count;j++){
            int id=candidates[j].id;treevorostressrk4particletype *q=rmem_particle[id];double b=scale*row[j];
            rmem_rate[2*owner]+=b*q->vx;rmem_rate[2*owner+1]+=b*q->vy;
            rmem_acc[2*id]-=b*rmem_state[owner].z[0];rmem_acc[2*id+1]-=b*rmem_state[owner].z[1];
            if(!gfs_kepler_frozen(q->x,q->y)){double entry=fabs(b)/sqrt((double)q->mass);row_norm+=entry;rmem_columns[id]+=entry;}
        }
        if(corotate){
            double spin=0;
            for(int j=0;j<count;j++){
                treevorostressrk4particletype *q=rmem_particle[candidates[j].id];
                spin+=.5*(gradient_x[j]*q->vy-gradient_y[j]*q->vx);
            }
            if(!isfinite(spin))row_errors|=16;
            else {
                /* Same-stage spatial-vector transport. This skew term does
                 * no auxiliary work, does not alter the force, and does not
                 * use a circular reference velocity or prescribed spin. */
                rmem_rate[2*owner]-=spin*rmem_state[owner].z[1];
                rmem_rate[2*owner+1]+=spin*rmem_state[owner].z[0];
                if(fabs(spin)>max_spin)max_spin=fabs(spin);
            }
        }
        if(row_norm>max_row)max_row=row_norm;
        }
        }
    }
    /* MPI_Abort stays outside the parallel region (MPI thread contract). */
    if(row_errors)rmem_fail(s,"parallel/serial row construction failed (coverage, capacity or rank)");
    double max_column=0;
    for(long i=0;i<rmem_n;i++){
        treevorostressrk4particletype *p=rmem_particle[i];
        if(gfs_kepler_frozen(p->x,p->y)){
            rmem_power-=rmem_acc[2*i]*p->vx+rmem_acc[2*i+1]*p->vy;
            rmem_acc[2*i]=rmem_acc[2*i+1]=0;
        } else {rmem_acc[2*i]/=p->mass;rmem_acc[2*i+1]/=p->mass;}
        if(rmem_columns[i]>max_column)max_column=rmem_columns[i];
    }
    double frequency=sqrt(max_row*max_column);
    if(corotate)frequency+=max_spin;
    double limit=frequency>0?.25/frequency:1e300;
    if(!isfinite(limit) || !(limit>0))rmem_fail(s,"invalid auxiliary timestep");
    if(corotate)fprintf(stderr,"[RMEM_SPIN_STAGE] max_spin=%.17g frequency_bound=%.17g dt_limit=%.17g scope=frozen_B_spin_block_not_full_nonlinear_bound\n",max_spin,frequency,limit);
    return limit;
}
static double rmem_a(long id,int component){return rmem_on()?rmem_acc[2*id+component]:0;}
static void rmem_record(SimParameters *s,int stage,double dt){
    if(!rmem_on())return;
    if(stage==1){rmem_trial_work=0;for(long i=0;i<rmem_n;i++)if(rms_begin(rmem_state+i,dt))rmem_fail(s,"trial begin failed");}
    for(long i=0;i<rmem_n;i++)if(rms_stage(rmem_state+i,stage,rmem_rate[2*i],rmem_rate[2*i+1]))rmem_fail(s,"stage update failed");
    rmem_trial_work+=dt*(stage==1 || stage==4?1.:2.)/6*rmem_power;
}
static void rmem_reject(SimParameters *s){
    if(!rmem_on())return;
    for(long i=0;i<rmem_n;i++)if(rms_reject(rmem_state+i))rmem_fail(s,"trial rollback failed");
    rmem_trial_work=0;
}
static void rmem_accept(SimParameters *s){
    if(!rmem_on())return;
    for(long i=0;i<rmem_n;i++)if(rms_accept(rmem_state+i))rmem_fail(s,"trial commit failed");
}
static void rmem_start(SimParameters *s){
    if(rmem_corotate() && !rmem_on())rmem_fail(s,"GFS_RMEM_COROTATE requires GFS_REVERSIBLE_MEMORY");
    if(!rmem_on())return;
    if(NID(s)!=1 || GAS_GPU_ENABLED(s) || GAS_AVMODE(s)!=5 || GAS_Kappa(s)!=0 || GAS_FCENTROID(s)!=0
       || GAS_XSPHEPS(s)!=0 || !gfs_step_retry_on() || !gfs_hydro_on() || !gfs_eaud_on()
       || gfs_dynw2_on() || gfs_subzone_on())rmem_fail(s,"unsupported experimental configuration");
    static int transport_banner=0;
    if(!transport_banner){fprintf(stderr,"[RMEM_TRANSPORT] corotate=%d scope=experimental_spatial_vector_transport_preserve_environment_on_restart\n",rmem_corotate());transport_banner=1;}
    rmem_map(s);rmem_aux_before=rmem_energy();
}

/* Existing CUDA geometry, shared C fit and central BE viscosity.
 * Initial coupled gate: fully owned one-rank cells and static wall faces.
 * CPU OpenMP solve is profiled separately; CUDA still builds the geometry
 * and pressure RHS. No particle ABI or alternate PDE solver is introduced.
 */
#include <float.h>
#include "../gfs_variational_core.h"
#include "../gfs_velocity_fit.h"
#include "../gfs_viscosity_implicit.h"
#include <stdint.h>
#include "gfs_unified_acoustic_audit.inc"
#include "gfs_unified_transport_audit.inc"
/* Diagnostic-only original-radius lookup. Missing files never alter physics. */
static double gfs_unified_original_radius(uint64_t id)
{
    static uint64_t *keys;static double *radius;static size_t mask;static int loaded=0;
    const char *path=getenv("GFS_UNIFIED_COHORT_FILE");if(!path)return NAN;
    if(!loaded){
        FILE *in=fopen(path,"rb");uint64_t count;
        if(!in||fread(&count,sizeof(count),1,in)!=1||count>10000000)exit(99);
        size_t capacity=1;while(capacity<2*count)capacity*=2;mask=capacity-1;
        keys=calloc(capacity,sizeof(*keys));radius=calloc(capacity,sizeof(*radius));
        if(!keys||!radius)exit(99);
        for(uint64_t k=0;k<count;k++){
            struct {uint64_t id;double r0;} row;
            if(fread(&row,sizeof(row),1,in)!=1||!isfinite(row.r0)||row.r0<0||row.id==UINT64_MAX)exit(99);
            size_t slot=(row.id*UINT64_C(11400714819323198485))&mask;
            while(keys[slot]){if(keys[slot]==row.id+1)exit(99);slot=(slot+1)&mask;}
            keys[slot]=row.id+1;radius[slot]=row.r0;
        }
        if(fgetc(in)!=EOF)exit(99);fclose(in);loaded=1;
    }
    size_t slot=(id*UINT64_C(11400714819323198485))&mask;
    while(keys[slot]){if(keys[slot]==id+1)return radius[slot];slot=(slot+1)&mask;}
    fprintf(stderr,"[UNIFIED_COHORT] unknown id=%llu\n",(unsigned long long)id);exit(99);
}
static int gfs_unified_interface_matches(int raw,int owner,double sx,double sy,
        const int *alias,const ParticleSoA *geometry,double tol)
{
    if(raw<0||alias[raw]!=owner)return 0;
    return fabs(geometry->x[raw]-geometry->x[owner]-sx)<=tol&&
           fabs(geometry->y[raw]-geometry->y[owner]-sy)<=tol;
}
static int gfs_unified_viscosity_apply(void *opaque,double gamma,double h,
        GfsViscosityStats *stats,int frozen_rate)
{
    SimParameters *simpar=opaque;GPUContext *ctx=gpu_get_context();
    int n=VORO_NP(simpar),nf=ctx->d_faces.n_faces_total,nt=n+VORO_NPAD(simpar);
    double lx=SIMBOX(simpar).x.max-SIMBOX(simpar).x.min,ly=SIMBOX(simpar).y.max-SIMBOX(simpar).y.min;
    if(!ctx->initialized||n<1||NID(simpar)!=1||!(gamma>1)||
       (nt!=n&&!gfs_unified_periodic_mode()))return 1;
    ParticleSoA *geometry=&ctx->h_parts;
    const int *alias=gpu_unified_owner_aliases(ctx,n,nt,lx,ly);
    double image_tolerance=128*DBL_EPSILON*(lx+ly);
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    FaceCSR *f=&ctx->h_faces;
    gpu_download_face_csr(ctx,n,nf,f->face_offset,f->c1x,f->c1y,f->c2x,f->c2y,
                         f->neighbor_idx,f->kp_idx,f->km_idx,f->is_ghost);
    static int capn=0,capf=0;static double *mass,*U,*x,*v,*G;
    static GfsViscosityEdge *edges;static unsigned char *valid;
    static GfsViscosityWorkspace *workspace;
    if(!workspace)workspace=gfs_viscosity_workspace_create();
    if(n>capn){
#define U_RESIZE(ptr,count,type) do{void *q=realloc(ptr,(size_t)(count)*sizeof(type));if(!q)return 1;ptr=q;}while(0)
        U_RESIZE(mass,n,double);U_RESIZE(U,n,double);U_RESIZE(x,2*n,double);
        U_RESIZE(v,2*n,double);U_RESIZE(G,4*n,double);capn=n;
    }
    if(nf>capf){U_RESIZE(edges,nf,GfsViscosityEdge);U_RESIZE(valid,nf,unsigned char);capf=nf;}
#undef U_RESIZE
    int bad=0;
#pragma omp parallel for schedule(static) reduction(max:bad)
    for(int i=0;i<n;i++){
        mass[i]=p[i].mass;U[i]=p[i].weight_thermal.U;
        x[2*i]=p[i].x;x[2*i+1]=p[i].y;v[2*i]=p[i].vx;v[2*i+1]=p[i].vy;
        if(p[i].weight.W!=0||p[i].weight.pi!=0||!(p[i].weight_thermal.V>0))bad=1;
    }
    if(bad){fprintf(stderr,"[UNIFIED_VISCOSITY] stage=%d input/geometry contract failed mask=%d\n",__LINE__,bad);return 1;}
#pragma omp parallel for schedule(static) reduction(max:bad)
    for(int i=0;i<n;i++){
        GfsVelocityFit fit;gfs_velocity_fit_init(&fit,2);double gradient[3][3];
        for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
            int raw=f->neighbor_idx[k],j=raw>=0&&raw<nt?alias[raw]:raw;double r[3]={0},dv[3]={0};
            if(raw>=0&&raw<nt&&j<n&&!f->is_ghost[k]){
                r[0]=geometry->x[raw]-x[2*i];r[1]=geometry->y[raw]-x[2*i+1];
                dv[0]=v[2*j]-v[2*i];dv[1]=v[2*j+1]-v[2*i+1];
            }else if(j<=-1000){
                double sx=f->c2y[k]-f->c1y[k],sy=f->c1x[k]-f->c2x[k],A=hypot(sx,sy);
                if(!(A>0)){bad=1;continue;}
                double nx=sx/A,ny=sy/A;
                double distance=.5*((f->c1x[k]+f->c2x[k])*nx+(f->c1y[k]+f->c2y[k])*ny);
                double vn=v[2*i]*nx+v[2*i+1]*ny;
                r[0]=2*distance*nx;r[1]=2*distance*ny;dv[0]=-2*vn*nx;dv[1]=-2*vn*ny;
            }else{bad=1;continue;}
            double area=hypot(f->c2y[k]-f->c1y[k],f->c1x[k]-f->c2x[k]);
            if(!gfs_velocity_fit_add_face(&fit,r,dv,area))bad=1;
        }
        if(!gfs_velocity_fit_solve(&fit,gradient,NULL)){bad=1;continue;}
        G[4*i]=gradient[0][0];G[4*i+1]=gradient[0][1];G[4*i+2]=gradient[1][0];G[4*i+3]=gradient[1][1];
    }
    if(bad){fprintf(stderr,"[UNIFIED_VISCOSITY] stage=%d input/geometry contract failed mask=%d\n",__LINE__,bad);return 1;}
    const char *planar_n=getenv("GFS_UNIFIED_PLANAR_AUDIT_N");
    if(planar_n){
        int side=atoi(planar_n);static int fit_call;static FILE *fit_trace;
        if(side<2||(long long)side*side!=n)return 1;
        double *low=malloc((size_t)side*4*sizeof(double)),*high=malloc((size_t)side*4*sizeof(double));
        if(!low||!high)return 1;
        for(int k=0;k<4*side;k++){low[k]=HUGE_VAL;high[k]=-HUGE_VAL;}
        int trace_calls=0;const char *trace_env=getenv("GFS_UNIFIED_PLANAR_TRACE_STEPS");
        if(trace_env)trace_calls=2*atoi(trace_env);
        ++fit_call;
        if(fit_call<=trace_calls&&!fit_trace){fit_trace=fopen("unified_planar_fit.csv","wx");if(!fit_trace)return 1;fprintf(fit_trace,"call,id,Gxx,Gxy,Gyx,Gyy\n");}
        double maxxy=0,maxyx=0,weighted_maxxy=0,weighted_maxxx_span=0;unsigned long long idxy=0;
        double *weighted_low=malloc((size_t)side*sizeof(double)),*weighted_high=malloc((size_t)side*sizeof(double));
        if(!weighted_low||!weighted_high)return 1;
        for(int c=0;c<side;c++){weighted_low[c]=HUGE_VAL;weighted_high[c]=-HUGE_VAL;}
        for(int i=0;i<n;i++){
            unsigned long long id=(unsigned long long)PINDX(p+i);int col=id%side;
            for(int a=0;a<4;a++){int k=4*col+a;low[k]=fmin(low[k],G[4*i+a]);high[k]=fmax(high[k],G[4*i+a]);}
            if(fabs(G[4*i+1])>maxxy){maxxy=fabs(G[4*i+1]);idxy=id;}
            maxyx=fmax(maxyx,fabs(G[4*i+2]));
            if(fit_call<=trace_calls)fprintf(fit_trace,"%d,%llu,%.17g,%.17g,%.17g,%.17g\n",fit_call,id,G[4*i],G[4*i+1],G[4*i+2],G[4*i+3]);
            /* Frozen-state causal comparison only. Production G is unchanged. */
            GfsVelocityFit weighted;gfs_velocity_fit_init(&weighted,2);
            for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
                int raw=f->neighbor_idx[k];if(raw<0||raw>=nt)return 1;int j=alias[raw];
                double r[3]={geometry->x[raw]-x[2*i],geometry->y[raw]-x[2*i+1],0};
                double dv[3]={v[2*j]-v[2*i],v[2*j+1]-v[2*i+1],0};
                double A=hypot(f->c2y[k]-f->c1y[k],f->c1x[k]-f->c2x[k]);
                if(!gfs_velocity_fit_add_weight(&weighted,r,dv,A/hypot(r[0],r[1])))return 1;
            }
            double Gw[3][3];if(!gfs_velocity_fit_solve(&weighted,Gw,NULL))return 1;
            weighted_maxxy=fmax(weighted_maxxy,fabs(Gw[0][1]));
            weighted_low[col]=fmin(weighted_low[col],Gw[0][0]);weighted_high[col]=fmax(weighted_high[col],Gw[0][0]);
        }
        double span[4]={0};for(int col=0;col<side;col++)for(int a=0;a<4;a++)span[a]=fmax(span[a],high[4*col+a]-low[4*col+a]);
        if(fit_trace&&fflush(fit_trace))return 1;
        for(int c=0;c<side;c++)weighted_maxxx_span=fmax(weighted_maxxx_span,weighted_high[c]-weighted_low[c]);
        fprintf(stderr,"[UNIFIED_PLANAR_FIT] call=%d column_span=(%.9e,%.9e,%.9e,%.9e) max_Gxy=%.9e id_Gxy=%llu max_Gyx=%.9e frozen_solid_angle_max_Gxy=%.9e frozen_solid_angle_Gxx_column_span=%.9e\n",fit_call,span[0],span[1],span[2],span[3],maxxy,idxy,maxyx,weighted_maxxy,weighted_maxxx_span);
        free(weighted_low);free(weighted_high);
        free(low);free(high);
    }
    memset(valid,0,nf);
    int unpaired=0;double max_unpaired_relative_area=0,max_area_mismatch=0,max_roundoff_ratio=0;
#pragma omp parallel for schedule(static) reduction(max:bad,max_unpaired_relative_area,max_area_mismatch,max_roundoff_ratio) reduction(+:unpaired)
    for(int i=0;i<n;i++)for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
        int raw=f->neighbor_idx[k],j=raw>=0&&raw<nt?alias[raw]:raw;
        if(raw>=nt){bad=1;continue;}
        double shift[3]={j>=0?geometry->x[raw]-x[2*j]:0,j>=0?geometry->y[raw]-x[2*j+1]:0,0};
        double A=hypot(f->c2y[k]-f->c1y[k],f->c1x[k]-f->c2x[k]);
        double c[3]={.5*(f->c1x[k]+f->c2x[k]),.5*(f->c1y[k]+f->c2y[k]),0};
        double V=p[i].weight_thermal.V,rhoi=mass[i]/V,Pi=(gamma-1)*U[i]/V,ci=sqrt(gamma*Pi/rhoi);
        double vi[3]={v[2*i],v[2*i+1],0},q=0,closing=0;GfsVariationalFace face;
        if(j>=0&&j<n){
            /* Canonical interface uses the sum of collinear segments on
             * each side. Degenerate clipping can split off a zero-length
             * segment with the same neighbor. Average the two side areas
             * before applying one exactly opposite pair impulse. */
            int duplicate=0,reverse=0;double Ai=0,Aj=0;
            for(int l=f->face_offset[i];l<f->face_offset[i+1];l++)if(gfs_unified_interface_matches(f->neighbor_idx[l],j,shift[0],shift[1],alias,geometry,image_tolerance)){
                if(l<k)duplicate=1;
                Ai+=hypot(f->c2y[l]-f->c1y[l],f->c1x[l]-f->c2x[l]);
            }
            if(duplicate)continue;
            for(int l=f->face_offset[j];l<f->face_offset[j+1];l++)if(gfs_unified_interface_matches(f->neighbor_idx[l],i,-shift[0],-shift[1],alias,geometry,image_tolerance)){
                reverse++;Aj+=hypot(f->c2y[l]-f->c1y[l],f->c1x[l]-f->c2x[l]);
            }
            if(reverse&&i>j)continue;
            double scale=fmax(sqrt(V),sqrt(p[j].weight_thermal.V));
            /* Coordinate differences inherit absolute-coordinate roundoff.
             * Scaling a clipping error only by a tiny cell length is
             * incorrect after subtraction of coordinates near x~3. */
            double coordinates=fabs(x[2*i])+fabs(x[2*i+1])+fabs(x[2*j])+fabs(x[2*j+1]);
            double bound=DBL_EPSILON*(1024*scale+64*coordinates);
            max_area_mismatch=fmax(max_area_mismatch,fabs(Ai-Aj));
            max_roundoff_ratio=fmax(max_roundoff_ratio,fabs(Ai-Aj)/bound);
            if(fabs(Ai-Aj)>bound){bad=16;continue;}
            if(!reverse){unpaired++;max_unpaired_relative_area=fmax(max_unpaired_relative_area,Ai/scale);}
            A=.5*(Ai+Aj);
            double r[3]={x[2*j]+shift[0]-x[2*i],x[2*j+1]+shift[1]-x[2*i+1],0},vj[3]={v[2*j],v[2*j+1],0};
            double jump[3]={vj[0]-vi[0]-.5*((G[4*i]+G[4*j])*r[0]+(G[4*i+1]+G[4*j+1])*r[1]),
                            vj[1]-vi[1]-.5*((G[4*i+2]+G[4*j+2])*r[0]+(G[4*i+3]+G[4*j+3])*r[1]),0};
            if(!gfs_variational_face(2,r,c,A,&face)){bad=1;continue;}
            double Vj=p[j].weight_thermal.V,rhoj=mass[j]/Vj,Pj=(gamma-1)*U[j]/Vj,cj=sqrt(gamma*Pj/rhoj);
            double sf=gfs_variational_compression_fraction(G[4*i]+G[4*i+3],G[4*i+2]-G[4*i+1],G[4*j]+G[4*j+3],G[4*j+2]-G[4*j+1]);
            q=gfs_variational_pair_pressure(&face,vi,vj,jump,rhoi,ci,rhoj,cj,sf,sqrt(V),sqrt(Vj));
            closing=(vi[0]-vj[0])*face.normal[0]+(vi[1]-vj[1])*face.normal[1];
        }else if(j<=-1000){
            double nx=(f->c2y[k]-f->c1y[k])/A,ny=(f->c1x[k]-f->c2x[k])/A;
            double distance=c[0]*nx+c[1]*ny,vn=vi[0]*nx+vi[1]*ny;
            double derivative=nx*(G[4*i]*nx+G[4*i+1]*ny)+ny*(G[4*i+2]*nx+G[4*i+3]*ny);
            double r[3]={2*distance*nx,2*distance*ny,0},vg[3]={vi[0]-2*vn*nx,vi[1]-2*vn*ny,0};
            double jump[3]={-2*(vn+distance*derivative)*nx,-2*(vn+distance*derivative)*ny,0};
            if(!gfs_variational_face(2,r,c,A,&face)){bad=1;continue;}
            double sf=gfs_variational_compression_fraction(G[4*i]+G[4*i+3],G[4*i+2]-G[4*i+1],G[4*i]+G[4*i+3],G[4*i+2]-G[4*i+1]);
            q=gfs_variational_pair_pressure(&face,vi,vg,jump,rhoi,ci,rhoi,ci,sf,sqrt(V),sqrt(V));
            closing=vn;j=-1;
        }else{bad=1;continue;}
        if(q>0&&closing>0){
            edges[k]=(GfsViscosityEdge){i,j,{face.normal[0],face.normal[1],0},q*A/closing,{shift[0],shift[1],0}};valid[k]=1;
        }
    }
    if(bad){fprintf(stderr,"[UNIFIED_VISCOSITY] stage=%d input/geometry contract failed mask=%d unpaired=%d max_relative_area=%.17g\n",__LINE__,bad,unpaired,max_unpaired_relative_area);return 1;}
    fprintf(stderr,"[UNIFIED_VISCOSITY_GRAPH] unpaired_roundoff_faces=%d max_relative_area=%.9e max_area_mismatch=%.9e max_roundoff_bound_ratio=%.9e\n",unpaired,max_unpaired_relative_area,max_area_mismatch,max_roundoff_ratio);
    int ne=0;for(int k=0;k<nf;k++)if(valid[k])edges[ne++]=edges[k];
    if(frozen_rate){
        if(gfs_unified_transport_audit(simpar,ctx,alias,gamma,G))return 1;
        const char *acoustic=getenv("GFS_UNIFIED_ACOUSTIC_AUDIT");
        if(acoustic&&!strcmp(acoustic,"1")&&gfs_unified_acoustic_audit(simpar,ctx,alias,gamma))return 1;
        /* Diagnostic rate on the actual pre-kick state. No BE solve or
         * state update. The same graph and constitutive law are used. */
        memset(U,0,(size_t)n*sizeof(double));long double total=0;
        for(int k=0;k<ne;k++){
            const GfsViscosityEdge *e=edges+k;double a=0;
            for(int d=0;d<2;d++)a+=e->normal[d]*(v[2*e->i+d]-(e->j>=0?v[2*e->j+d]:0));
            double rate=e->coefficient*a*a;total+=rate;
            U[e->i]+=e->j>=0?.5*rate:rate;if(e->j>=0)U[e->j]+=.5*rate;
        }
        FILE *out=fopen("unified_initial_heat.csv","wx");if(!out)return 1;
        fprintf(out,"id,x,y,V,U,ax,ay,Gxx,Gxy,Gyx,Gyy,heat_rate\n");
        for(int i=0;i<n;i++)fprintf(out,"%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            (unsigned long long)PINDX(p+i),p[i].x,p[i].y,p[i].weight_thermal.V,p[i].weight_thermal.U,p[i].ax,p[i].ay,G[4*i],G[4*i+1],G[4*i+2],G[4*i+3],U[i]);
        if(ferror(out)||fclose(out))return 1;
        fprintf(stderr,"[UNIFIED_INITIAL_HEAT_RATE] edges=%d rate=%.17Lg no_state_update=1\n",ne,total);return 0;
    }
    int rc=gfs_viscosity_be(workspace,2,n,ne,edges,mass,x,v,U,h,stats);
    if(rc){
        fprintf(stderr,"[UNIFIED_VISCOSITY] common solve failed rc=%d edges=%d h=%.17g\n",rc,ne,h);
        FILE *out=fopen("unified_viscosity_failure.bin","wbx");
        uint64_t header[4]={UINT64_C(0x564953434641494c),(uint64_t)n,(uint64_t)ne,2};
        if(!out||fwrite(header,sizeof(header),1,out)!=1||fwrite(&h,sizeof(h),1,out)!=1||
           fwrite(mass,sizeof(double),n,out)!=(size_t)n||fwrite(x,sizeof(double),2*n,out)!=(size_t)(2*n)||
           fwrite(v,sizeof(double),2*n,out)!=(size_t)(2*n)||fwrite(U,sizeof(double),n,out)!=(size_t)n||
           fwrite(edges,sizeof(*edges),ne,out)!=(size_t)ne||fclose(out)){
            fprintf(stderr,"[UNIFIED_VISCOSITY] failure-state write error\n");MPI_Abort(MPI_COMM(simpar),99);
        }
        return rc;
    }
/* Initialize the diagnostic hash serially before parallel readers. */
    int cohorts=getenv("GFS_UNIFIED_COHORT_FILE")!=NULL;
    if(cohorts)(void)gfs_unified_original_radius((uint64_t)PINDX(p));
    static int rbins=0,rthreads=0;static double rmax=0,*radial_partials,*radial_sum;
    const char *radial_env=getenv("GFS_UNIFIED_RADIAL_BINS");
    if(radial_env&&!cohorts)return 1;
    if(radial_env&&!rbins){
        rbins=atoi(radial_env);const char *limit=getenv("GFS_UNIFIED_RADIAL_MAX");rmax=limit?atof(limit):0;
        if(rbins<1||rbins>64||!(rmax>0)||!isfinite(rmax))return 1;
        rthreads=1;
#ifdef _OPENMP
        rthreads=omp_get_max_threads();
#endif
        radial_partials=calloc((size_t)rthreads*rbins,sizeof(double));radial_sum=calloc(rbins,sizeof(double));
        if(!radial_partials||!radial_sum)return 1;
    }
    if(rbins)memset(radial_partials,0,(size_t)rthreads*rbins*sizeof(double));
    double Qcore=0,Qdisk=0,Qambient=0;
#pragma omp parallel for schedule(static) reduction(+:Qcore,Qdisk,Qambient)
    for(int i=0;i<n;i++){
        if(cohorts){
            double r0=gfs_unified_original_radius((uint64_t)PINDX(p+i));
            double Q=U[i]-p[i].weight_thermal.U;
            if(rbins){
                int tid=0,bin=(int)floor(r0/rmax*rbins);if(bin>=rbins)bin=rbins-1;
#ifdef _OPENMP
                tid=omp_get_thread_num();
#endif
                radial_partials[tid*rbins+bin]+=Q;
            }else if(r0<.5)Qcore+=Q;else if(r0<=2)Qdisk+=Q;else Qambient+=Q;
        }
        p[i].vx=v[2*i];p[i].vy=v[2*i+1];p[i].weight_thermal.U=U[i];p[i].ie=(float)U[i];
        double V=p[i].weight_thermal.V,rho=mass[i]/V,P=(gamma-1)*U[i]/V;
        p[i].den=(float)rho;p[i].pressure=(float)P;p[i].csound=(float)sqrt(gamma*P/rho);
    }
    if(rbins){
        for(int b=0;b<rbins;b++)for(int t=0;t<rthreads;t++)radial_sum[b]+=radial_partials[t*rbins+b];
        fprintf(stderr,"[UNIFIED_RADIAL_HEAT] interval=%.17g bins=%d r0_max=%.17g segment_cumulative=(",h,rbins,rmax);
        for(int b=0;b<rbins;b++)fprintf(stderr,"%s%.17g",b?",":"",radial_sum[b]);fprintf(stderr,")\n");
    }else if(cohorts)fprintf(stderr,"[UNIFIED_COHORT_HEAT] interval=%.17g core=%.17g original_disk=%.17g original_ambient=%.17g\n",h,Qcore,Qdisk,Qambient);
    return 0;
}
int gpu_unified_viscosity_half_step(void *opaque,double gamma,double h,GfsViscosityStats *stats){return gfs_unified_viscosity_apply(opaque,gamma,h,stats,0);}
int gpu_unified_viscosity_rate_audit(void *opaque,double gamma){return gfs_unified_viscosity_apply(opaque,gamma,0,NULL,1);}

/* Exact owned-state and immutable-W checks before reusing local geometry.
 * Image coordinates stay cached after the native wrapper frees CPU padding.
 * The caller separately checks time, boundary, and geometry parameters. */
int gpu_unified_geometry_reusable(void *opaque)
{
    SimParameters *simpar=opaque;GPUContext *ctx=gpu_get_context();
    int n=VORO_NP(simpar),nt=n+VORO_NPAD(simpar);
    if(!ctx->initialized||NID(simpar)!=1||n<1||ctx->d_faces.n_particles!=n||
       ctx->h_parts.n_total!=nt||ctx->weight_geometry_count!=nt||
       (nt!=n&&!gfs_unified_periodic_mode()))return 0;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);ParticleSoA *a=&ctx->h_parts;
    int valid=1;
#pragma omp parallel for schedule(static) reduction(&:valid)
    for(int i=0;i<n;i++){
        if(a->indx[i]!=(long long)PINDX(p+i)||a->x[i]!=p[i].x||a->y[i]!=p[i].y||
           a->mass[i]!=p[i].mass||a->weight_volume64[i]!=p[i].weight_thermal.V||
           a->weight_W[i]!=p[i].weight.W||ctx->weight_geometry_snapshot[i]!=p[i].weight.W)valid=0;
    }
    return valid;
}

/* Reuse the current face CSR after a velocity/thermal update. No new mesh. */
double gpu_unified_pressure_refresh(void *opaque)
{
    SimParameters *simpar=opaque;GPUContext *ctx=gpu_get_context();
    int n=VORO_NP(simpar),nt=n+VORO_NPAD(simpar);
    if(!ctx->initialized||NID(simpar)!=1||(nt!=n&&!gfs_unified_periodic_mode()))return 0;
    if(nt==n)fillParticleSoA(simpar,&ctx->h_parts,1);
    else{
        /* The normal force wrapper frees CPU padding. Keep its immutable
         * image coordinates; refresh owned states and alias them on device. */
        treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
        ParticleSoA *a=&ctx->h_parts;
        for(int i=0;i<n;i++){
            if(a->indx[i]!=(long long)PINDX(p+i)||a->x[i]!=p[i].x||a->y[i]!=p[i].y){
                fprintf(stderr,"[UNIFIED_REFRESH] geometry no longer matches owned state\n");MPI_Abort(MPI_COMM(simpar),99);
            }
            a->vx[i]=p[i].vx;a->vy[i]=p[i].vy;a->mass[i]=p[i].mass;
            a->weight_U64[i]=p[i].weight_thermal.U;a->weight_W[i]=p[i].weight.W;
            a->ie[i]=p[i].ie;
        }
    }
    gpu_upload_particles(ctx,nt,1);
    GPUPhysicsParams params={0};params.Gamma=GAS_GAMMA(simpar);
    params.Courant=GAS_COURANT(simpar);params.phase1=1;
    params.domain_x=SIMBOX(simpar).x.max-SIMBOX(simpar).x.min;
    params.domain_y=SIMBOX(simpar).y.max-SIMBOX(simpar).y.min;
    double dt=gpu_launch_unified_force_kernel(ctx,n,nt,&params);
    gpu_download_results(ctx,n);
    writeBackResults_range(simpar,&ctx->h_parts,1,0,n);
    return dt;
}
#include "gfs_unified_acoustic_cfl.inc"

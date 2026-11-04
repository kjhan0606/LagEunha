/* Same W0 pressure, face-angular LSQ, positive central BE viscosity,
 * and P-D-drift-D-P update as the 1D and native CUDA 2D candidates.
 * Only existing linked-cell Voro3D/OpenMP geometry and LAG3DV1 IO are reused.
 * No legacy Riemann, RK4, floor, jitter, or finite face-area cutoff enters.
 */
#include <float.h>
#include <limits.h>
#include "../gfs_variational_core.h"
#include "../gfs_velocity_fit.h"
#include "../gfs_viscosity_implicit.h"
typedef struct {
    double *x,*v;long *offset;GfsViscosityEdge *edges;unsigned char *valid;
    long capn,capf;GfsViscosityWorkspace *solve;
} G3UnifiedWorkspace;
static double g3_u_area(const g3face *f){return sqrt(f->S[0]*f->S[0]+f->S[1]*f->S[1]+f->S[2]*f->S[2]);}
static double g3_u_clock(void){
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock()/CLOCKS_PER_SEC;
#endif
}
static int g3_u_primitives(g3sim *s){
    int bad=0;
#pragma omp parallel for schedule(static) reduction(max:bad)
    for(long i=0;i<s->n;i++){
        double m=s->m[i],V=s->vol[i],U=s->ie[i];
        if(!(m>0)||!(V>0)||!(U>0)||s->w2[i]!=0||!isfinite(m)||!isfinite(V)||!isfinite(U)||
           !isfinite(s->vx[i])||!isfinite(s->vy[i])||!isfinite(s->vz[i])){bad=1;continue;}
        s->den[i]=m/V;s->P[i]=(s->gamma-1)*U/V;s->cs[i]=sqrt(s->gamma*s->P[i]/s->den[i]);
        if(!isfinite(s->cs[i]))bad=1;
    }
    return bad;
}
static int g3_u_fit(g3sim *s){
    int bad=0;
#pragma omp parallel for schedule(static) reduction(max:bad)
    for(long i=0;i<s->n;i++){
        GfsVelocityFit fit;gfs_velocity_fit_init(&fit,3);g3face *f=g3_faces(s,i);
        for(int k=0;k<s->fcnt[i];k++){
            int j=f[k].j;if(j<0||j>=s->n||f[k].mf){bad=1;continue;}
            double dv[3]={s->vx[j]-s->vx[i],s->vy[j]-s->vy[i],s->vz[j]-s->vz[i]};
            if(!gfs_velocity_fit_add_face(&fit,f[k].rel,dv,g3_u_area(f+k)))bad=1;
        }
        double G[3][3];if(!gfs_velocity_fit_solve(&fit,G,NULL)){bad=1;continue;}
        for(int a=0;a<3;a++)for(int b=0;b<3;b++)s->grad[15*i+6+3*a+b]=G[a][b];
    }
    return bad;
}
static int g3_u_pressure(g3sim *s,double *dt){
    if(g3_u_primitives(s))return 1;
    int bad=0;double minimum=HUGE_VAL;
#pragma omp parallel for schedule(static) reduction(max:bad) reduction(min:minimum)
    for(long i=0;i<s->n;i++){
        double vi[3]={s->vx[i],s->vy[i],s->vz[i]},F[3]={0},sweep=0,h=cbrt(s->vol[i]);
        double local=s->courant*h/s->cs[i];g3face *faces=g3_faces(s,i);
        for(int k=0;k<s->fcnt[i];k++){
            g3face *f=faces+k;int j=f->j;GfsVariationalFace b;
            if(j<0||j>=s->n||f->mf||!gfs_variational_face(3,f->rel,f->c,g3_u_area(f),&b)){bad=1;continue;}
            double vj[3]={s->vx[j],s->vy[j],s->vz[j]},fp[3];
            gfs_variational_pressure(&b,s->P[i],s->P[j],fp);
            double relative=0;for(int a=0;a<3;a++){F[a]+=fp[a];relative+=(vj[a]-vi[a])*b.normal[a];}
            sweep+=fabs(gfs_variational_sweep(&b,vi,vj));
            local=fmin(local,s->courant*h/(s->cs[i]+fabs(relative)));
            if(relative<0)local=fmin(local,.2*b.separation/(-relative));
        }
        if(sweep>0)local=fmin(local,s->courant*s->vol[i]/sweep);
        s->ax[i]=F[0]/s->m[i];s->ay[i]=F[1]/s->m[i];s->az[i]=F[2]/s->m[i];
        if(!(local>0)||!isfinite(local)||!isfinite(F[0])||!isfinite(F[1])||!isfinite(F[2]))bad=1;
        minimum=fmin(minimum,local);
    }
    *dt=minimum;return bad;
}
static double g3_u_curl(const double *G){
    double a=G[7]-G[5],b=G[2]-G[6],c=G[3]-G[1];return sqrt(a*a+b*b+c*c);
}
static int g3_u_match(const g3sim *s,long i,const g3face *f,int owner,const double shift[3]){
    if(f->j!=owner||f->mf)return 0;
    double delta[3]={s->x[owner]-s->x[i],s->y[owner]-s->y[i],s->z[owner]-s->z[i]};
    for(int a=0;a<3;a++){
        double tol=128*DBL_EPSILON*(s->L[a]+fabs(delta[a])+fabs(f->rel[a]));
        if(fabs(f->rel[a]-delta[a]-shift[a])>tol)return 0;
    }
    return 1;
}
static int g3_u_viscosity(g3sim *s,G3UnifiedWorkspace *w,double h,GfsViscosityStats *stats){
    if(g3_u_primitives(s)||g3_u_fit(s))return 1;
    long n=s->n,nf=0;
    if(n>w->capn){
        free(w->x);free(w->v);free(w->offset);
        w->x=g3_alloc(3*n*sizeof(double));w->v=g3_alloc(3*n*sizeof(double));w->offset=g3_alloc((n+1)*sizeof(long));w->capn=n;
    }
    for(long i=0;i<n;i++){w->offset[i]=nf;nf+=s->fcnt[i];}w->offset[n]=nf;
    if(nf>INT_MAX)return 1;
    if(nf>w->capf){free(w->edges);free(w->valid);w->edges=g3_alloc(nf*sizeof(*w->edges));w->valid=g3_alloc(nf);w->capf=nf;}
    memset(w->valid,0,nf);int bad=0;long unpaired=0;double mismatch=0;
#pragma omp parallel for schedule(static) reduction(max:bad,mismatch) reduction(+:unpaired)
    for(long i=0;i<n;i++){
        w->x[3*i]=s->x[i];w->x[3*i+1]=s->y[i];w->x[3*i+2]=s->z[i];
        w->v[3*i]=s->vx[i];w->v[3*i+1]=s->vy[i];w->v[3*i+2]=s->vz[i];
        g3face *faces=g3_faces(s,i);
        for(int k=0;k<s->fcnt[i];k++){
            g3face *f=faces+k;int j=f->j;if(j<0||j>=n||f->mf){bad=1;continue;}
            if(i==j){bad=1;continue;}if(i>j)continue;
            double shift[3]={f->rel[0]-s->x[j]+s->x[i],f->rel[1]-s->y[j]+s->y[i],f->rel[2]-s->z[j]+s->z[i]};
            for(int a=0;a<3;a++){
                double exact=s->L[a]*nearbyint(shift[a]/s->L[a]);
                if(fabs(shift[a]-exact)>128*DBL_EPSILON*(s->L[a]+fabs(f->rel[a])))bad=1;
                shift[a]=exact;
            }
            int duplicate=0,reverse=0;double Ai=0,Aj=0,opposite[3]={-shift[0],-shift[1],-shift[2]};
            for(int l=0;l<s->fcnt[i];l++)if(g3_u_match(s,i,faces+l,j,shift)){if(l<k)duplicate=1;Ai+=g3_u_area(faces+l);}
            if(duplicate)continue;
            g3face *other=g3_faces(s,j);
            for(int l=0;l<s->fcnt[j];l++)if(g3_u_match(s,j,other+l,i,opposite)){reverse++;Aj+=g3_u_area(other+l);}
            double scale=fmax(cbrt(s->vol[i]),cbrt(s->vol[j]));
            double bound=2048*DBL_EPSILON*scale*scale;
            mismatch=fmax(mismatch,fabs(Ai-Aj));if(fabs(Ai-Aj)>bound){bad=2;continue;}
            if(!reverse)unpaired++;
            double area=.5*(Ai+Aj);if(!(area>0))continue;
            double *Gi=s->grad+15*i+6,*Gj=s->grad+15*j+6;
            double vi[3]={s->vx[i],s->vy[i],s->vz[i]},vj[3]={s->vx[j],s->vy[j],s->vz[j]},jump[3];
            for(int a=0;a<3;a++){jump[a]=vj[a]-vi[a];for(int b=0;b<3;b++)jump[a]-=.5*(Gi[3*a+b]+Gj[3*a+b])*f->rel[b];}
            GfsVariationalFace b;if(!gfs_variational_face(3,f->rel,f->c,area,&b)){bad=1;continue;}
            double sf=gfs_variational_compression_fraction(Gi[0]+Gi[4]+Gi[8],g3_u_curl(Gi),Gj[0]+Gj[4]+Gj[8],g3_u_curl(Gj));
            double q=gfs_variational_pair_pressure(&b,vi,vj,jump,s->den[i],s->cs[i],s->den[j],s->cs[j],sf,cbrt(s->vol[i]),cbrt(s->vol[j]));
            double closing=0;for(int a=0;a<3;a++)closing+=(vi[a]-vj[a])*b.normal[a];
            if(q>0&&closing>0){long at=w->offset[i]+k;w->edges[at]=(GfsViscosityEdge){(int)i,j,{b.normal[0],b.normal[1],b.normal[2]},q*area/closing,{shift[0],shift[1],shift[2]}};w->valid[at]=1;}
        }
    }
    if(bad){
        fprintf(stderr,"[UNIFIED_3D_GRAPH] invalid mask=%d mismatch=%.17g\n",bad,mismatch);
        FILE *out=fopen("unified3d_geometry_failure.bin","wbx");
        uint64_t header[4]={UINT64_C(0x473347454f464149),(uint64_t)n,(uint64_t)nf,1};
        if(!out||fwrite(header,sizeof(header),1,out)!=1||fwrite(&s->gamma,sizeof(double),1,out)!=1)return 1;
        for(long i=0;i<n;i++){
            double row[9]={s->x[i],s->y[i],s->z[i],s->vx[i],s->vy[i],s->vz[i],s->m[i],s->ie[i],s->vol[i]};
            if(fwrite(row,sizeof(row),1,out)!=1||fwrite(s->id+i,sizeof(int64_t),1,out)!=1||fwrite(s->fcnt+i,sizeof(int),1,out)!=1||fwrite(g3_faces(s,i),sizeof(g3face),s->fcnt[i],out)!=(size_t)s->fcnt[i])return 1;
        }
        if(fclose(out))return 1;
        return 1;
    }
    int ne=0;for(long k=0;k<nf;k++)if(w->valid[k])w->edges[ne++]=w->edges[k];
    if(!w->solve)w->solve=gfs_viscosity_workspace_create();
    int rc=gfs_viscosity_be(w->solve,3,(int)n,ne,w->edges,s->m,w->x,w->v,s->ie,h,stats);
    if(rc){fprintf(stderr,"[UNIFIED_3D_GRAPH] solve failed rc=%d edges=%d\n",rc,ne);return rc;}
#pragma omp parallel for schedule(static)
    for(long i=0;i<n;i++){s->vx[i]=w->v[3*i];s->vy[i]=w->v[3*i+1];s->vz[i]=w->v[3*i+2];}
    fprintf(stderr,"[UNIFIED_3D_GRAPH] edges=%d unpaired_roundoff=%ld max_area_mismatch=%.9e periodic_torque=(%.9e,%.9e,%.9e)\n",ne,unpaired,mismatch,stats->periodic_torque[0],stats->periodic_torque[1],stats->periodic_torque[2]);
    return 0;
}
static int g3_u_geometry_audit(g3sim *s){
    double closure=0,volume_error=0,parallel_error=0;long nf=0;long double V=0,work=0,scale_work=0;
    for(long i=0;i<s->n;i++){
        double sum[3]={0},area=0,Vs=0,p=1+.1*sin((double)s->id[i]);
        double vi[3]={sin(.31*s->id[i]),cos(.47*s->id[i]),sin(.19*s->id[i])};
        g3face *faces=g3_faces(s,i);
        for(int k=0;k<s->fcnt[i];k++){
            g3face *f=faces+k;double A=g3_u_area(f);GfsVariationalFace b;
            if(!gfs_variational_face(3,f->rel,f->c,A,&b)||f->mf)return 1;
            double vj[3]={sin(.31*s->id[f->j]),cos(.47*s->id[f->j]),sin(.19*s->id[f->j])};
            double pj=1+.1*sin((double)s->id[f->j]),F[3];gfs_variational_pressure(&b,p,pj,F);
            double rate=gfs_variational_sweep(&b,vi,vj);work-=p*rate;scale_work+=fabs(p*rate);
            double pe=0;for(int a=0;a<3;a++){sum[a]+=f->S[a];Vs+=f->c[a]*f->S[a]/3;pe+=pow(f->S[a]-A*b.normal[a],2);work+=F[a]*vi[a];scale_work+=fabs(F[a]*vi[a]);}
            parallel_error=fmax(parallel_error,sqrt(pe)/pow(s->vol[i],2./3));area+=A;nf++;
        }
        closure=fmax(closure,sqrt(sum[0]*sum[0]+sum[1]*sum[1]+sum[2]*sum[2])/area);
        volume_error=fmax(volume_error,fabs(Vs/s->vol[i]-1));V+=s->vol[i];
    }
    double sumV_error=fabs((double)V/(s->L[0]*s->L[1]*s->L[2])-1);
    double adjoint=fabs((double)work)/fmax(1,(double)scale_work);
    int pass=closure<4096*DBL_EPSILON&&volume_error<4096*DBL_EPSILON&&parallel_error<4096*DBL_EPSILON&&sumV_error<4096*DBL_EPSILON&&adjoint<4096*DBL_EPSILON&&s->n_jitter==0;
    fprintf(stderr,"[UNIFIED_3D_GEOMETRY] N=%ld faces=%ld closure=%.17g volume_identity=%.17g parallel_scaled=%.17g sumV_error=%.17g adjoint_work=%.17g jitter=%ld endpoint_clamps=%lld max_endpoint_residual_over_bound=%.17g pass=%d\n",s->n,nf,closure,volume_error,parallel_error,sumV_error,adjoint,s->n_jitter,s->strict_endpoint_clamps,s->strict_endpoint_max_ratio,pass);
    return !pass;
}
static int g3_u_volume_fd(g3sim *s){
    /* Diagnostic perturbations move actual generators and rebuild the same
     * strict mesh. Restore the exact original positions before any evolution. */
    long n=s->n;double *original=g3_alloc(3*n*sizeof(double)),*direction=g3_alloc(3*n*sizeof(double));
    double *prediction=g3_alloc(n*sizeof(double)),*plus=g3_alloc(n*sizeof(double));
    double h=cbrt(s->L[0]*s->L[1]*s->L[2]/n),etas[]={1e-3,1e-4,1e-5,1e-6,1e-7,1e-8};
    for(long i=0;i<n;i++){
        original[3*i]=s->x[i];original[3*i+1]=s->y[i];original[3*i+2]=s->z[i];
        direction[3*i]=sin(.31*s->id[i]);direction[3*i+1]=cos(.47*s->id[i]);direction[3*i+2]=sin(.19*s->id[i]);
    }
    for(long i=0;i<n;i++){
        g3face *faces=g3_faces(s,i);
        for(int k=0;k<s->fcnt[i];k++){
            g3face *f=faces+k;GfsVariationalFace b;
            if(!gfs_variational_face(3,f->rel,f->c,g3_u_area(f),&b))return 1;
            prediction[i]+=gfs_variational_sweep(&b,direction+3*i,direction+3*f->j);
        }
    }
    int failed=0;double best=HUGE_VAL;
    for(int q=0;q<6&&!failed;q++){
        double eps=etas[q]*h;
        for(int sign=1;sign>=-1;sign-=2){
            for(long i=0;i<n;i++){s->x[i]=original[3*i]+sign*eps*direction[3*i];s->y[i]=original[3*i+1]+sign*eps*direction[3*i+1];s->z[i]=original[3*i+2]+sign*eps*direction[3*i+2];}
            if(g3_tessellate(s)){failed=1;break;}
            if(sign==1)for(long i=0;i<n;i++)plus[i]=s->vol[i];
        }
        if(failed)break;
        double error=0;for(long i=0;i<n;i++)error=fmax(error,fabs((plus[i]-s->vol[i])/(2*eps)-prediction[i])/(h*h));
        best=fmin(best,error);fprintf(stderr,"[UNIFIED_3D_VOLUME_FD] eta=%.9e max_derivative_error_over_h2=%.17g\n",etas[q],error);
    }
    for(long i=0;i<n;i++){s->x[i]=original[3*i];s->y[i]=original[3*i+1];s->z[i]=original[3*i+2];}
    if(g3_tessellate(s))failed=1;
    free(original);free(direction);free(prediction);free(plus);
    fprintf(stderr,"[UNIFIED_3D_VOLUME_FD] best=%.17g pass=%d\n",best,!failed&&best<1e-6);
    return failed||!(best<1e-6);
}
static double g3_u_step(g3sim *s,G3UnifiedWorkspace *w,double cap){
    double begin=g3_u_clock(),dt;
    if(g3_u_pressure(s,&dt))return -1;
    double factor=g3_env_d("GFS_UNIFIED_DT_FACTOR",1);if(!(factor>0)||factor>1)return -1;
    dt=fmin(cap,dt*factor);if(!(dt>0)||!isfinite(dt))return -1;
    long n=s->n;
    /* Preserve the entire accepted dynamical state for a failed step. */
    for(long i=0;i<n;i++){
        double *q=s->k+8*i;q[0]=s->x[i];q[1]=s->y[i];q[2]=s->z[i];q[3]=s->vx[i];q[4]=s->vy[i];q[5]=s->vz[i];q[6]=s->ie[i];q[7]=s->vol[i];
        s->vx[i]+=.5*dt*s->ax[i];s->vy[i]+=.5*dt*s->ay[i];s->vz[i]+=.5*dt*s->az[i];
    }
    GfsViscosityStats D0={0},D1={0};int bad=0;
    if(g3_u_viscosity(s,w,.5*dt,&D0))goto fail;
#pragma omp parallel for schedule(static)
    for(long i=0;i<n;i++){
        s->ie0[i]=s->ie[i];s->V0[i]=s->vol[i];
        s->x[i]+=dt*s->vx[i];s->y[i]+=dt*s->vy[i];s->z[i]+=dt*s->vz[i];
    }
    g3_wrap(s);if(s->n!=n||g3_tessellate(s))goto fail;
    if(g3_env_on("GFS_UNIFIED_GEOMETRY_EACH_STAGE")&&g3_u_geometry_audit(s))goto fail;
#pragma omp parallel for schedule(static) reduction(max:bad)
    for(long i=0;i<n;i++)if(!gfs_variational_adiabatic(s->ie0[i],s->V0[i],s->vol[i],s->gamma,s->ie+i))bad=1;
    if(bad||g3_u_viscosity(s,w,.5*dt,&D1)||g3_u_pressure(s,&cap))goto fail;
#pragma omp parallel for schedule(static)
    for(long i=0;i<n;i++){s->vx[i]+=.5*dt*s->ax[i];s->vy[i]+=.5*dt*s->ay[i];s->vz[i]+=.5*dt*s->az[i];}
    if(g3_u_primitives(s))goto fail;
    s->time+=dt;s->step++;
    long double E=0,px=0,py=0,pz=0;double minU=HUGE_VAL,minV=HUGE_VAL,maxvy=0,maxvz=0;
    for(long i=0;i<n;i++){
        E+=s->ie[i]+.5L*s->m[i]*(s->vx[i]*s->vx[i]+s->vy[i]*s->vy[i]+s->vz[i]*s->vz[i]);
        px+=s->m[i]*s->vx[i];py+=s->m[i]*s->vy[i];pz+=s->m[i]*s->vz[i];
        minU=fmin(minU,s->ie[i]);minV=fmin(minV,s->vol[i]);maxvy=fmax(maxvy,fabs(s->vy[i]));maxvz=fmax(maxvz,fabs(s->vz[i]));
    }
    fprintf(stderr,"[UNIFIED_3D_KDK] step=%ld t=%.17g dt=%.17g E=%.17Lg dErel=%.9e momentum=(%.9Lg,%.9Lg,%.9Lg) minU=%.9e minV=%.9e max_vy=%.9e max_vz=%.9e heat=%.9e viscosity_energy_defect=%.9e CG=%d/%d active=%d/%d wall_s=%.6f\n",s->step,s->time,dt,E,((double)E-s->E0tot)/fabs(s->E0tot),px,py,pz,minU,minV,maxvy,maxvz,D0.heat+D1.heat,D0.energy_defect+D1.energy_defect,D0.iterations,D1.iterations,D0.active_iterations,D1.active_iterations,g3_u_clock()-begin);
    return dt;
fail:
    for(long i=0;i<n;i++){
        double *q=s->k+8*i;s->x[i]=q[0];s->y[i]=q[1];s->z[i]=q[2];s->vx[i]=q[3];s->vy[i]=q[4];s->vz[i]=q[5];s->ie[i]=q[6];s->vol[i]=q[7];
    }
    s->geom_valid=0;return -1;
}
static int g3_unified_run(g3sim *s,const lag3d_header *h0){
    const char *boundary=getenv("GFS_UNIFIED_BOUNDARY");
    if(!boundary||strcmp(boundary,"periodic")||s->n<1||s->n>INT_MAX||s->gamma<=1||s->courant<=0||s->kappa!=0||s->alphavis!=0||s->betavis!=0||s->gravity||s->pm_on||s->acc_on||s->de_eta||s->es_on||s->half_on||s->stage_reset){fprintf(stderr,"[UNIFIED_3D] unsupported extra physics or boundary\n");return 2;}
    for(int a=0;a<3;a++)if(s->bc[a]!=BC_PERIODIC)return 2;
    fprintf(stderr,"[UNIFIED_3D] common W0 fixed-mass pressure/face-angular LSQ/BE/KDK; strict double geometry with bounded endpoints, no jitter/floor/finite area cutoff; CPU/OpenMP only\n");
    g3_wrap(s);if(g3_tessellate(s)||g3_u_primitives(s)||g3_u_geometry_audit(s))return 4;
    if(g3_env_on("GFS_UNIFIED_GEOMETRY_FD")&&(g3_u_volume_fd(s)||g3_u_primitives(s)))return 4;
    long double E=0;for(long i=0;i<s->n;i++)E+=s->ie[i]+.5L*s->m[i]*(s->vx[i]*s->vx[i]+s->vy[i]*s->vy[i]+s->vz[i]*s->vz[i]);s->E0tot=E;
    int isnap=0;if(g3_write_snap(s,isnap++,h0))return 6;
    if(g3_env_on("GFS_UNIFIED_GEOMETRY_ONLY"))return 0;
    G3UnifiedWorkspace w={0};double next=s->dumpdt>0?fmin(s->tend,s->time+s->dumpdt):s->tend;
    while(s->time<s->tend&&s->step<s->maxsteps){
        if(g3_u_step(s,&w,next-s->time)<0){fprintf(stderr,"[UNIFIED_3D] failed, restoring last accepted time %.17g\n",s->time);g3_write_snap(s,isnap++,h0);return 5;}
        if(s->time>=next){if(g3_write_snap(s,isnap++,h0))return 6;next=s->dumpdt>0?fmin(s->tend,next+s->dumpdt):s->tend;}
    }
    if(s->time!=s->tend){fprintf(stderr,"[UNIFIED_3D] did not reach exact target\n");return 5;}
    fprintf(stderr,"[UNIFIED_3D] complete step=%ld t=%.17g jitter=%ld\n",s->step,s->time,s->n_jitter);
    gfs_viscosity_workspace_free(w.solve);free(w.x);free(w.v);free(w.offset);free(w.edges);free(w.valid);return 0;
}

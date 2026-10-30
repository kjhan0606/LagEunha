/* Default-off one-rank bounded-W0 native RK4 event transaction.
 * Rejected attempts never reach permanent RK4 ledgers or driver time.
 */
#ifdef USE_CUDA
#include <limits.h>
#include <stdint.h>
#include "kepler_slip_wall.h"
#include "gfs_wall_dense.h"
typedef struct { double x,y,m; unsigned int present,seen; } GfsWallBase;
typedef struct {
    int on,tries,maxtries,search,id,side,hits;
    ptrdiff_t np,total;size_t bytes,extent;unsigned int generation;void *backup;GfsWallBase *base;
    float dtold;double floors;int nfloors;long long pairs;
    double audit[5],external,ceiling,lo,hi,glo,ghi;
    double pressure[3],gravity[2],collision[3],roundoff;
    double start[7];const char *reason;
} GfsWallEventStep;
static double gfs_wall_nx[1024],gfs_wall_ny[1024];
static const double gfs_wall_band=1.e-6;
static double gfs_wall_apothem(void){return 2.75*cos(M_PI/1024.);}
static double gfs_wall_tolerance(double x,double y){
    return 32.*DBL_EPSILON*fmax(2.75,fabs(x-3.)+fabs(y-3.));
}
static void gfs_wall_totals(SimParameters *s,double z[7]){
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
    long double a[7]={0,0,0,0,0,0,0};
    for(ptrdiff_t i=0;i<VORO_NP(s);i++){
        double m=p[i].mass,x=p[i].x-3.,y=p[i].y-3.,vx=p[i].vx,vy=p[i].vy;
        a[0]+=m*vx;a[1]+=m*vy;a[2]+=m*(x*vy-y*vx);
        a[3]+=p[i].ie+.5*m*(vx*vx+vy*vy)+m*kepler_phi(p[i].x,p[i].y);
        a[4]+=m;a[5]+=p[i].ie;a[6]+=.5*m*(vx*vx+vy*vy);
    }
    for(int k=0;k<7;k++)z[k]=(double)a[k];
}
static void gfs_wall_restore(SimParameters *s,GfsWallEventStep *w){
    gpu_geometry_trial_end();
    gfs_ssp_restore(s,w->backup,w->np,w->total,w->dtold,w->floors,w->nfloors,w->pairs);
    gfs_eaud_hyd=w->audit[0];gfs_eaud_hyd_abs=w->audit[1];gfs_eaud_wg=w->audit[2];
    gfs_eaud_eh0=w->audit[3];gfs_eaud_ep0=w->audit[4];
    double restored_audit[5]={gfs_eaud_hyd,gfs_eaud_hyd_abs,gfs_eaud_wg,gfs_eaud_eh0,gfs_eaud_ep0};
    int cache_invalidated=!gpu_get_context()->initialized;
    int exact=memcmp(VORORK4_TBP(s),w->backup,w->bytes)==0 &&
        VORO_NP(s)==w->np && VORO_TNP(s)==w->total && GAS_dtold(s)==w->dtold &&
        NPSUM(s)==w->total && VORO_NPAD(s)==0 && cache_invalidated &&
        memcmp(restored_audit,w->audit,sizeof(restored_audit))==0 &&
        gfs_stage_floor_l==w->floors && gfs_stage_floor_n_l==w->nfloors &&
        gfs_pair_face_ends==w->pairs;
    fprintf(stderr,"[WALL_ROLLBACK] step=%d attempt=%d reason=%s payload_counts_dtold_ledgers_exact=%d device_cache_invalidated=%d\n",
        gfs_rk4_context_step,w->tries+1,w->reason?w->reason:"trial",exact,cache_invalidated);
    if(!exact)MPI_Abort(MPI_COMM(s),99);
}
static void gfs_wall_stop(SimParameters *s,GfsWallEventStep *w,const char *reason){
    w->reason=reason;gfs_wall_restore(s,w);
    gfs_rk4_failure_snapshot(s,w->backup,(int)w->np,"wall_event_start",w->dtold,w->dtold);
    fprintf(stderr,"[WALL_EVENT_STOP] step=%d t=%.17g reason=%s accepted start preserved\n",
        gfs_rk4_context_step,gfs_rk4_context_time,reason);fflush(stderr);
    gfs_wall_runtime_active=0;MPI_Abort(MPI_COMM(s),99);
}
static void gfs_wall_begin(SimParameters *s,GfsWallEventStep *w,
                           void (*post)(SimParameters *)){
    memset(w,0,sizeof(*w));w->id=w->side=-1;
    const char *setting=getenv("GFS_KEPLER_WALL_EVENTS");
    w->on=setting && strcmp(setting,"1")==0;if(!w->on)return;
    gfs_es_config();
    if(NID(s)!=1 || post!=fixed_polygon_postStage_blend || !GAS_GPU_ENABLED(s) ||
       !gpu_varpress_w0_enabled() || GAS_AVMODE(s)!=5 || GAS_Kappa(s)!=0 ||
       GAS_ENTROPY_MODE(s)!=0 || !sedov_phase1_on() || GAS_XSPHEPS(s)!=0 ||
       GAS_FCENTROID(s)!=0 || GAS_VISCOSITY(s)!=0 || GAS_AlphaVis(s)!=0 ||
       GAS_BetaVis(s)!=0 || GAS_CDAMAX(s)!=0 || GAS_HYPERVALPHA(s)!=0 ||
       GAS_ACCX(s)!=0 || GAS_ACCY(s)!=0 || gfs_dual_energy_eta()>0 || gfs_es_on ||
       gfs_half_on || gfs_dynw2_on() || gfs_w2lin_on() || sedov_lagvol_on() ||
       gfs_vaud_on() || gfs_thermal64_requested() ||
       (getenv("GFS_RK4_TRANSACTION") && atoi(getenv("GFS_RK4_TRANSACTION"))) ||
       (getenv("GFS_RK4_DIRECT_U") && atoi(getenv("GFS_RK4_DIRECT_U"))) ||
       !getenv("GFS_GPU_RAW_FACES") || atoi(getenv("GFS_GPU_RAW_FACES"))!=1 ||
       !getenv("GFS_GPU_FINAL_VOLUME") || atoi(getenv("GFS_GPU_FINAL_VOLUME"))!=1 ||
       TVORORK4_DDINFO(s)[0].n_size!=sizeof(treevorostressrk4particletype)){
        fprintf(stderr,"[WALL_EVENTS] requires one-rank native bounded smooth W0, fixed mass, Re=inf and legacy thermal RK4\n");
        MPI_Abort(MPI_COMM(s),99);
    }
    static int ready=0;
    if(!ready){for(int k=0;k<1024;k++){double a=(k+.5)*2.*M_PI/1024.;gfs_wall_nx[k]=cos(a);gfs_wall_ny[k]=sin(a);}ready=1;}
    w->np=VORO_NP(s);w->total=VORO_TNP(s);
    if(w->np<=0 || w->np>INT_MAX || (size_t)w->np>SIZE_MAX/sizeof(treevorostressrk4particletype)){
        fprintf(stderr,"[WALL_EVENTS] invalid owned count for bounded event transaction\n");MPI_Abort(MPI_COMM(s),99);
    }
    w->bytes=(size_t)w->np*sizeof(treevorostressrk4particletype);
    w->backup=malloc(w->bytes);
    if(!w->backup)MPI_Abort(MPI_COMM(s),99);
    memcpy(w->backup,VORORK4_TBP(s),w->bytes);
    w->dtold=GAS_dtold(s);w->floors=gfs_stage_floor_l;w->nfloors=gfs_stage_floor_n_l;w->pairs=gfs_pair_face_ends;
    w->audit[0]=gfs_eaud_hyd;w->audit[1]=gfs_eaud_hyd_abs;w->audit[2]=gfs_eaud_wg;
    w->audit[3]=gfs_eaud_eh0;w->audit[4]=gfs_eaud_ep0;
    w->external=gfs_rk4_step_ceiling;w->ceiling=w->external>0?w->external:1.e100;
    w->maxtries=64;const char *limit=getenv("GFS_WALL_MAX_RETRIES");
    if(limit){char *end=NULL;long n=strtol(limit,&end,10);if(*end || n<1 || n>128)MPI_Abort(MPI_COMM(s),99);w->maxtries=(int)n;}
    treevorostressrk4particletype *p=w->backup;
    /* Preserve global IDs. Bound the direct lookup allocation before indexing;
     * the intended truncated Q1 and remapped Q2 fixtures use less than 3*N.
     */
    size_t id_limit=(size_t)w->np<=(SIZE_MAX-1024)/8?(size_t)w->np*8+1024:SIZE_MAX;
    for(ptrdiff_t i=0;i<w->np;i++){
        size_t id=(size_t)PINDX((treevorork4particletype*)(p+i));
        if(id>(size_t)INT_MAX || id>=id_limit)gfs_wall_stop(s,w,"unsupported_accepted_ID_range");
        if(id+1>w->extent)w->extent=id+1;
    }
    if(w->extent>SIZE_MAX/sizeof(*w->base))gfs_wall_stop(s,w,"ID_lookup_allocation_overflow");
    w->base=calloc(w->extent,sizeof(*w->base));
    if(!w->base)gfs_wall_stop(s,w,"ID_lookup_allocation_failed");
    for(ptrdiff_t i=0;i<w->np;i++){
        size_t id=(size_t)PINDX((treevorork4particletype*)(p+i));
        if(w->base[id].present || !isfinite(p[i].mass) || !(p[i].mass>0) || !(p[i].ie>0) ||
           !(p[i].volume>0) || p[i].w2!=0 || !isfinite(p[i].ie) || !isfinite(p[i].volume) ||
           !isfinite(p[i].x) || !isfinite(p[i].y) || !isfinite(p[i].vx) || !isfinite(p[i].vy))
            gfs_wall_stop(s,w,"invalid_accepted_ID_mass_W0_state");
        w->base[id]=(GfsWallBase){p[i].x,p[i].y,p[i].mass,1,0};
        double rx=p[i].x-3.,ry=p[i].y-3.,tol=gfs_wall_tolerance(p[i].x,p[i].y);
        if(hypot(rx,ry)>=gfs_wall_apothem())for(int k=0;k<1024;k++)
            if(rx*gfs_wall_nx[k]+ry*gfs_wall_ny[k]-gfs_wall_apothem()>tol)
                gfs_wall_stop(s,w,"accepted_start_outside_domain");
    }
    gfs_wall_totals(s,w->start);
    const char *dump=getenv("GFS_WALL_DUMP_ACCEPTED_START");
    if(dump && atoi(dump)>0)
        gfs_rk4_failure_snapshot(s,w->backup,(int)w->np,"wall_event_baseline",w->dtold,w->dtold);
}
static void gfs_wall_new_generation(GfsWallEventStep *w){
    if(++w->generation==0){
        for(size_t i=0;i<w->extent;i++)w->base[i].seen=0;
        w->generation=1;
    }
}
static int gfs_wall_trial_state_valid(treevorostressrk4particletype *p,int rhs){
    return p->ie>0 && p->volume>0 && isfinite(p->mass) && p->mass>0 && p->w2==0 &&
        isfinite(p->ie) && isfinite(p->volume) && isfinite(p->x) && isfinite(p->y) &&
        isfinite(p->vx) && isfinite(p->vy) &&
        (!rhs || (isfinite(p->ax) && isfinite(p->ay) && isfinite(p->die)));
}
#define GFS_WALL_SCAN_PARTICLE_TYPE treevorostressrk4particletype
#define GFS_WALL_SCAN_ID(p) PINDX((treevorork4particletype*)(p))
#define GFS_WALL_SCAN_MASS(p) ((p)->mass)
#define GFS_WALL_SCAN_STATE_VALID(p,rhs) gfs_wall_trial_state_valid((p),(rhs))
#include "gfs_wall_validation_scan.h"
#undef GFS_WALL_SCAN_PARTICLE_TYPE
#undef GFS_WALL_SCAN_ID
#undef GFS_WALL_SCAN_MASS
#undef GFS_WALL_SCAN_STATE_VALID
static int gfs_wall_identity(GfsWallEventStep *w,treevorostressrk4particletype *p){
    size_t id=(size_t)PINDX((treevorork4particletype*)p);
    if(id>=w->extent || !w->base[id].present || w->base[id].seen==w->generation || p->mass!=w->base[id].m){
        w->reason="trial_ID_or_mass_changed";return 0;
    }
    w->base[id].seen=w->generation;return 1;
}
static void gfs_wall_restart(GfsWallEventStep *w){
    if(!w->on)return;
    memset(w->pressure,0,sizeof(w->pressure));memset(w->gravity,0,sizeof(w->gravity));
    memset(w->collision,0,sizeof(w->collision));w->hits=0;w->roundoff=0;
    gfs_wall_stage_failed=0;gfs_wall_runtime_active=1;
    gfs_rk4_step_ceiling=w->ceiling;gpu_geometry_trial_begin();
}
static int gfs_wall_valid(SimParameters *s,GfsWallEventStep *w,int rhs){
    if(!w->on)return 1;
    if(gpu_geometry_trial_failed() || gfs_wall_stage_failed){w->reason="geometry_or_band";return 0;}
    if(VORO_NP(s)!=w->np){w->reason="owned_count";return 0;}
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
    gfs_wall_new_generation(w);
    const char *parallel=getenv("GFS_WALL_PARALLEL_VALIDATION");
    if(parallel && strcmp(parallel,"1")==0){
        if(gfs_wall_parallel_validation_scan(w->np,w->extent,w->generation,w->base,p,rhs))return 1;
        gfs_wall_new_generation(w);
    }
    for(ptrdiff_t i=0;i<w->np;i++){
        /* N distinct members of the immutable N-member initial ID set imply
         * exact set equality, even after native particle reordering.
         */
        if(!gfs_wall_identity(w,p+i))return 0;
        if(!gfs_wall_trial_state_valid(p+i,rhs)){
            w->reason="nonpositive_or_nonfinite_trial";return 0;
        }
    }
    return 1;
}
static int gfs_wall_endpoint_inside(SimParameters *s,GfsWallEventStep *w,double h){
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
    for(ptrdiff_t i=0;i<w->np;i++){
        double x=p[i].x-3.,y=p[i].y-3.,tol=gfs_wall_tolerance(p[i].x,p[i].y);
        if(hypot(x,y)<gfs_wall_apothem()-tol)continue;
        for(int k=0;k<1024;k++)if(x*gfs_wall_nx[k]+y*gfs_wall_ny[k]-gfs_wall_apothem()>tol){
            w->reason="unresolved_outside_endpoint";w->ceiling=.5*h;w->search=0;return 0;
        }
    }
    return 1;
}
static void gfs_wall_pressure_stage(SimParameters *s,GfsWallEventStep *w,double h,double weight){
    if(!w->on)return;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
    long double fx=0,fy=0,torque=0,gxsum=0,gysum=0;
    for(ptrdiff_t i=0;i<VORO_NP(s);i++){
        postype gx=0,gy=0;kepler_accel(p[i].x,p[i].y,&gx,&gy);
        double m=p[i].mass;
        fx+=m*p[i].ax;fy+=m*p[i].ay;
        torque+=m*((p[i].x-3.)*p[i].ay-(p[i].y-3.)*p[i].ax);
        gxsum+=m*gx;gysum+=m*gy;
    }
    w->pressure[0]-=(double)fx*h*weight;w->pressure[1]-=(double)fy*h*weight;
    w->pressure[2]-=(double)torque*h*weight;
    w->gravity[0]+=(double)gxsum*h*weight;w->gravity[1]+=(double)gysum*h*weight;
}
static double gfs_wall_proposal(SimParameters *s,GfsWallEventStep *w){
    double next=w->lo-w->glo*(w->hi-w->lo)/(w->ghi-w->glo);
    if(!(next>w->lo && next<w->hi))next=.5*(w->lo+w->hi);
    const double minimum=nextafter(gfs_rk4_context_time,INFINITY)-gfs_rk4_context_time;
    if(next<minimum)next=minimum;
    if(!(next>0) || !isfinite(next) || next>w->hi ||
       !(gfs_rk4_context_time+next>gfs_rk4_context_time))
        gfs_wall_stop(s,w,"event_time_resolution_limit");
    return next;
}
/* 1 commits this completed trial; 0 restores and retries at w->ceiling. */
static int gfs_wall_decide(SimParameters *s,GfsWallEventStep *w,double h){
    if(!w->on)return 1;
    static int forced_once=0;
    const char *force=getenv("GFS_WALL_FORCE_REJECT_ONCE"),*always=getenv("GFS_WALL_FORCE_REJECT_ALWAYS");
    if((always && atoi(always)>0) || (force && atoi(force)>0 && !forced_once)){
        forced_once=1;w->reason="verification_rejection";w->ceiling=h;return 0;
    }
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
    double earliest=2.;int first_id=-1,first_side=-1;
    if(VORO_NP(s)!=w->np)gfs_wall_stop(s,w,"owned_count");
    gfs_wall_new_generation(w);
    for(ptrdiff_t i=0;i<w->np;i++){
        size_t id=(size_t)PINDX((treevorork4particletype*)(p+i));
        if(!gfs_wall_identity(w,p+i))gfs_wall_stop(s,w,"trial_ID_or_mass_changed");
        double x=w->base[id].x-3.,y=w->base[id].y-3.;
        double bound=(5./24.)*hypot(p[i].rk4.k1x,p[i].rk4.k1y)+
            (1./3.)*(hypot(p[i].rk4.k2x,p[i].rk4.k2y)+hypot(p[i].rk4.k3x,p[i].rk4.k3y))+
            (1./6.)*hypot(p[i].rk4.k4x,p[i].rk4.k4y);
        if(hypot(x,y)+bound<gfs_wall_apothem())continue;
        for(int k=0;k<1024;k++){
            double nx=gfs_wall_nx[k],ny=gfs_wall_ny[k],theta=0.;
            GfsWallCubic poly=gfs_wall_cubic(nx*x+ny*y-gfs_wall_apothem(),
                nx*p[i].rk4.k1x+ny*p[i].rk4.k1y,nx*p[i].rk4.k2x+ny*p[i].rk4.k2y,
                nx*p[i].rk4.k3x+ny*p[i].rk4.k3y,nx*p[i].rk4.k4x+ny*p[i].rk4.k4y);
            int found=gfs_wall_first_outward(poly,&theta);
            if(found<0)gfs_wall_stop(s,w,"nonfinite_dense_path");
            if(found && theta<earliest){earliest=theta;first_id=(int)id;first_side=k;}
        }
    }
    if(w->search){
        ptrdiff_t index=-1;
        for(ptrdiff_t i=0;i<w->np;i++)if((int)PINDX((treevorork4particletype*)(p+i))==w->id){index=i;break;}
        if(index<0)gfs_wall_stop(s,w,"event_owner_missing");
        double nx=gfs_wall_nx[w->side],ny=gfs_wall_ny[w->side];
        double gap=nx*(p[index].x-3.)+ny*(p[index].y-3.)-gfs_wall_apothem();
        double tol=gfs_wall_tolerance(p[index].x,p[index].y);
        double vn=nx*p[index].vx+ny*p[index].vy;
        int competing=first_id>=0 && (first_id!=w->id || first_side!=w->side) && earliest*h<h-2.*tol/fmax(fabs(vn),1.e-8);
        if(!competing && fabs(gap)<=tol){
            if(!(vn>1.e-8))gfs_wall_stop(s,w,"unresolved_grazing_contact");
            /* Every contact is checked before any impulse, including corners. */
            for(ptrdiff_t i=0;i<w->np;i++){
                double x=p[i].x-3.,y=p[i].y-3.,t=gfs_wall_tolerance(p[i].x,p[i].y);
                int contacts=0;
                if(hypot(x,y)<gfs_wall_apothem()-t)continue;
                for(int k=0;k<1024;k++){
                    double g=gfs_wall_nx[k]*x+gfs_wall_ny[k]*y-gfs_wall_apothem();
                    if(g>t){w->reason="competing_outside_endpoint";w->ceiling=.5*h;w->search=0;return 0;}
                    if(fabs(g)<=t)contacts++;
                }
                if(contacts>1)gfs_wall_stop(s,w,"unsupported_corner_contact");
            }
            for(ptrdiff_t i=0;i<w->np;i++){
                double x=p[i].x-3.,y=p[i].y-3.,t=gfs_wall_tolerance(p[i].x,p[i].y);
                if(hypot(x,y)<gfs_wall_apothem()-t)continue;
                for(int k=0;k<1024;k++){
                    double g=gfs_wall_nx[k]*x+gfs_wall_ny[k]*y-gfs_wall_apothem();
                    if(fabs(g)>t || gfs_wall_nx[k]*p[i].vx+gfs_wall_ny[k]*p[i].vy<=1.e-8)continue;
                    GfsSlipImpulse impulse;
                    if(gfs_slip_plane_impact(p[i].x,p[i].y,3.,3.,gfs_wall_nx[k],gfs_wall_ny[k],
                        gfs_wall_apothem(),p[i].mass,p[i].vx,p[i].vy,&impulse)!=GFS_SLIP_REFLECTED)
                        gfs_wall_stop(s,w,"invalid_contact_impulse");
                    p[i].vx=impulse.vx_after;p[i].vy=impulse.vy_after;
                    w->collision[0]+=impulse.wall_jx;w->collision[1]+=impulse.wall_jy;
                    w->collision[2]+=impulse.wall_torque_impulse;w->roundoff+=impulse.fluid_kinetic_roundoff;w->hits++;
                    fprintf(stderr,"[WALL_IMPACT] step=%d t=%.17g id=%zu side=%d gap=%.17g tolerance=%.17g vn_before=%.17g kinetic_roundoff=%.17g\n",
                        gfs_rk4_context_step,gfs_rk4_context_time+h,(size_t)PINDX((treevorork4particletype*)(p+i)),k,g,t,
                        (impulse.wall_jx*gfs_wall_nx[k]+impulse.wall_jy*gfs_wall_ny[k])/(2.*p[i].mass),impulse.fluid_kinetic_roundoff);
                    break;
                }
            }
            if(!w->hits)gfs_wall_stop(s,w,"missing_resolved_contact");
            return 1;
        }
        if(!competing){
            if(gap < -tol && first_id>=0 && earliest*h<h-2.*tol/fmax(fabs(vn),1.e-8)){
                w->search=0;w->reason="out_and_back_subdivision";
                w->ceiling=h*fmax(.01,fmin(.5,.5*earliest));return 0;
            }
            if(gap>0 && h<w->hi){w->hi=h;w->ghi=gap;}
            if(gap<0 && h>w->lo){w->lo=h;w->glo=gap;}
            if(h<=w->lo && first_id<0 && gap < -tol){w->search=0;return gfs_wall_endpoint_inside(s,w,h);}
            w->ceiling=gfs_wall_proposal(s,w);w->reason="actual_trial_root_refinement";return 0;
        }
        w->search=0;
    }
    if(first_id<0)return gfs_wall_endpoint_inside(s,w,h);
    ptrdiff_t index=-1;
    for(ptrdiff_t i=0;i<w->np;i++)if((int)PINDX((treevorork4particletype*)(p+i))==first_id){index=i;break;}
    double nx=gfs_wall_nx[first_side],ny=gfs_wall_ny[first_side];
    double gap=nx*(p[index].x-3.)+ny*(p[index].y-3.)-gfs_wall_apothem();
    if(gap<=0){w->reason="out_and_back_subdivision";w->ceiling=h*fmax(.01,fmin(.5,.5*earliest));return 0;}
    w->search=1;w->id=first_id;w->side=first_side;w->lo=0.;w->hi=h;
    w->glo=nx*(w->base[first_id].x-3.)+ny*(w->base[first_id].y-3.)-gfs_wall_apothem();w->ghi=gap;
    if(w->glo>gfs_wall_tolerance(w->base[first_id].x,w->base[first_id].y))gfs_wall_stop(s,w,"invalid_initial_event_gap");
    w->ceiling=gfs_wall_proposal(s,w);w->reason="first_outward_root_bracket";return 0;
}
static void gfs_wall_reject(SimParameters *s,GfsWallEventStep *w,double h){
    if(!w->reason || strcmp(w->reason,"geometry_or_band")==0 || strcmp(w->reason,"nonpositive_or_nonfinite_trial")==0){
        if(!(h>0))gfs_wall_stop(s,w,"invalid_first_RHS_at_accepted_start");
        w->ceiling=.5*h;w->search=0;
    }
    gfs_wall_restore(s,w);w->tries++;
    if(w->tries>=w->maxtries || !(w->ceiling>0) ||
       !(gfs_rk4_context_time+w->ceiling>gfs_rk4_context_time))gfs_wall_stop(s,w,"retry_or_time_progress_limit");
}
static void gfs_wall_commit(SimParameters *s,GfsWallEventStep *w,double h){
    if(!w->on)return;
    double end[7];gfs_wall_totals(s,end);
    double px=end[0]-w->start[0]+w->pressure[0]+w->collision[0]-w->gravity[0];
    double py=end[1]-w->start[1]+w->pressure[1]+w->collision[1]-w->gravity[1];
    double angular=end[2]-w->start[2]+w->pressure[2]+w->collision[2];
    fprintf(stderr,"[WALL_ACCEPT] step=%d t_start=%.17g h=%.17g attempts=%d impacts=%d wall_pressure_J=(%.17g,%.17g,%.17g) wall_collision_J=(%.17g,%.17g,%.17g) momentum_residual=(%.9e,%.9e) angular_residual=%.9e energy_change=%.9e kinetic_roundoff=%.9e\n",
        gfs_rk4_context_step,gfs_rk4_context_time,h,w->tries+1,w->hits,
        w->pressure[0],w->pressure[1],w->pressure[2],w->collision[0],w->collision[1],w->collision[2],
        px,py,angular,end[3]-w->start[3],w->roundoff);
    gpu_geometry_trial_end();gfs_rk4_step_ceiling=w->external;gfs_wall_runtime_active=0;
    free(w->backup);free(w->base);
}
#endif

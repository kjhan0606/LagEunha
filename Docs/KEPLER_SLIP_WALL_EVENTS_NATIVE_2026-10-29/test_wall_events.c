/* Execute the actual native event include with a one-particle analytic ODE.
 * MPI/cache/geometry adapters are test doubles; real native GPU rollback is
 * checked separately on the frozen Q2 checkpoint. This is not a gas solver.
 */
#define _POSIX_C_SOURCE 200809L
#define USE_CUDA 1
#include <assert.h>
#include <float.h>
#include <math.h>
#include <setjmp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
typedef double postype;
typedef struct {
    double k1x,k1y,k1vx,k1vy,k2x,k2y,k2vx,k2vy;
    double k3x,k3y,k3vx,k3vy,k4x,k4y,k4vx,k4vy;
} TestRK;
typedef struct {
    size_t id;double x,y,vx,vy,mass,ie,volume,w2,ax,ay,die;TestRK rk4;
} treevorostressrk4particletype;
typedef treevorostressrk4particletype treevorork4particletype;
typedef struct {treevorostressrk4particletype *p;ptrdiff_t np,total;float dtold;struct {size_t n_size;} info[1];} SimParameters;
#define VORORK4_TBP(s) ((s)->p)
#define VORO_NP(s) ((s)->np)
#define VORO_TNP(s) ((s)->total)
#define NPSUM(s) ((s)->total)
#define VORO_NPAD(s) 0
#define GAS_dtold(s) ((s)->dtold)
#define TVORORK4_DDINFO(s) ((s)->info)
#define PINDX(p) ((p)->id)
#define NID(s) 1
#define MPI_COMM(s) 0
#define GAS_GPU_ENABLED(s) 1
#define GAS_AVMODE(s) 5
#define GAS_Kappa(s) 0
#define GAS_ENTROPY_MODE(s) 0
#define GAS_XSPHEPS(s) 0
#define GAS_FCENTROID(s) 0
#define GAS_VISCOSITY(s) 0
#define GAS_AlphaVis(s) 0
#define GAS_BetaVis(s) 0
#define GAS_CDAMAX(s) 0
#define GAS_HYPERVALPHA(s) 0
#define GAS_ACCX(s) 0
#define GAS_ACCY(s) 0
static int gfs_es_on=0,gfs_half_on=0,gfs_wall_runtime_active=0,gfs_wall_stage_failed=0;
static int gfs_rk4_context_step=0,gfs_stage_floor_n_l=0;
static double gfs_rk4_context_time=0,gfs_rk4_step_ceiling=0,gfs_stage_floor_l=0;
static long long gfs_pair_face_ends=0;
static double gfs_eaud_hyd=0,gfs_eaud_hyd_abs=0,gfs_eaud_wg=0,gfs_eaud_eh0=0,gfs_eaud_ep0=0;
static int gpu_failed=0,restore_count=0,expect_abort=0;
static int field_mode=0;static double acceleration=.2,harmonic=.3;
static jmp_buf failure_jump;
static treevorostressrk4particletype saved_failure;
static double saved_time=0;
static int no(void){return 0;}
static int yes(void){return 1;}
#define gpu_varpress_w0_enabled yes
#define sedov_phase1_on yes
#define gfs_dynw2_on no
#define gfs_w2lin_on no
#define sedov_lagvol_on no
#define gfs_vaud_on no
#define gfs_thermal64_requested no
static double gfs_dual_energy_eta(void){return 0.;}
static void gfs_es_config(void){}
static void gpu_geometry_trial_begin(void){gpu_failed=0;}
static void gpu_geometry_trial_end(void){}
static int gpu_geometry_trial_failed(void){return gpu_failed;}
static struct {int initialized;} mock_context={0};
#define gpu_get_context() (&mock_context)
static void MPI_Abort(int comm,int code);
static void kepler_accel(double x,double y,double *ax,double *ay){
    double theta=M_PI/1024.;
    *ax=field_mode==1?acceleration*cos(theta):(field_mode==2?-harmonic*(x-3.):0.);
    *ay=field_mode==1?acceleration*sin(theta):(field_mode==2?-harmonic*(y-3.):0.);
}
static double kepler_phi(double x,double y){
    return field_mode==1?-acceleration*((x-3.)*cos(M_PI/1024.)+(y-3.)*sin(M_PI/1024.)):
        (field_mode==2?.5*harmonic*((x-3.)*(x-3.)+(y-3.)*(y-3.)):0.);
}
static void fixed_polygon_postStage_blend(SimParameters *s){
    double x=s->p->x-3.,y=s->p->y-3.,a=2.75*cos(M_PI/1024.);
    for(int k=0;k<1024;k++){
        double t=(k+.5)*2.*M_PI/1024.;
        if(x*cos(t)+y*sin(t)-a>1.e-6)gfs_wall_stage_failed=1;
    }
}
static void gfs_ssp_restore(SimParameters *s,const void *b,ptrdiff_t np,ptrdiff_t total,
                           float dtold,double floors,int nfloors,long long pairs){
    free(s->p);s->p=malloc((size_t)np*sizeof(*s->p));memcpy(s->p,b,(size_t)np*sizeof(*s->p));
    s->np=np;s->total=total;s->dtold=dtold;
    gfs_stage_floor_l=floors;gfs_stage_floor_n_l=nfloors;gfs_pair_face_ends=pairs;
    mock_context.initialized=0;restore_count++;
}
static void gfs_rk4_failure_snapshot(SimParameters *s,const treevorostressrk4particletype *p,
                                    int n,const char *label,double dt,float dtold){
    (void)s;(void)n;(void)label;(void)dt;(void)dtold;saved_failure=p[0];saved_time=gfs_rk4_context_time;
}
#include "Exam/gfs_wall_events.inc"
static GfsWallEventStep *active=NULL;
static void MPI_Abort(int comm,int code){
    (void)comm;
    if(expect_abort){if(active){free(active->backup);free(active->base);}longjmp(failure_jump,1);}
    fprintf(stderr,"Unexpected test abort %d\n",code);exit(code);
}

static int stage(SimParameters *s,GfsWallEventStep *w){
    fixed_polygon_postStage_blend(s);return gfs_wall_valid(s,w,0);
}
static int trial(SimParameters *s,GfsWallEventStep *w,double maximum,double *actual){
    double h=fmin(maximum,gfs_rk4_step_ceiling);*actual=h;
    mock_context.initialized=1;
    treevorostressrk4particletype *p=s->p;double x=p->x,y=p->y,vx=p->vx,vy=p->vy;
    double kx[4],ky[4],kvx[4],kvy[4];
    const double weights[4]={1./6.,1./3.,1./3.,1./6.};
    for(int k=0;k<4;k++){
        if(k){double f=k==3?1.:.5;int j=k-1;
            p->x=x+f*kx[j];p->y=y+f*ky[j];p->vx=vx+f*kvx[j];p->vy=vy+f*kvy[j];}
        if(!stage(s,w))return 0;
        double ax,ay;kepler_accel(p->x,p->y,&ax,&ay);
        kx[k]=h*p->vx;ky[k]=h*p->vy;kvx[k]=h*ax;kvy[k]=h*ay;
        gfs_wall_pressure_stage(s,w,h,weights[k]);
    }
    p->x=x;p->y=y;p->vx=vx;p->vy=vy;
    for(int k=0;k<4;k++){p->x+=weights[k]*kx[k];p->y+=weights[k]*ky[k];p->vx+=weights[k]*kvx[k];p->vy+=weights[k]*kvy[k];}
    p->rk4=(TestRK){kx[0],ky[0],kvx[0],kvy[0],kx[1],ky[1],kvx[1],kvy[1],
                    kx[2],ky[2],kvx[2],kvy[2],kx[3],ky[3],kvx[3],kvy[3]};
    return stage(s,w);
}
static int hit_count=0;static double hit_times[32];
static double advance(SimParameters *s,double maximum,double remaining){
    gfs_rk4_step_ceiling=remaining;GfsWallEventStep w;gfs_wall_begin(s,&w,fixed_polygon_postStage_blend);active=&w;
    double h=0.;
    for(;;){
        gfs_wall_restart(&w);
        if(trial(s,&w,maximum,&h) && gfs_wall_decide(s,&w,h))break;
        gfs_wall_reject(s,&w,h);
    }
    if(w.hits){assert(hit_count<32);hit_times[hit_count++]=gfs_rk4_context_time+h;}
    gfs_wall_commit(s,&w,h);active=NULL;
    assert(h>0 && gfs_rk4_context_time+h>gfs_rk4_context_time);
    gfs_rk4_context_time+=h;gfs_rk4_context_step++;s->dtold=(float)h;
    return h;
}
static SimParameters initialize(double gap,double speed,int mode){
    SimParameters s={0};s.p=calloc(1,sizeof(*s.p));s.np=s.total=1;s.info[0].n_size=sizeof(*s.p);
    double t=M_PI/1024.,r=2.75*cos(t)-gap;
    s.p->x=3.+r*cos(t);s.p->y=3.+r*sin(t);s.p->vx=speed*cos(t);s.p->vy=speed*sin(t);
    s.p->mass=1;s.p->ie=.1;s.p->volume=1;field_mode=mode;
    gfs_rk4_context_time=0;gfs_rk4_context_step=0;hit_count=0;restore_count=0;
    return s;
}
static void run(SimParameters *s,double h,double end){
    for(int k=0;gfs_rk4_context_time<end;k++){
        assert(k<20000);advance(s,h,end-gfs_rk4_context_time);
    }
}
int main(void){
    setenv("GFS_KEPLER_WALL_EVENTS","1",1);setenv("GFS_GPU_RAW_FACES","1",1);setenv("GFS_GPU_FINAL_VOLUME","1",1);
    static SimParameters s;s=initialize(.001,.1,0);run(&s,.02,.03);
    assert(hit_count==1 && fabs(hit_times[0]-.01)<1.e-10);
    assert(fabs(hypot(s.p->vx,s.p->vy)-.1)<1.e-13);free(s.p);
    double errors[3],a=2.75*cos(M_PI/1024.),x0=a-.001,v0=.1;
    double omega=sqrt(harmonic),amp=hypot(x0,v0/omega);
    double exact=(atan2(v0/omega,x0)-acos(a/amp))/omega;
    for(int k=0;k<3;k++){
        s=initialize(.001,.1,2);run(&s,.01/pow(2.,k),.03);
        assert(hit_count==1);errors[k]=fabs(hit_times[0]-exact);assert(errors[k]<1.e-7);free(s.p);
    }
    assert(errors[2]<=.2*errors[0]+3.e-12);
    /* Turn before the wall, despite an outward Euler predictor within band. */
    s=initialize(1.e-7,.001,1);acceleration=-10.;run(&s,.0003,.0003);
    assert(hit_count==0);free(s.p);acceleration=.2;
    /* Actual out-and-back event, with both endpoints of the initial trial inside. */
    s=initialize(1.e-8,.001,1);acceleration=-10.;run(&s,.0003,.0003);
    assert(hit_count==1 && fabs(hit_times[0]-(.001-sqrt(8.e-7))/10.)<1.e-9);
    free(s.p);acceleration=.2;
    s=initialize(.01,1.,0);run(&s,.7,12.);
    assert(hit_count==3 && fabs(hypot(s.p->vx,s.p->vy)-1.)<1.e-12);free(s.p);
    setenv("GFS_WALL_FORCE_REJECT_ONCE","1",1);
    s=initialize(.1,.1,0);run(&s,.001,.001);assert(restore_count==1 && hit_count==0);free(s.p);
    unsetenv("GFS_WALL_FORCE_REJECT_ONCE");
    setenv("GFS_WALL_FORCE_REJECT_ALWAYS","1",1);setenv("GFS_WALL_MAX_RETRIES","2",1);
    s=initialize(.1,.1,0);treevorostressrk4particletype before=*s.p;expect_abort=1;
    if(setjmp(failure_jump)==0){advance(&s,.001,.001);assert(0);}
    expect_abort=0;assert(saved_time==0. && memcmp(&saved_failure,&before,sizeof(before))==0);
    active=NULL;assert(gfs_rk4_context_time==0.);free(s.p);
    unsetenv("GFS_WALL_FORCE_REJECT_ALWAYS");unsetenv("GFS_WALL_MAX_RETRIES");
    s=initialize(.001,.1,0);
    s.p->x=3.+2.75-1.e-7;s.p->y=3.;s.p->vx=.1;s.p->vy=0.;before=*s.p;expect_abort=1;
    if(setjmp(failure_jump)==0){advance(&s,.001,.001);assert(0);}
    expect_abort=0;active=NULL;assert(saved_time==0. && memcmp(&saved_failure,&before,sizeof(before))==0);free(s.p);
    printf("PASS actual native event include: constant flight, accelerated errors %.3e %.3e %.3e, noncollision grazing, true out-and-back, repeated impacts, retry, corner and exact failure rollback\n",errors[0],errors[1],errors[2]);
    return 0;
}

"""Compile extracted native wall guard with synthetic invalid states."""
from pathlib import Path
import os
import subprocess
import tempfile

src=Path('/gpfs/kjhan/LagForce/Kepler_WallGuardOMP_20261029/project/Exam/gfs_wall_events.inc').read_text()
start=src.index('static void gfs_wall_new_generation(')
end=src.index('static int gfs_wall_endpoint_inside(',start)
extract=src[start:end]
prefix=r'''
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { double x,y,m; unsigned present,seen; } GfsWallBase;
typedef struct {int on; ptrdiff_t np; size_t extent; unsigned generation; GfsWallBase *base; const char *reason; double pressure[3],gravity[2],collision[3],roundoff,ceiling; int hits;} GfsWallEventStep;
typedef struct {size_t id; double mass,ie,volume,w2,x,y,vx,vy,ax,ay,die,U,V,anchor,rhs; int tag;} treevorostressrk4particletype;
typedef treevorostressrk4particletype treevorork4particletype;
typedef struct {ptrdiff_t np; treevorostressrk4particletype *p;} SimParameters;
#define PINDX(p) ((p)->id)
#define VORO_NP(s) ((s)->np)
#define VORORK4_TBP(s) ((s)->p)
#define GFS_T64_TAGGED(p) ((p).tag)
#define GFS_T64_U(p) ((p).U)
#define GFS_T64_V(p) ((p).V)
#define GFS_T64_ANCHOR(p) ((p).anchor)
#define GFS_T64_RHS(p) ((p).rhs)
static int thermal=1, geom_failed=0, gfs_wall_stage_failed=0, gfs_wall_runtime_active=0;
static double gfs_rk4_step_ceiling=0;
static int gfs_thermal64_requested(void){return thermal;}
static int gpu_geometry_trial_failed(void){return geom_failed;}
static void gpu_geometry_trial_begin(void){}
'''
suffix=r'''
int main(void){
    enum {N=256}; treevorostressrk4particletype p[N]; GfsWallBase base[N];
    SimParameters s={N,p}; GfsWallEventStep w={0};
    int failed=0;
    for(int mode=0;mode<2;mode++)for(int trial=0;trial<13;trial++){
        setenv("GFS_WALL_OMP_VALIDATE",mode?"1":"0",1);
        memset(p,0,sizeof(p));memset(base,0,sizeof(base));memset(&w,0,sizeof(w));
        w.on=1;w.np=N;w.extent=N;w.base=base;s.np=N;thermal=1;geom_failed=0;gfs_wall_stage_failed=0;
        for(int i=0;i<N;i++){
            p[i].id=(i*17)%N;p[i].mass=base[p[i].id].m=1.;base[p[i].id].present=1;
            p[i].ie=p[i].volume=p[i].U=p[i].V=p[i].tag=1.;
        }
        int rhs=1, expected=trial==0;
        if(trial==1)p[5].id=p[7].id;
        if(trial==2)p[5].id=N;
        if(trial==3)p[5].mass=2;
        if(trial==4)p[5].ie=0;
        if(trial==5)p[5].volume=NAN;
        if(trial==6)p[5].w2=1;
        if(trial==7)p[5].ax=NAN;
        if(trial==8)p[5].tag=0;
        if(trial==9)p[5].U=0;
        if(trial==10)p[5].rhs=NAN;
        if(trial==11)geom_failed=1;
        if(trial==12)s.np=N-1;
        int actual=gfs_wall_valid(&s,&w,rhs);
        if(actual!=expected){fprintf(stderr,"mode=%d trial=%d got=%d expected=%d\n",mode,trial,actual,expected);failed=1;}
    }
    if(failed)return 1;puts("PASS: serial and OMP guards reject 12 invalid-state classes");return 0;
}
'''
with tempfile.TemporaryDirectory(prefix='wallguard_omp_') as d:
    c=Path(d)/'test.c';exe=Path(d)/'test';c.write_text(prefix+extract+suffix)
    subprocess.run(['gcc','-O2','-std=gnu11','-fopenmp',str(c),'-lm','-o',str(exe)],check=True)
    for threads in ('2','8'):
        env=dict(os.environ,OMP_NUM_THREADS=threads,OMP_DYNAMIC='FALSE')
        subprocess.run([str(exe)],env=env,check=True)

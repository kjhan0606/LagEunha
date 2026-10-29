"""Compile the candidate's actual guard with minimal particle adapters."""
import json
from pathlib import Path
import subprocess
import tempfile

src=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64EOSGuard_20261029/project/Exam/gfs_wall_events.inc').read_text()
start=src.index('static int gfs_wall_t64_eos_valid(')
end=src.index('static int gfs_wall_endpoint_inside(',start)
helper=src[start:end]
code='''#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <float.h>
typedef struct { double u,v,mass,tag; } treevorostressrk4particletype;
typedef struct { treevorostressrk4particletype *p; double gamma; } SimParameters;
typedef struct { int on; ptrdiff_t np; const char *reason; } GfsWallEventStep;
#define VORORK4_TBP(s) ((s)->p)
#define GAS_GAMMA(s) ((s)->gamma)
#define GFS_T64_U(p) ((p).u)
#define GFS_T64_V(p) ((p).v)
#define GFS_T64_TAGGED(p) ((p).tag==42.)
static int enabled=1;
static int gfs_thermal64_requested(void){return enabled;}
''' + helper + '''
static int check(treevorostressrk4particletype p,int want){
 SimParameters s={&p,5./3.};GfsWallEventStep w={1,1,0};
 int got=gfs_wall_t64_eos_valid(&s,&w);
 return got==want && (want || w.reason!=0);
}
int main(void){
 treevorostressrk4particletype p={1.,1.,1.,42.};
 if(!check(p,1))return 1;
 p.u=DBL_MAX;p.v=.5;if(!check(p,0))return 2;
 p.u=1e100;p.v=1.;if(!check(p,0))return 3;
 p.u=1e-100;if(!check(p,0))return 4;
 p.u=1.;p.v=0.;if(!check(p,0))return 5;
 p.v=1.;p.tag=0.;if(!check(p,0))return 6;
 p.tag=42.;p.mass=1e-300;p.v=1e300;if(!check(p,0))return 7;
 p.mass=1.;p.v=1.;enabled=0;if(!check(p,1))return 8;
 puts("valid, double overflow, float overflow/underflow, zero volume, tag, density underflow and disabled mode passed");return 0;
}'''
with tempfile.TemporaryDirectory(prefix='eos-guard-') as d:
    path=Path(d);(path/'guard.c').write_text(code)
    subprocess.run(['cc','-std=c11','-O2','-Wall','-Werror',str(path/'guard.c'),'-lm','-o',str(path/'guard')],check=True)
    run=subprocess.run([str(path/'guard')],capture_output=True,text=True,check=True)
print(json.dumps({'passed':True,'result':run.stdout.strip(),
                  'scope':'Compiled actual guard with stub state/solver adapters; not a GPU trajectory or real rollback test.'},indent=2))

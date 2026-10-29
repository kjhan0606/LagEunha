"""Exercise the actual dispatch guard with a minimal non-CUDA host harness."""
import os
from pathlib import Path
import subprocess
import tempfile

src=Path(__file__).resolve().parent.parent/'Exam/gfs_thermal64_native.inc'
text=src.read_text()
begin=text.index('static void gfs_t64_dispatch_guard(')
end=text.index('static void gfs_t64_eos(',begin)
guard=text[begin:end]
harness='''#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct { int tagged; } treevorostressrk4particletype;
typedef struct { treevorostressrk4particletype p; } SimParameters;
#define VORORK4_TBP(s) (&(s)->p)
#define VORO_NP(s) (1)
#define GFS_T64_TAGGED(p) ((p).tagged)
#define MPI_COMM(s) (0)
static int gfs_thermal64_requested(void) { const char *s=getenv("GFS_THERMAL64"); return s && atoi(s)>0; }
static void gfs_es_config(void) {}
static void MPI_Abort(int comm,int code) { (void)comm; exit(code); }
'''+guard+'''
int main(void) { SimParameters s={{0}}; gfs_t64_dispatch_guard(&s); return 0; }
'''
with tempfile.TemporaryDirectory(prefix='lagforce_guard_') as temp:
    c=Path(temp)/'guard.c';binary=Path(temp)/'guard';c.write_text(harness)
    subprocess.run(['cc','-std=c99','-Wall','-Wextra','-Werror',str(c),'-o',str(binary)],check=True)
    def run(t64,moment):
        env=os.environ.copy();env['GFS_THERMAL64']=t64;env['GFS_MOMENT_INIT_ONLY']=moment
        return subprocess.run([str(binary)],env=env,capture_output=True,text=True)
    both=run('1','1');disabled=run('0','1');ordinary=run('1','0')
    assert both.returncode==99 and 'moment initialization cannot share' in both.stderr
    assert disabled.returncode==0
    assert ordinary.returncode==99 and 'unsupported or mixed payload contract' in ordinary.stderr
    print('Actual source guard: joint mode rejects, thermal-off moment path remains available, ordinary thermal dispatch continues: PASS')

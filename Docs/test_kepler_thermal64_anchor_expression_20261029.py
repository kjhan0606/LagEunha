"""Compile the actual extracted RK4 expression and check precision and legacy modes."""
import json
from pathlib import Path
import re
import subprocess
import tempfile

root = Path('/gpfs/kjhan/LagForce/Kepler_Thermal64AnchorFix_20261029')
src = (root/'project/Exam/exam.c').read_text()
hits = re.findall(r'double u_final=.*?;', src)
if len(hits) != 1 or 'thermal64?GFS_T64_ANCHOR' not in hits[0]:
    raise RuntimeError('Unexpected candidate RK4 expression')
code = '''#include <math.h>
#include <stdio.h>
#define GFS_T64_ANCHOR(p) ((p).anchor)
struct P { double anchor; float ie; struct { double E_inv_xx; } stress; };
static double eval(struct P *sbp,int thermal64,int direct_u,double u_increment){int i=0;
''' + hits[0] + '''
return u_final; }
int main(void){
 struct P p={1.00000003, 1.0f, {2.0}};
 double a=eval(&p,1,0,1e-8), legacy=eval(&p,0,0,1e-8), direct=eval(&p,0,1,1e-8);
 if (a != p.anchor+1e-8 || legacy != (double)p.ie+1e-8 || direct != p.stress.E_inv_xx+1e-8 || !(a>legacy)) return 1;
 printf("anchor %.17g legacy %.17g direct %.17g\\n",a,legacy,direct);
 return 0;
}'''
with tempfile.TemporaryDirectory(prefix='thermal64-anchor-') as d:
    path = Path(d)
    (path/'test.c').write_text(code)
    subprocess.run(['cc','-std=c11','-O2','-Wall','-Werror',str(path/'test.c'),'-o',str(path/'test')],check=True)
    out = subprocess.run([str(path/'test')],text=True,capture_output=True,check=True).stdout.strip()
print(json.dumps({'passed': True, 'extracted_expression': hits[0], 'output': out,
                  'scope': 'Standalone compiled arithmetic regression, not a native evolution or wall-event test.'}, indent=2))

"""Patch an isolated v3 source tree to handle a zero pressure gradient."""

import hashlib
import json
from pathlib import Path
import shutil
import sys


EXPECTED="85dbe967c2d289380711dad0bfed2ebf43d4cbe7b12739af29c9ffef0e3111c6"
HEADER="gfs_pressure_neutral_projection_20261027.h"
OLD='''    {
        long double dot=0.,norm=0.;
        for(int i=0;i<n_particles;++i) {
            if(!host_active[i]) continue;
            if(!isfinite(host_gradient[i]) || !isfinite(host_rate[i])) {
                status=-2; goto done;
            }
            dot+=(long double)host_gradient[i]*host_rate[i];
            norm+=(long double)host_gradient[i]*host_gradient[i];
        }
        if(!(norm>0.)) { status=-3; goto done; }
        double alpha=(double)(dot/norm);
        for(int i=0;i<n_particles;++i) {
            if(host_active[i]) host_rate[i]-=alpha*host_gradient[i];
        }
        *host_scalar=alpha;
    }
'''
NEW='''    status=gfs_pressure_neutral_projection(host_rate,host_gradient,
        host_active,n_particles,host_scalar);
    if(status) goto done;
'''


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if len(sys.argv)!=2:
        raise SystemExit('Expected isolated source root')
    root=Path(sys.argv[1]).resolve()
    source=root/'Exam/gfs_energy_neutral_probe.cu'
    makefile=root/'Exam/Makefile'
    header=root/'Exam'/HEADER
    origin=Path(__file__).resolve().parent/HEADER
    if digest(source)!=EXPECTED or header.exists():
        raise ValueError('Unexpected source hash or preexisting projection header')
    text=source.read_text()
    if text.count(OLD)!=1 or text.count('#include <float.h>')!=1:
        raise ValueError('Projection source anchor changed')
    rule='gfs_energy_neutral_probe.o: gfs_energy_neutral_probe.cu kepler_gpu_jacobi_probe_20261027.cu exam_gpu.h'
    make=makefile.read_text()
    if make.count(rule)!=1:
        raise ValueError('Projection Makefile dependency anchor changed')
    shutil.copy2(origin,header)
    source.write_text(text.replace('#include <float.h>',
        '#include <float.h>\n#include "'+HEADER+'"').replace(OLD,NEW))
    makefile.write_text(make.replace(rule,rule+' '+HEADER))
    print(json.dumps({'parent_probe_sha256':EXPECTED,'probe_sha256':digest(source),
                      'header_sha256':digest(header),'makefile_sha256':digest(makefile),
                      'scope':'Exact-zero gradient leaves rates unchanged; nonzero '
                              'projection retains legacy accumulation order.'},indent=2))


if __name__=='__main__':
    main()

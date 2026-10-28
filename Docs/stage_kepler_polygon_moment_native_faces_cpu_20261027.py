"""Stage polygon-moment comparisons on saved W0 and weighted native states."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/polygon_moment_native_faces_cpu_v1')


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if RUN.exists():raise FileExistsError(RUN)
    sources={'job.slurm':DOCS/'kepler_polygon_moment_native_faces_cpu_20261027.slurm',
             'w0_jitter_t65.json':DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json',
             'weighted_jitter_t10.json':DOCS/'KEPLER_Q2_WEIGHTED_T10_EQUILIBRIUM_2026-10-27/jitter_corrected_audit.json'}
    for name in ('gfs_polygon_moment_20261027.h','gfs_polygon_moment_cpu_bridge_20261027.cpp',
                 'audit_kepler_polygon_moment_native_faces_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        sources[name]=DOCS/name
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir()
    for name,source in sources.items():
        shutil.copy2(source,RUN/name)
        if digest(RUN/name)!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as f:
        for name in sorted(sources):f.write(f'{digest(RUN/name)}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as f:
        json.dump({'run':str(RUN),'inputs':{name:{'source':str(sources[name]),'sha256':digest(RUN/name)}
                                           for name in sorted(sources)},
                   'scope':'CPU comparison of CUDA-compatible moments with independent '
                           'triangle integration on reconstructed native checkpoint faces.'},f,indent=2)
        f.write('\n')
    print(json.dumps({'run':str(RUN),'input_count':len(sources)}))


if __name__=='__main__':main()

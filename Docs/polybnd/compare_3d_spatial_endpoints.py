"""Independent-resolution saved-state comparison at equal physical time."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from dense_record_3d import load
from power_cells_3d import Geometry3D
from check_3d_vortex_evolution import score,atomic


def run(paths,library):
    rows=[];time=None;amplitude=None
    for path in map(Path,paths):
        cfg=json.loads((path.parent/'configuration.json').read_text())
        initial=json.loads((path.parent/'initial_state.json').read_text())
        record,endpoint=load(path);t=record['t']+record['dt']
        if time is not None and abs(time-t)>1e-15:raise ValueError('different comparison times')
        if amplitude is not None and amplitude!=cfg['amplitude']:raise ValueError('different reference amplitudes')
        time=t;amplitude=cfg['amplitude']
        np.testing.assert_array_equal(initial['mass'],endpoint['raw'][:,0])
        rows.append(dict(n=cfg['n'],dt=record['dt'],input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            configuration=cfg,initial_fields=json.loads((path.parent/'initial_score.json').read_text()),
            evolved_fields=score(Geometry3D(library),endpoint['sites'],endpoint['weights'],endpoint['raw'],amplitude)))
    report=dict(t=time,rows=rows,full_validation=False,scope='Equal-time saved-endpoint spatial comparison; no evolution or cross-mesh array subtraction.')
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');p.add_argument('paths',nargs='+');a=p.parse_args();run(a.paths,a.library)

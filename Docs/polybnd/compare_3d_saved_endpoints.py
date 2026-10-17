"""Read-only equal-time physical field comparison of saved 3D steps."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from dense_record_3d import load
from power_cells_3d import Geometry3D
from check_3d_vortex_evolution import score,atomic


def run(coarse,fine,library,amplitude=.3):
    records=[load(p) for p in (coarse,fine)]
    times=[r['t']+r['dt'] for r,e in records]
    if abs(times[0]-times[1])>1e-15:raise ValueError('comparison requires equal physical time')
    a,b=[e for r,e in records];np.testing.assert_array_equal(a['raw'][:,0],b['raw'][:,0])
    geo=Geometry3D(library)
    scores=[]
    for label,(record,endpoint),path in zip(('coarse','fine'),records,(coarse,fine)):
        scores.append(dict(label=label,dt=record['dt'],input_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            fields=score(geo,endpoint['sites'],endpoint['weights'],endpoint['raw'],amplitude)))
    report=dict(t=times[0],scores=scores,maximum_raw_difference=float(np.max(abs(a['raw']-b['raw']))),
        maximum_site_difference=float(np.max(abs(a['sites']-b['sites']))),full_validation=False)
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('coarse');p.add_argument('fine');p.add_argument('library');a=p.parse_args();run(a.coarse,a.fine,a.library)

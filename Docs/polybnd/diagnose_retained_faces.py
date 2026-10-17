"""Check every saved failed capacity face with analytic interval clipping."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from saved_curved_contact import checked_json,initialize_saved
from curved_transport_reference import numerical_affine_flow
from transported_density_capacity import PolygonDensity,capacity_edges
from active_power_cells import power_cells_active
import moment_ale_evolve as model


def run():
    cfg=json.loads(Path('curved_transport_config.json').read_text())
    ic={k:cfg[k] for k in ('initial_state','initial_state_sha256','reference','reference_sha256')}
    ic.update(kind='saved_curved_contact',n=cfg['n']);mass,_=initialize_saved(ic)
    source=checked_json(cfg['source'],cfg['source_sha256']);a=np.asarray(source['state'])[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],source['t'])
    b,B,P,diag=numerical_affine_flow(source['state'],mass,np.asarray([Polygon(v).area for v in cells]))
    saved=json.loads(Path('t0.01_pending.json').read_text());t=saved['t'];dt=t-source['t'];F=np.eye(2)+dt*B;shift=dt*b
    data=np.load('retained_field.npz');field=PolygonDensity(data['vertices'],data['densities']).affine_image(F,shift)
    sites=a[:,:2]@F.T+shift;domain=model.domain(source['t'])@F.T+shift;weights=np.asarray(saved['weights'])
    cells=power_cells_active(sites,weights,domain);failures=[];max_gap=0.;count=0
    for i,j,start,end in capacity_edges(sites,weights,cells):
        count+=1
        try:
            integral,gap=field.line_integral(start,end);max_gap=max(max_gap,gap)
        except Exception as exc:failures.append(dict(i=i,j=j,error=repr(exc)))
    result=dict(passed=not failures,faces=count,max_gap=max_gap,failures=failures,jacobian_enclosures=field.jacobian_enclosures,
        note='All faces at the preserved failed trial. Local long-double halfplane intervals; unchanged coverage guard and exact overlap mass acceptance.')
    atomic('retained_face_diagnosis.json',result);print(result,flush=True)
    raise SystemExit(0 if result['passed'] else 2)


if __name__=='__main__':run()

"""Compare actual late native event trajectories at three timestep limits."""
import argparse
import json
from pathlib import Path
import re
import numpy as np
from audit_kepler_native_radial_20261026 import ordered


def analyze(root):
    paths=[root.parent/'crossing_to484',root/'half',root/'quarter']
    audits=[json.loads((r/'event_audit.json').read_text()) for r in paths]
    states=[ordered(Path(a['endpoint']))[0] for a in audits]
    times=[a['impacts'][0]['time'] for a in audits]
    time_scale=sum(a['impacts'][0]['tolerance'] for a in audits)/min(a['impacts'][0]['vn_before'] for a in audits)
    time_scale+=4*float(np.spacing(max(times)))
    fields=['x','y','vx','vy','volume','den','ie','pressure']
    differences=[]
    for i in [0,1]:
        a,b=states[i],states[2];weights=b['mass'].astype(float)
        rows={}
        for name in fields:
            delta=a[name].astype(float)-b[name].astype(float)
            scale=max(float(np.sqrt(np.sum(weights*b[name].astype(float)**2)/sum(weights))),1.e-300)
            rows[name]={'max_absolute':float(np.max(abs(delta))),
                        'normalized_mass_RMS':float(np.sqrt(np.sum(weights*delta**2)/sum(weights))/scale)}
        differences.append({'case':str(paths[i]),'reference':str(paths[2]),'fields':rows})
    impulses=[]
    for r in paths:
        rows=[line for line in (r/'log').read_text(errors='replace').splitlines() if '[WALL_ACCEPT]' in line]
        impulse=np.zeros(3)
        for line in rows:
            match=re.search(r'wall_collision_J=\(([^,]+),([^,]+),([^\)]+)\)',line)
            impulse+=np.array([float(match[k]) for k in [1,2,3]])
        impulses.append(impulse.tolist())
    checks={'all_individual_predeclared_gates_passed':all(a['passed'] for a in audits),
            'one_collision_same_owner_and_facet':all(len(a['impacts'])==1 and (a['impacts'][0]['id'],a['impacts'][0]['side'])==(55283,543) for a in audits),
            'same_final_time':all(a['time']==48.4 for a in audits),
            'same_fixed_masses':all(np.array_equal(states[0]['mass'],p['mass']) for p in states[1:]),
            'event_time_spread_within_contact_resolution':max(times)-min(times)<=time_scale}
    return {'passed':all(checks.values()),'checks':checks,'cases':[str(p) for p in paths],
            'event_times':times,'event_time_spread':max(times)-min(times),
            'combined_spatial_to_time_resolution':time_scale,
            'restart_relative_energy_errors':[a['restart_relative_total_energy_change'] for a in audits],
            'maximum_angular_balance_residuals':[a['maximum_angular_residual'] for a in audits],
            'collision_impulses_wall_Jx_Jy_torque':impulses,'endpoint_field_comparisons':differences,
            'scope':'Three actual native late-state trajectories. The contact-resolution comparison is a reported numerical scale, not a rigorous continuum error bound. Endpoint field differences are reported without fitted pass tolerances. No full fourth-order event convergence, Q1 accuracy or long-time stability claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=analyze(args.root)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r,indent=2))
    if not r['passed']:raise SystemExit(1)

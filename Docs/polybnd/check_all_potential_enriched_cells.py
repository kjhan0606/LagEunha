"""All-active-cell fixed-geometry potential enrichment representability audit."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from adaptive_coupled_mass import restore
from gravity_coupled_2d import hydro_dense_state
from audit_kepler_energy_step import configure,phi,atomic
from audit_density_potential_dual import basis
from kepler_disk_reference import polar_rule
from check_potential_enriched_density import check


def run(seed,tangent):
    root=Path(seed);tangent=Path(tangent)
    names=('source.json','trial.json','configuration.json')
    data={n:(root/n).read_bytes() for n in names}
    data['enrichment_cells.json']=(tangent/'enrichment_cells.json').read_bytes()
    identity={n:hashlib.sha256(v).hexdigest() for n,v in data.items()}
    identity['source_manifest_sha256']=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest()
    prior=json.loads((tangent/'report.json').read_text())['identity']
    if any(identity[n]!=prior[n] for n in names):raise ValueError('tangent inputs differ')
    configure(json.loads(data['configuration.json']))
    mass=np.asarray(json.loads(data['source.json'])['mass'])
    record=restore(json.loads(data['trial.json'])['record']);y,L,dy=hydro_dense_state(record,mass,0.)
    a=y[:-4].reshape(-1,20);mask=np.asarray(record['frozen_first'],dtype=bool)
    tangent_data=json.loads(data['enrichment_cells.json'])
    np.testing.assert_array_equal(mass,tangent_data['mass']);np.testing.assert_array_equal(mask,tangent_data['frozen'])
    cells,_=base.geometry(a[:,:2],a[:,2],record['t']);active=np.flatnonzero(~mask)
    saved=dict(identity=identity,next_index=0,rows=[])
    if Path('progress.json').exists():
        saved=json.loads(Path('progress.json').read_text())
        if saved['identity']!=identity:raise ValueError('audit cache identity changed')
    for index in range(saved['next_index'],len(active)):
        i=int(active[index]);raw=np.r_[mass[i],a[i,3:]];density=base.decode(cells[i],raw)['density']
        x,w=polar_rule(cells[i],16,radii=(.25,.5));dm=w*density.value(x);b=basis(x,density)
        coefficient=np.linalg.solve(b.T@((dm/mass[i])[:,None]*b),b.T@(dm*phi(x)/mass[i]))
        remainder=phi(x)-b@coefficient;sigma=float(np.sqrt(dm@(remainder**2)/mass[i]))
        fixture=dict(cell=i,polygon=cells[i],raw=raw,coefficient=coefficient,sigma=sigma,
            potential_mean=float(dm@phi(x)/mass[i]),missing_rate=tangent_data['missing'][i],
            density_quadrature_order=density.quadrature_order)
        try:result=check(fixture,steps=(.00025,))
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            result=dict(cell=i,rows=[dict(converged=False,error=str(exc))])
        if not all(r['converged'] for r in result['rows']):atomic(f'failed_cell_{i}.json',fixture)
        saved['rows'].append(result);saved['next_index']=index+1
        if (index+1)%16==0 or index+1==len(active):
            atomic('progress.json',saved);print('enriched cells',index+1,'/',len(active),flush=True)
    failures=[row['cell'] for row in saved['rows'] if not all(r['converged'] for r in row['rows'])]
    checks=[c for row in saved['rows'] for r in row['rows'] if r['converged'] for c in r['checks']]
    report=dict(identity=identity,cells=len(active),failures=failures,
        max_mass_relative=max((abs(c['mass_relative']) for c in checks),default=None),
        max_retained_moment_error=max((c['retained_moment_error'] for c in checks),default=None),
        max_potential_mean_error=max((abs(c['potential_mean_error']) for c in checks),default=None),
        full_validation=False,scope='all active cells, fixed geometry and m,q,Q; not coupled evolution')
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');p.add_argument('tangent');a=p.parse_args();run(a.seed,a.tangent)

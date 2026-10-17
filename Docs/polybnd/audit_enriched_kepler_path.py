"""Independent mass and potential path audit of a saved enriched hydro step."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import restore
from potential_moment_coupled_2d import dense_state,provider
from audit_kepler_energy_step import configure,phi,atomic
from kepler_disk_reference import polar_rule,frozen
from smooth_mass_path import evaluate_smooth_path
from check_kepler_dense_step import simpson


def run(seed,configuration):
    data=Path(seed).read_bytes();cfg=Path(configuration).read_bytes()
    identity=dict(trial_sha256=hashlib.sha256(data).hexdigest(),configuration_sha256=hashlib.sha256(cfg).hexdigest(),
        source_manifest_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest())
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('audit identity changed')
    else:atomic('identity.json',identity)
    configure(json.loads(cfg));trial=json.loads(data);record=restore(trial['hydro_record']);mass=np.asarray(trial['mass'])
    mask=np.asarray(record['frozen_first'],dtype=bool);quad=lambda p:polar_rule(p,16,radii=(.25,.5))
    last,_,_,lastU,_=dense_state(record,mass,1.)
    np.testing.assert_allclose(last[:-4].reshape(-1,20)[:,:9],np.asarray(trial['state'])[:-4].reshape(-1,20)[:,:9],atol=2e-14,rtol=0)
    np.testing.assert_allclose(lastU,trial['potential_moments'],atol=2e-14,rtol=0)
    rows=json.loads(Path('nodes.json').read_text()) if Path('nodes.json').exists() else []
    for k in range(len(rows),17):
        s=k/16;y,L,dy,U,dU=dense_state(record,mass,s);a=dy[:-4].reshape(-1,20)
        np.testing.assert_array_equal(frozen(record['t']+s*record['dt'],y[:-4].reshape(-1,20)[:,:2]),mask)
        flux,d=evaluate_smooth_path(record['t']+s*record['dt'],y,L,mass,(a[:,:2],a[:,2]),
            frozen_mask=mask,potential=phi,density_provider=provider(U,phi,quad,mask))
        rows.append(dict(flux=flux,generator=d['generator_defect'],potential_boundary_rate=d['active_potential_mass_exchange']))
        atomic('nodes.json',rows);print('enriched path node',k,'/16',flush=True)
    rates=np.asarray([r['flux'] for r in rows]);dt=record['dt']
    absolute=simpson(abs(rates),dt);signed=simpson(rates,dt)
    report=dict(identity=identity,dt=dt,
        absolute_relative_mass=float(np.max(absolute[~mask]/mass[~mask])),
        signed_relative_mass=float(np.max(abs(signed[~mask])/mass[~mask])),
        nested_difference=float(np.max(abs(signed-simpson(rates[::2],dt))[~mask]/mass[~mask])),
        maximum_generator_defect=max(r['generator'] for r in rows),
        integrated_potential_boundary=float(simpson(np.asarray([r['potential_boundary_rate'] for r in rows]),dt)),
        full_validation=False,scope='single saved hydro subflow; no root solve at audit nodes')
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');p.add_argument('configuration');a=p.parse_args();run(a.seed,a.configuration)

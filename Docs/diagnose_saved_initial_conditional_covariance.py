"""Read the saved initial quadratic split; no fits, fluxes or evolution.

Quantify the covariance still represented by the four residual samples.
This is a scale diagnostic, not a new closure or angular counterfactual.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with np.load(p) as f:return {k:f[k] for k in f.files}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    b=Path('/gpfs/kjhan/LagForce');icfile=b/'Codex_full_curved_kepler_IC_20261019/initial.npz';initial=read(icfile)
    folder=b/'Codex_full_disk_initial_velocity_gravity_20261019';prepfile=folder/'prepared.npz';prep=read(prepfile)
    report=json.loads((folder/'report.json').read_text());assert report['status']=='passed_initial_velocity_gravity_preparation' and report['initial_sha256']==sha(icfile)
    assert all(np.array_equal(prep[key],initial[key]) for key in ('ids','mass','state','factors','Psi'))
    ids=initial['ids'];mass=initial['mass'];n=len(ids);raw=np.column_stack([mass,initial['state'][:-4].reshape(n,20)[:,3:]])
    T=raw[:,17]-.5*(raw[:,13]+raw[:,16]);assert (T>0).all()
    H=initial['factors']@prep['residual_map'];sing=np.linalg.svd(H,compute_uv=False)
    residual=.5*np.sum(sing*sing,axis=1);total=.5*np.sum(initial['factors']**2,axis=(1,2));assert (total>0).all()
    physical_map=H/np.sqrt(mass)[:,None,None]
    fields=dict(conditional_residual_energy=residual,conditional_fraction_of_total_R=residual/total,
        conditional_energy_over_thermal=residual/T,conditional_maximum_velocity_dispersion=sing[:,0]/np.sqrt(mass),
        physical_residual_map_asymmetry=np.max(abs(physical_map-physical_map.transpose(0,2,1)),axis=(1,2)))
    cavity_file=b/'Codex_full_disk_checkpoint_error_diagnostics_20261019/saved_density_volume_drift/initial_density_volume.npz'
    cavity=read(cavity_file);assert np.array_equal(ids,cavity['ids'])
    radius=np.linalg.norm(initial['state'][:-4].reshape(n,20)[:,:2],axis=1)
    regions=dict(all=np.ones(n,bool),inner_original_sites_lt_0p4=radius<.4,exact_initial_uniform_cavity=cavity['initial_uniform_cavity'])
    rows={}
    for name,mask in regions.items():
        w=mass[mask].astype(np.longdouble);terms={}
        for key,value in fields.items():
            x=value[mask].astype(np.longdouble);j=int(np.argmax(abs(x)))
            terms[key]=dict(mass_weighted_mean=float(np.sum(w*x)/w.sum()),mass_weighted_RMS=float(np.sqrt(np.sum(w*x*x)/w.sum())),unweighted_median=float(np.median(x)),maximum=float(x[j]),maximum_id=int(ids[mask][j]))
        rows[name]=dict(cells=int(mask.sum()),terms=terms)
    np.savez(args.output/'conditional_covariance.npz',ids=ids,**fields,**regions)
    result=dict(status='completed_saved_initial_conditional_covariance_scale',time=0.,regions=rows,
        source_sha256={str(f):sha(f) for f in (icfile,prepfile,folder/'report.json',cavity_file)},
        physical_providers=0,density_fits=0,velocity_fits=0,geometry_calls=0,accepted_time_steps=0,
        limits='Initial saved conditional covariance scale only, not its evolved value, actual face-flux error, angular invariance certification or a causal attribution. No alternative sample/closure evaluated.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],regions=rows)),flush=True)


if __name__=='__main__':main()

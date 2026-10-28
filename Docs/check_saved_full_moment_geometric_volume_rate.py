"""Compare actual captured face wn integral with exact power-plane area rate."""
import argparse,json
from pathlib import Path
import numpy as np


def read(p):
    with np.load(p) as f:return {k:f[k] for k in f.files}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--capture',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists();base=Path('/gpfs/kjhan/LagForce')
    cp=base/'Codex_full_disk_continuous_two_dt_to_t1_20261019/coarse/checkpoints/coarse_0005'
    state=read(cp/'accepted.npz');mesh=read(cp/'geometry/mesh.npz');stage=read(args.capture/'stage_compact.npz');captured=read(args.capture/'signed_flux_volume.npz')
    ld=np.longdouble;ids=state['ids'];n=len(ids);lookup={int(k):i for i,k in enumerate(ids)}
    assert all(np.array_equal(ids,d['ids']) for d in (mesh,stage,captured))
    x=state['state'][:-4].reshape(n,20)[:,:2].astype(ld);g=stage['generator'].astype(ld);z=stage['weight_rates'].astype(ld)
    i,j=np.array([(lookup[int(a)],lookup[int(b)]) for a,b in mesh['face_keys']]).T
    f=mesh['faces'].astype(ld);mid,length,normal=f[:,:2],f[:,4],f[:,5:7]
    delta=x[j]-x[i];dg=g[j]-g[i];distance=np.sqrt(np.sum(delta*delta,axis=1))
    normal_error=float(abs(delta/distance[:,None]-normal).max())
    wn=(.5*np.sum((x[j]+x[i]-2*mid)*dg,axis=1)+.5*np.sum(delta*(g[j]+g[i]),axis=1)-.5*(z[j]-z[i]))/distance
    area_rate=np.zeros(n,ld);np.add.at(area_rate,i,length*wn);np.add.at(area_rate,j,-length*wn)
    difference=captured['Vdot']-area_rate;area=mesh['moments'][:,0].astype(ld)
    report=dict(status='completed_saved_geometric_volume_rate_crosscheck',time=.025,
        normal_error_maximum=normal_error,exact_power_plane_total_volume_rate=float(area_rate.sum()),
        captured_total_volume_rate=float(captured['Vdot'].sum()),
        maximum_cell_volume_rate_difference=float(abs(difference).max()),
        maximum_fractional_volume_rate_difference=float(abs(difference/area).max()),
        maximum_fractional_volume_rate=float(abs(area_rate/area).max()),
        physical_providers=0,accepted_time_steps=0,thresholds_added=False,
        limits='Independent saved straight-face affine speed integral, stationary true circular wall. This checks instantaneous Vdot only, not finite-time volume-map consistency.')
    args.output.mkdir();np.savez(args.output/'volume_rates.npz',ids=ids,exact_plane_Vdot=area_rate,captured_Vdot=captured['Vdot'],difference=difference)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()

"""Compare the captured instantaneous tangent to saved later error changes.

No new RHS. A one-time tangent cannot attribute accumulated finite-step errors.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from compare_full_moment_checkpoint_errors import load,REGIONS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--analysis',type=Path,required=True)
    p.add_argument('--initial-cavity',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists()
    root=Path('/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019/coarse/checkpoints')
    a,ra,pa=load(root/'coarse_0005/physics');b,rb,pb=load(root/'coarse_0010/physics')
    assert np.array_equal(a['ids'],b['ids']) and np.array_equal(a['actual_mass'],b['actual_mass'])
    with np.load(args.analysis/'signed_rates.npz') as f:rates={k:f[k] for k in f.files}
    with np.load(args.initial_cavity) as f:uniform=f['initial_uniform_cavity'];assert np.array_equal(f['ids'],a['ids'])
    assert np.array_equal(rates['ids'],a['ids']);dt=rb['time']-ra['time'];assert dt==.025
    mass=a['actual_mass'].astype(np.longdouble)
    def metrics(observed,tangent,mask):
        x=observed[mask].astype(np.longdouble);y=(dt*tangent[mask]).astype(np.longdouble);w=mass[mask]
        rms=lambda z:float(np.sqrt(np.sum(w*z*z)/w.sum()))
        denom=np.sqrt(np.sum(w*x*x)*np.sum(w*y*y))
        return dict(cells=int(mask.sum()),observed_change_RMS=rms(x),tangent_prediction_RMS=rms(y),
                    residual_RMS=rms(x-y),alignment=float(np.sum(w*x*y)/denom) if denom>0 else None,
                    signed_mean_observed=float(np.sum(w*x)/w.sum()),signed_mean_predicted=float(np.sum(w*y)/w.sum()))
    pressure={}
    for name in REGIONS:
        mask=a[name]&b[name]
        if not mask.any():continue
        pressure[name]={k:metrics(b['relative_pressure_error']-a['relative_pressure_error'],rates[k],mask)
                        for k in ('relative_pressure_volume','relative_pressure_raw_thermal','relative_pressure_factor_minus_raw','relative_pressure_factor_total')}
    mask=uniform&a['cavity']&b['cavity']
    # Only here is the reference density a fixed floor in every compared cell.
    density_rate=(1+a['relative_density_reference_error'])*rates['cell_mean_density_fractional_rate']
    density=metrics(b['relative_density_reference_error']-a['relative_density_reference_error'],density_rate,mask)
    report=dict(status='completed_saved_tangent_vs_later_error_growth',start=.025,end=.05,
        source_A=pa,source_B=pb,pressure=pressure,uniform_cavity_density=density,
        rates_sha256=hashlib.sha256((args.analysis/'signed_rates.npz').read_bytes()).hexdigest(),
        physical_providers=0,accepted_time_steps=0,
        limits='Frozen instantaneous tangent prediction, not integrated along the trajectory. Common masks only. Pressure components can cancel, so individual correlations are not unique causes. No convergence or long-time extrapolation.')
    args.output.mkdir();(args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],density=density,inner_contact_pressure=pressure['inner_contact'])),flush=True)


if __name__=='__main__':main()

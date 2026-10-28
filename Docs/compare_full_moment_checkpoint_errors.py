"""Compare errors relative to each checkpoint's own actual-cell steady field.

Read-only diagnostics, no RHS or evolution. Primary differences use the same
mass weighting and wrapped phase convention as the approved two-dt monitor.
No new acceptance tolerance is introduced by this script.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

REGIONS=('all','inner','cavity','disk_core','inner_contact','outer_contact','dynamic_buffer','wall')
FIELDS=('mass','density','centroid','velocity','radial_velocity','azimuthal_velocity',
        'intrinsic_ell','origin_L','R','thermal','specific_thermal','pressure','Psi')


def sha(file):
    h=hashlib.sha256()
    with Path(file).open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()


def load(folder,allow_failed_monitor=False):
    folder=Path(folder).resolve()
    report=json.loads((folder/'report.json').read_text())
    assert report['status']=='completed_physical_state_comparison'
    gate=json.loads((folder.parent/'physical_gate.json').read_text())
    passed=gate['status']=='passed_checkpoint_physical_monitor'
    assert passed or (allow_failed_monitor and gate['status']=='failed_checkpoint_physical_monitor')
    assert gate['time']==report['time'] and gate['regions']
    assert not passed or all(row['passed'] for row in gate['regions'].values())
    density=json.loads((folder.parent/'density/report.json').read_text())
    reference=json.loads((folder.parent/'reference/report.json').read_text())
    assert density['status']=='passed_full_stage_Psi_density_preparation'
    assert reference['status']=='passed_steady_reference_on_moving_cells'
    assert density['state_sha256']==report['state_sha256']
    assert reference['actual_state_sha256']==report['state_sha256']
    assert sha(reference['actual_state_path'])==report['state_sha256']
    assert reference['time']==report['time']
    file=folder/'fields.npz'
    with np.load(file) as f:data={key:f[key] for key in f.files}
    n=len(data['ids']);assert n and np.all(np.diff(data['ids'])>0)
    assert data['actual_mass'].shape==(n,) and (data['actual_mass']>0).all()
    assert all(np.isfinite(value).all() for value in data.values())
    for key in REGIONS:assert data[key].shape==(n,) and data[key].dtype==np.dtype(bool)
    assert data['all'].all() and data['phase_defined'].dtype==np.dtype(bool)
    for key in FIELDS:assert data['actual_'+key].shape==data['reference_'+key].shape
    return data,report,dict(folder=str(folder),fields_sha256=sha(file),
        physical_monitor_passed=passed,physical_monitor_status=gate['status'],
        physics_sha256=hashlib.sha256(json.dumps(reference['physics'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        physics_report_sha256=sha(folder/'report.json'),
        physical_gate_sha256=sha(folder.parent/'physical_gate.json'),
        density_report_sha256=sha(folder.parent/'density/report.json'),
        reference_report_sha256=sha(folder.parent/'reference/report.json'))


def compare(a,b):
    assert np.array_equal(a['ids'],b['ids']) and np.array_equal(a['actual_mass'],b['actual_mass'])
    mass=a['actual_mass'].astype(np.longdouble);regions={}
    for name in REGIONS:
        mask=a[name]|b[name]
        if not mask.any():continue
        weight=mass[mask]
        def rms(value):
            value=np.asarray(value,dtype=np.longdouble).reshape(len(weight),-1)
            return float(np.sqrt(np.sum(weight*np.sum(value*value,axis=1))/weight.sum()))
        pressure=b['relative_pressure_error'][mask]-a['relative_pressure_error'][mask]
        density=b['relative_density_reference_error'][mask]-a['relative_density_reference_error'][mask]
        fields={}
        for key in FIELDS:
            ea=a['actual_'+key][mask].astype(np.longdouble)-a['reference_'+key][mask].astype(np.longdouble)
            eb=b['actual_'+key][mask].astype(np.longdouble)-b['reference_'+key][mask].astype(np.longdouble)
            difference=eb-ea;reference_rms=rms(a['reference_'+key][mask])
            fields[key]=dict(error_difference_RMS=rms(difference),
                maximum_absolute_error_difference=float(abs(difference).max()),
                reference_A_RMS=reference_rms,
                error_difference_relative_RMS=rms(difference)/reference_rms if reference_rms>0 else None)
        valid=mask&a['phase_defined']&b['phase_defined']
        phase=None
        if valid.any():
            d=b['velocity_direction_error_degrees'][valid]-a['velocity_direction_error_degrees'][valid]
            d=(d+180)%360-180
            phase=dict(cells=int(valid.sum()),RMS=float(np.sqrt(np.average(d*d,weights=mass[valid]))),maximum=float(abs(d).max()))
        regions[name]=dict(cells=int(mask.sum()),cells_A=int(a[name].sum()),cells_B=int(b[name].sum()),
            intersection_cells=int((a[name]&b[name]).sum()),
            only_A_ids=a['ids'][a[name]&~b[name]].tolist(),
            only_B_ids=b['ids'][b[name]&~a[name]].tolist(),
            phase_excluded_cells=int((mask&~valid).sum()),
            relative_pressure_error_difference_RMS=rms(pressure),
            relative_pressure_error_difference_maximum=float(abs(pressure).max()),
            relative_density_error_difference_RMS=rms(density),
            velocity_error_difference_relative_RMS=fields['velocity']['error_difference_relative_RMS'],
            phase_error_difference_degrees=phase,fields=fields)
    return regions


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--a',type=Path,required=True);p.add_argument('--b',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--require-identical',action='store_true',help='Diagnostic self-comparison sanity check only')
    p.add_argument('--include-failed-monitor',action='store_true',
        help='Read failed monitor candidates for diagnostics; never treat them as accepted states')
    args=p.parse_args()
    a,ra,pa=load(args.a,args.include_failed_monitor);b,rb,pb=load(args.b,args.include_failed_monitor)
    assert ra['time']==rb['time']
    assert pa['physics_sha256']==pb['physics_sha256']
    regions=compare(a,b)
    if args.require_identical:
        assert all(np.array_equal(a[key],b[key]) for key in a)
        for row in regions.values():
            assert row['relative_pressure_error_difference_RMS']==row['relative_density_error_difference_RMS']==0
            assert all(field['error_difference_RMS']==0 for field in row['fields'].values())
            assert row['phase_error_difference_degrees'] is None or row['phase_error_difference_degrees']['RMS']==0
    report=dict(status='completed_checkpoint_actual_cell_error_difference_diagnostics',time=ra['time'],
        source_A=pa,source_B=pb,regions=regions,
        existing_physical_monitors_both_passed=pa['physical_monitor_passed'] and pb['physical_monitor_passed'],
        failed_monitor_diagnostics_requested=args.include_failed_monitor,
        self_comparison_sanity_checked=args.require_identical,thresholds_added=False,
        physical_providers=0,accepted_time_steps=0,
        limits='Diagnostic only. Each reference is projected onto its own actual moving cell; same IDs do not imply identical cell geometry. No cross-host roundoff bound, convergence order or long-time stability is certified.')
    args.output.mkdir(exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],time=report['time'],all=regions['all'])),flush=True)


if __name__=='__main__':main()

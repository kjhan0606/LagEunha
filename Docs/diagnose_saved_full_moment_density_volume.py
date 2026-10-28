"""Density drift in cells wholly inside the uniform-density cavity.

Reuse saved IC quadrature errors and saved geometric Vdot; no new fluid RHS.
Separate IC error, exact fixed-mass volume change, and the initial tangent.
This identifies a kinematic contribution, not a unique dynamical cause.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from compare_full_moment_checkpoint_errors import load


def sha(file):return hashlib.sha256(Path(file).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    base=Path('/gpfs/kjhan/LagForce');ic=base/'Codex_full_curved_kepler_IC_20261019/initial.npz'
    mesh_file=base/'Codex_full_curved_disk_mesh_20261019/mesh.npz'
    native_file=base/'Codex_circular_sparse_support_full_20261019/full_exact_supports.json'
    tangent_root=base/'Codex_full_disk_saved_pressure_tangent_20261019'
    metadata=json.loads((tangent_root/'report.json').read_text())['stages']['initial']
    assert metadata['state_sha256']==sha(ic) and Path(metadata['geometry'])==mesh_file
    with np.load(ic) as f:initial={key:f[key] for key in f.files}
    with np.load(mesh_file) as f:mesh={key:f[key] for key in f.files}
    with np.load(tangent_root/'initial.npz') as f:tangent={key:f[key] for key in f.files}
    errors_file=ic.parent/'projection_errors.npz'
    with np.load(errors_file) as f:errors={key:f[key] for key in f.files}
    ids=initial['ids'];mass=initial['mass'];area=mesh['moments'][:,0]
    assert all(np.array_equal(ids,d['ids']) for d in (mesh,tangent,errors))
    assert float(initial['time'])==0 and not initial['frozen'].any()
    floor=json.loads((ic.parent/'report.json').read_text())['physics']['floor']
    sys.path.insert(0,str(base/'Codex_full_disk_arc_roots_modules_v7_20261019'))
    from curved_disk_cell import DiskCell
    native=json.loads(native_file.read_text());cavity=np.zeros(len(ids),bool)
    # W=0 Voronoi cells contain their sites. Exact endpoint radius tests then
    # establish that the whole initial convex cell lies strictly inside r=.5.
    assert not mesh['weights'].any()
    for i in np.flatnonzero(np.linalg.norm(mesh['sites'],axis=1)<.5000001):
        cell=DiskCell.from_native(native['cells'][i],2.5)
        cavity[i]=not cell.arcs and all(sum(x*x for x in point)<Fraction(1,4)
                                      for line in cell.lines for point in (line.A,line.B))
    assert cavity.any()
    initial_error=mass/(floor*area)-1
    density_rate=-tangent['Vdot'].astype(np.longdouble)/area.astype(np.longdouble)
    quadrature_mass_error=errors['raw_q24_minus_q48'][:,0]/mass
    def metrics(value,mask):
        value=np.asarray(value,dtype=np.longdouble)[mask];w=mass[mask].astype(np.longdouble)
        return dict(RMS=float(np.sqrt(np.sum(w*value*value)/w.sum())),maximum=float(abs(value).max()),signed_mean=float(np.sum(w*value)/w.sum()))
    root=base/'Codex_full_disk_continuous_two_dt_to_t1_20261019';checkpoints=[]
    args.output.mkdir(exist_ok=False)
    np.savez(args.output/'initial_density_volume.npz',ids=ids,initial_uniform_cavity=cavity,
        initial_area=area,initial_density_relative_error=initial_error,density_relative_rate=density_rate)
    for branch,steps in [('coarse',(2,5)),('fine',(4,10))]:
        for step in steps:
            folder=root/branch/'checkpoints'/f'{branch}_{step:04d}'/'physics'
            data,report,provenance=load(folder)
            assert np.array_equal(data['ids'],ids) and np.array_equal(data['actual_mass'],mass)
            mask=cavity&data['cavity'];assert mask.any()
            observed=data['relative_density_reference_error']
            geometric=mass/(floor*data['area'])-1
            predictor=initial_error+np.longdouble(report['time'])*density_rate
            difference=observed-predictor
            w=mass[mask].astype(np.longdouble)
            a=observed[mask].astype(np.longdouble);b=predictor[mask].astype(np.longdouble)
            norm=np.sqrt(np.sum(w*a*a)*np.sum(w*b*b))
            alignment=float(np.sum(w*a*b)/norm) if norm>0 else None
            checkpoints.append(dict(branch=branch,time=report['time'],cells=int(mask.sum()),source=provenance,
                observed_density_error=metrics(observed,mask),fixed_mass_volume_identity_error=metrics(observed-geometric,mask),
                current_reference_mass_vs_floor_area_error=metrics(data['reference_mass']/(floor*data['area'])-1,mask),
                initial_tangent_prediction=metrics(predictor,mask),observed_minus_initial_tangent=metrics(difference,mask),
                initial_tangent_alignment=alignment))
    result=dict(status='completed_saved_density_volume_drift_diagnostic',initial_uniform_cavity_cells=int(cavity.sum()),
        IC_density_error_in_uniform_cavity=metrics(initial_error,cavity),
        IC_q24_vs_q48_relative_mass_error=metrics(quadrature_mass_error,np.ones(len(ids),bool)),
        initial_density_rate_in_uniform_cavity=metrics(density_rate,cavity),checkpoints=checkpoints,
        source_sha256={str(file):sha(file) for file in (ic,mesh_file,native_file,tangent_root/'initial.npz',errors_file)},
        physical_providers=0,accepted_time_steps=0,
        limits='Initial cavity classified by exact rational vertices. Density-volume identity is kinematic; tangent is a short-time predictor, not a certified error bound or unique cause. Saved states/weights/masses unchanged.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],initial_cavity_error=result['IC_density_error_in_uniform_cavity'],checkpoints=[dict(branch=r['branch'],time=r['time'],observed=r['observed_density_error'],identity_error=r['fixed_mass_volume_identity_error'],alignment=r['initial_tangent_alignment']) for r in checkpoints])),flush=True)


if __name__=='__main__':main()

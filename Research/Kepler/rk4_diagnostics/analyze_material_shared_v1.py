"""Read native shared-kinematics MPI checks; no evolution."""
import json,hashlib,struct
from pathlib import Path

root=Path(__file__).resolve().parent
out=dict(scope='shared H1 material reference and source-locality diagnostics, not an evolved ALE method',
    previous_goal_turn='progress: local compatible RHS, finite affine remap and CUDA reference checks',
    current_goal_turn='progress: globally shared native MPI kinematics exposes local physical Euler ledger mismatch; source controls localize it',
    accepted_evolution=False,cases={},MPI_comparison={},native_handles=[])
for version in ['v1','v2']:
    d=root/('material_shared_'+version)
    assert (d/'status.txt').read_text().strip()=='AUDIT_COMPLETE'
    one=json.loads((d/'mpi1/contract.json').read_text())['cases']
    four=json.loads((d/'mpi4/contract.json').read_text())['cases']
    for a,b in zip(one,four):
        name=a['case'];assert name==b['case'];out['cases'][name]=a
        delta={k:abs(a[k]-b[k]) for k in ['global_energy_residual_over_U','global_angular_residual_over_2K',
                'max_cell_energy_residual_over_U','max_cell_angular_residual_over_2K']}
        def state(path):
            raw=path.read_bytes();magic,n=struct.unpack_from('<QQ',raw)
            assert magic==0x47465348414c3031
            return n,struct.unpack_from('<'+'d'*(6*n),raw,16)
        n,x=state(d/f'mpi1/{name}_state.bin');nn,y=state(d/f'mpi4/{name}_state.bin');assert n==nn
        delta['maximum_raw_coefficient_difference']=max(abs(v-w) for v,w in zip(x,y))
        delta['maximum_nodal_acceleration_difference']=max(abs(x[k]-y[k]) for block in [4,5] for k in range(block*n,block*n+a['nodes']))
        out['MPI_comparison'][name]=delta
out['interpretation']=[
    'Shared position/velocity/acceleration traces and global conservation alone do not establish the original per-cell Euler flux/work ledger.',
    'The local energy residual is not Pdot or heat injection: physical e is still paired with pressure work.',
    'Gravity-off and representable linear-source controls reduce the local defect; actual Kepler source projection yields ~16 U/time.',
    'The exact factorization v^T(M_cell*a-f_cell) and its angular counterpart identify variational intercell mass/load coupling.',
    'Whether to derive legitimate hybrid interface tractions or revise the kinematic/momentum structure is a bundled Astra evaluation dependency.',
    'Curved Kepler remap, shock stress, complete update stability and t6-10 remain unverified.']
out['sources']={}
for version in ['v1','v2']:
    d=root/('material_shared_'+version)
    out['sources'][version]=dict(directory=str(d),executable_sha256=hashlib.sha256((d/'contract.exe').read_bytes()).hexdigest())
p=root/'MATERIAL_SHARED_BUNDLE_20261021.json'
with p.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(dict(summary=str(p),MPI_comparison=out['MPI_comparison']),indent=2))

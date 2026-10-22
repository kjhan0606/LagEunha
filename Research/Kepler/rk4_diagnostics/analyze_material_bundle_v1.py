"""Read frozen native material-map/remap/CUDA checks; no dynamics."""
import hashlib,json
from pathlib import Path

root=Path(__file__).resolve().parent
def packet(name):
    d=root/name
    assert (d/'status.txt').read_text().strip()=='AUDIT_COMPLETE'
    return json.loads((d/'contract.json').read_text())

coarse=packet('material_rhs_v1');fine=packet('material_rhs_v3');gpu=packet('material_cuda_v1')
remap_dir=root/'material_remap_v4'
assert (remap_dir/'status.txt').read_text().strip()=='AUDIT_COMPLETE'
raw=(remap_dir/'contract.json').read_text()
# Preserve the native output verbatim. Its one JSON-format typo .125 is
# normalized in this derived record only; no value or physics result changes.
needle='"finite_shear_strain":.125'
assert raw.count(needle)==1
remap=json.loads(raw.replace(needle,'"finite_shear_strain":0.125'))
(remap_dir/'analysis.json').write_text(json.dumps(dict(
    source='contract.json',source_sha256=hashlib.sha256(raw.encode()).hexdigest(),
    serialization_normalization='one numeric literal .125 -> 0.125; raw native artifact retained',
    native_packet=remap),indent=2)+'\n')
indexed={(x['id'],x['case']):x for x in coarse['cells']}
quantities=['initial_velocity_error','coordinate_error','K_minus_reference_over_U',
  'energy_work_residual','energy_work_relative','angular_residual_over_2K',
  'max_Kepler_divergence_error','max_abs_Pdot_over_P0',
  'volume_vs_boundary_rate_over_area','finite_map_volume_rate_error_over_area',
  'native_initial_mass_relative_error','native_initial_energy_error_over_U',
  'mass_matrix_solve_relative','mass_solve_kinetic_power_residual']
cases={}
for name in ['Kepler','Kepler_perturb']:
    rows=[x for x in fine['cells'] if x['case']==name]
    cases[name]=dict(cells=len(rows),maximum_absolute={k:max(abs(x[k]) for x in rows) for k in quantities})
report=dict(scope='partial first material-ALE bundle: cell-local smooth material RHS, affine conservative remap controls, native CUDA helper; no evolved Kepler state',
    previous_goal_turn='progress: rejected nodal operators by measured thermal, angular and Jacobian evidence',
    current_goal_turn='progress: kinematic obstruction, compatible native map/stress/moments, finite mass-weight remap controls and CUDA matching',
    final_goal='t6-10 and longest GIZMO Kepler stability; not achieved',
    cases=cases,
    quadrature=dict(coarse='32x20',fine='64x28',
      max_mass_relative_change=max(abs(x['mass']-indexed[x['id'],x['case']]['mass'])/x['mass'] for x in fine['cells']),
      max_energy_change_over_thermal=max(abs(x['K']+x['U']-indexed[x['id'],x['case']]['K']-indexed[x['id'],x['case']]['U'])/x['U'] for x in fine['cells'])),
    intrinsic_rotation=dict(scope='actual represented initial field and material RHS; no separately advanced spin reservoir or relaxation',
      minimum_resolved_strain_remainder=min(x['resolved_strain_remainder'] for x in fine['cells'])),
    remap=remap,cuda=gpu,
    accepted_evolution=False,
    incomplete=['Global continuous/shared kinematic space and mass assembly; current RHS is cell-local at matching initial traces.',
      'Curved material faces and remap of the actual nonuniform-density Kepler field.',
      'Cold-conservative finite representation/fragment compression for repeated remaps.',
      'Fixed-target density-wave composed-update response and historical cs/h sliding mode.',
      'Multidimensional symmetric Riemann compression stress, physical shock heating and positivity.',
      'Actual native CUDA+MPI runtime/halo/migration/checkpoint/time stepping.',
      'Historical t6-10 range and longest matched-GIZMO duration.'],
    native_handles=[],artifacts={})
for name in ['affine_face_kinematics_v1','material_rhs_v1','material_rhs_v2','material_rhs_v3','material_remap_v3','material_remap_v4','material_cuda_v1']:
    d=root/name
    assert (d/'status.txt').read_text().strip()=='AUDIT_COMPLETE'
    report['artifacts'][name]=dict(path=str(d),executable_sha256=hashlib.sha256((d/'contract.exe').read_bytes()).hexdigest(),
        raw_contract_sha256=hashlib.sha256((d/'contract.json').read_bytes()).hexdigest())
dest=root/'MATERIAL_ALE_BUNDLE_20261021.json'
with dest.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(dict(summary=str(dest),cases=cases,quadrature=report['quadrature'],cuda=gpu),indent=2))

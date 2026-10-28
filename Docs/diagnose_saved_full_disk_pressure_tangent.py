"""Pressure tangent from saved geometry/rates, without a fluid RHS or evolution.

The exact moving power-plane speed gives Vdot. The saved raw moment RHS and
the infinitesimal exponential/QR covariance map give two thermal derivatives.
These are diagnostics, not equilibrium corrections or a new acceptance gate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, required=True)
args = p.parse_args()
base = Path('/gpfs/kjhan/LagForce')
ld = np.longdouble
icpath = base/'Codex_full_curved_kepler_IC_20261019/initial.npz'
ic = np.load(icpath)
ids, mass = ic['ids'], ic['mass']
n = len(ids)
physics = json.loads((icpath.parent/'report.json').read_text())['physics']
gm, P0 = ld(physics['gamma'])-1, ld(physics['P0'])
r0 = np.linalg.norm(ic['state'][:-4].reshape(n, 20)[:, :2], axis=1)
regions = dict(all=np.ones(n, bool), inner=r0 < .4,
               disk_initial_radius_0p5_to_2=(r0 > .5) & (r0 < 2))
lookup = {int(g): i for i, g in enumerate(ids)}
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()

def stats(v, mask):
    v, w = v[mask], mass[mask].astype(ld)
    return dict(RMS=float(np.sqrt(np.average(v*v, weights=w))),
                maximum_absolute=float(abs(v).max()),
                signed_mass_mean=float(np.average(v, weights=w)))

rows = [
 ('initial', icpath, 'Codex_full_curved_disk_mesh_20261019',
  'Codex_full_disk_initial_stage0_arc_roots_retry_20261019'),
 ('coarse_midpoint', base/'Codex_full_disk_midpoint_predictors_20261019/coarse_midpoint/state.npz',
  'Codex_full_disk_coarse_midpoint_geometry_20261019', 'Codex_full_disk_coarse_midpoint_stage_20261019'),
 ('fine_first_midpoint', base/'Codex_full_disk_midpoint_predictors_20261019/fine_first_midpoint/state.npz',
  'Codex_full_disk_fine_midpoint_geometry_20261019', 'Codex_full_disk_fine_midpoint_stage_20261019'),
 ('fine_second_stage0', base/'Codex_full_disk_fine_first_accepted_20261019/accepted.npz',
  'Codex_full_disk_fine_first_endpoint_geometry_20261019', 'Codex_full_disk_fine_second_stage0_20261019'),
 ('fine_second_midpoint', base/'Codex_full_disk_fine_second_midpoint_predictor_20261019/candidate/state.npz',
  'Codex_full_disk_fine_second_midpoint_geometry_20261019', 'Codex_full_disk_fine_second_midpoint_stage_20261019')]
args.output.mkdir(exist_ok=False)
report = {}
for label, statepath, geomname, stagename in rows:
    data = np.load(statepath)
    meshpath = base/geomname/'mesh.npz'
    mesh = np.load(meshpath)
    stagepath = base/stagename/'stage0.npz'
    saved = np.load(stagepath)
    a = data['state'][:-4].reshape(n, 20)
    assert np.array_equal(data['ids'], ids) and np.array_equal(data['mass'], mass)
    assert np.array_equal(mesh['ids'], ids) and np.array_equal(mesh['sites'], a[:, :2])
    assert np.array_equal(mesh['weights'], a[:, 2])
    assert saved['initial_state_sha256'].item() == hashlib.sha256(np.ascontiguousarray(data['state']).view(np.uint8)).hexdigest()
    st = saved['stage'].astype(ld)
    x, gen, z = a[:, :2].astype(ld), saved['generator'].astype(ld), saved['weight_rates'].astype(ld)
    ij = np.array([(lookup[int(i)], lookup[int(j)]) for i, j in mesh['face_keys']])
    i, j = ij.T
    face = mesh['faces'].astype(ld)
    mid, length, normal = face[:, :2], face[:, 4], face[:, 5:7]
    delta, dgen = x[j]-x[i], gen[j]-gen[i]
    distance = np.sqrt(np.sum(delta*delta, axis=1))
    normal_error = np.max(abs(delta/distance[:, None]-normal))
    assert normal_error < 1e-12, float(normal_error)
    # d/dt[(xj-xi).x=(|xj|²-|xi|²-wj+wi)/2].
    # This centered expression avoids cancellation of large position products.
    wn = (.5*np.sum((x[j]+x[i]-2*mid)*dgen, axis=1)
          + .5*np.sum(delta*(gen[j]+gen[i]), axis=1)-.5*(z[j]-z[i]))/distance
    face_rate = length*wn
    Vdot = np.zeros(n, dtype=ld)
    np.add.at(Vdot, i, face_rate)
    np.add.at(Vdot, j, -face_rate)
    # Outer circular wall is stationary; no arc contribution.
    volume_sum = Vdot.sum()
    volume_bound = 64*np.finfo(ld).eps*abs(face_rate).sum()
    assert abs(volume_sum) <= volume_bound
    raw = np.column_stack([mass, a[:, 3:]]).astype(ld)
    m, q, mom = raw[:, 0], raw[:, 1:3], raw[:, 7:9]
    outer = lambda u,v: u[:, :, None]*v[:, None, :]
    S = raw[:, 3:7].reshape(n, 2, 2)-outer(q,q)/m[:, None, None]
    C = raw[:, 9:13].reshape(n, 2, 2)-outer(q,mom)/m[:, None, None]
    B = np.linalg.solve(S.astype(float), C.astype(float)).transpose(0, 2, 1).astype(ld)
    # Stored mass is fixed by the integrator, not advanced by tiny RHS residuals.
    qdot, pdot = st[:, 1:3], st[:, 7:9]
    Sdot = st[:, 3:7].reshape(n, 2, 2)-(outer(qdot,q)+outer(q,qdot))/m[:, None, None]
    Cdot = st[:, 9:13].reshape(n, 2, 2)-(outer(qdot,mom)+outer(q,pdot))/m[:, None, None]
    L = data['factors']
    u,s,vh = np.linalg.svd(L)
    scaled = s/s.max(axis=1)[:, None]
    ratio = scaled[:, :, None]*scaled[:, None, :]/(scaled[:, :, None]**2+scaled[:, None, :]**2)
    K = np.asarray(st[:, 22:26].reshape(n, 2, 2), float)
    loss = (u@((vh@K@vh.transpose(0,2,1))*ratio)@u.transpose(0,2,1)).astype(ld)
    R = L.astype(ld)@L.astype(ld).transpose(0,2,1)
    Ac = st[:, 56:60].reshape(n,2,2)-loss
    Rdot = Ac@R+R@Ac.transpose(0,2,1)+st[:, 26:30].reshape(n,2,2)
    bc = B@Cdot
    Pi_factor_dot = ((outer(pdot,mom)+outer(mom,pdot))/m[:,None,None]
                     +bc+bc.transpose(0,2,1)-B@Sdot@B.transpose(0,2,1)+Rdot)
    trace = lambda t: np.trace(t, axis1=1, axis2=2)
    Pi_raw_dot = st[:,13:17].reshape(n,2,2)
    V = mesh['moments'][:,0].astype(ld)
    thermal = raw[:,17]-.5*(raw[:,13]+raw[:,16])
    scale = gm/(P0*V)
    fields = dict(
        volume_contribution=-scale*thermal*Vdot/V,
        raw_thermal_transport=scale*(st[:,17]-.5*trace(Pi_raw_dot)),
        factor_minus_raw_thermal=-.5*scale*trace(Pi_factor_dot-Pi_raw_dot),
        direct_gravity_thermal_subset=scale*(st[:,47]-.5*(st[:,43]+st[:,46])))
    fields['raw_relative_pressure_rate'] = fields['volume_contribution']+fields['raw_thermal_transport']
    fields['factor_relative_pressure_rate'] = fields['raw_relative_pressure_rate']+fields['factor_minus_raw_thermal']
    assert all(np.isfinite(v).all() for v in fields.values())
    np.savez(args.output/(label+'.npz'), ids=ids, Vdot=Vdot, **fields)
    report[label] = dict(time=float(data['time']), state=str(statepath), state_sha256=sha(statepath),
        stage=str(stagepath), geometry=str(meshpath), maximum_normal_error=float(normal_error),
        total_volume_rate=float(volume_sum), volume_rate_roundoff_bound=float(volume_bound),
        ignored_fixed_mass_RHS_residual=float(abs(st[:,0]).max()),
        regions={name:{k:stats(v,mask) for k,v in fields.items()} for name,mask in regions.items()})
result = dict(status='completed_saved_pressure_tangent_diagnostic', stages=report,
    provider_calls=0, evolution_steps=0, runner_sha256=sha(Path(__file__)),
    limits='Analytic tangent of saved numerical maps. Batched binary64 SVD/solve with longdouble assembly; no certified tangent error envelope. No physical correction, state reset, or unique-cause conclusion. Outer wall stationary; fixed-mass derivative uses mdot=0 exactly. Component RMS values can cancel.')
(args.output/'report.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(status=result['status'], regions={k:v['regions']['inner'] for k,v in report.items()}),indent=2))

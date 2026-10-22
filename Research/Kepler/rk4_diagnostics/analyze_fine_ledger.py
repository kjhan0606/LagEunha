"""Read-only native failure ledger and independent Voronoi Jacobian checks.

Run in a Slurm allocation. No fluid evolution, orbit forcing, or state repair.
"""
from pathlib import Path
import csv, ctypes, hashlib, json, sys, time
import numpy as np
from scipy.spatial import Voronoi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from audit_snapshot import read_snapshot

BASE = ROOT / 'rk4_diagnostics/fine_cold_interval'
OUT = BASE / 'ledger'
OUT.mkdir(exist_ok=False)
manifest = json.loads((BASE / 'manifest.json').read_text())
p, meta = read_snapshot(manifest['source'])
ic, _ = read_snapshot(ROOT / 'core_refinement/fine_integrated/khout.000001.dat')
ids = (p['u4if'] & np.uint64((1 << 56) - 1)).astype(int)
ic_ids = (ic['u4if'] & np.uint64((1 << 56) - 1)).astype(int)
lookup = {int(v): k for k, v in enumerate(ids)}
iclookup = {int(v): k for k, v in enumerate(ic_ids)}
N = len(p)
gamma = float(np.float32(1.666667))
xy = np.column_stack((p['x'], p['y']))
vel = np.column_stack((p['vx'], p['vy']))
box = np.array([meta['Lx'], meta['Ly']])
minimum_image = lambda a: a - box * np.rint(a / box)

with (BASE / 'full/rk4_stage_r0.tsv').open() as fh:
    rows = list(csv.DictReader(fh, delimiter='\t'))
assert len(rows) == N and all(r['stage'] == 'k1_rhs' for r in rows)
ordered = [None] * N
for row in rows:
    ordered[lookup[int(row['id'])]] = row
arr = lambda key: np.array([float(r[key]) for r in ordered])
P, V, U, mass = arr('P_stored'), arr('V'), arr('U_raw'), arr('mass')
rho = arr('rho_stored')
native_F = mass[:, None] * np.column_stack((arr('ax_hydro'), arr('ay_hydro')))
native_Edot, native_Udot = arr('die_stored'), arr('Udot')
assert np.max(abs(arr('x') - xy[:, 0])) < 1e-13
assert np.max(abs(arr('y') - xy[:, 1])) < 1e-13
assert np.max(abs(arr('vx') - vel[:, 0])) < 1e-13
assert np.max(abs(arr('vy') - vel[:, 1])) < 1e-13

f = np.loadtxt(BASE / 'full/face_dump_r0_call1.tsv', skiprows=1)
assert not np.any(f[:, 6])
i = np.array([lookup[int(v)] for v in f[:, 4]])
j = np.array([lookup[int(v)] for v in f[:, 5]])
S = f[:, 10:12]
corner1 = minimum_image(f[:, 12:14] - xy[i])
corner2 = minimum_image(f[:, 14:16] - xy[i])
d = minimum_image(xy[j] - xy[i])
d2 = np.sum(d*d, axis=1)
assert np.all(d2 > 0)
c = .5*(corner1 + corner2 - d)
dv = vel[j] - vel[i]
dot = lambda a, b: np.sum(a*b, axis=1)
D = dot(d, S) / d2
geom = .5*dot(vel[i] + vel[j], S) - dot(dv, c)*D
work_velocity = f[:, 8]
viS = dot(vel[i], S)
s_geom, s_used = geom - viS, work_velocity - viS
reduce = lambda x: np.bincount(i, weights=x, minlength=N)
closure = np.column_stack((reduce(S[:, 0]), reduce(S[:, 1])))
Vdot = reduce(geom)
recon_F = np.column_stack((reduce(-f[:, 7]*S[:, 0]), reduce(-f[:, 7]*S[:, 1])))
recon_Edot = reduce(f[:, 9])
recon_Udot = reduce(f[:, 9] + f[:, 7]*viS)
residual = recon_Udot + P*Vdot
pressure_residual = reduce((P[i] - f[:, 7])*s_geom)
velocity_residual = reduce(-f[:, 7]*(work_velocity - geom))
closure_residual = P*dot(vel, closure)

# Exact production pair-pressure helper. Pair input U and dtold are not changed.
dist = np.sqrt(d2)
q = dot(dv, d) / dist
cs = np.sqrt(gamma*P/rho).astype(np.float32).astype(float)
inputs = np.column_stack((dist, V[i], V[j], rho[i], rho[j], cs[i], cs[j],
                         np.maximum(-q, 0), np.linalg.norm(S, axis=1),
                         np.full(len(f), meta['dtold']), mass[i], mass[j], q,
                         U[i], U[j], s_used, dot(dv, S)-s_used, np.zeros(len(f))))
pair_output = np.empty((len(f), 2))
lib = ctypes.CDLL(str(ROOT / 'affine_limiter_intervention/pair_wrapper.so'))
lib.eval_pairs.argtypes = [ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')]
lib.eval_pairs(len(f), np.ascontiguousarray(inputs), pair_output)
pair = pair_output[:, 1]
pair_Udot = reduce(-pair*s_used)

# Fixed-weight true-volume Jacobian: J_ii(face)=S/2+c*(d.S)/d^2.
Jself = .5*S + c*D[:, None]
adj_F = np.column_stack((reduce((P[i]-P[j])*Jself[:, 0]),
                         reduce((P[i]-P[j])*Jself[:, 1])))
adj_Udot = -P*Vdot
adj_power = float(np.sum(dot(vel, adj_F)) + np.sum(adj_Udot))
adj_force_scale = float(np.sum(np.linalg.norm(adj_F, axis=1)))

focus_ids = sorted({int(r.split('id=')[1].split()[0])
                    for r in manifest['cases'][0]['invalid_stages']})
selected_ids = set(manifest['cases'][0]['selected_ids_for_subdivision_stages'])
for k in focus_ids:
    selected_ids.update(f[f[:, 4] == k, 5].astype(int).tolist())
selected = np.array([lookup[k] for k in sorted(selected_ids)])

def cell_row(k):
    a = iclookup[int(ids[k])]
    Uad = float(ic['ie'][a]) * (float(ic['volume'][a])/V[k])**(gamma-1)
    return dict(id=int(ids[k]), radius=float(np.linalg.norm(xy[k]-12)),
        initial_radius=float(np.hypot(ic['x'][a]-12, ic['y'][a]-12)),
        U=float(U[k]), U_adiabatic_from_IC=Uad, entropy_ratio=float(U[k]/Uad),
        V=float(V[k]), Vdot_geometric=float(Vdot[k]), P=float(P[k]),
        native_Udot=float(native_Udot[k]), reversible_Udot=float(adj_Udot[k]),
        residual_Udot_plus_PVdot=float(residual[k]),
        pressure_difference_contribution=float(pressure_residual[k]),
        work_velocity_contribution=float(velocity_residual[k]),
        closure_contribution=float(closure_residual[k]), pair_Udot=float(pair_Udot[k]),
        native_force=native_F[k].tolist(), adjoint_force=adj_F[k].tolist(),
        FE_U_positive_time=float(U[k]/-native_Udot[k]) if native_Udot[k] < 0 else None)

# Preserve each face contribution with signs; no division by almost-zero PdV.
face_keys = ['id_i', 'id_j', 'P_i', 'P_j', 'Pstar', 'Ppair', 'distance',
             'Sx', 'Sy', 'W_geom', 'W_used', 's_geom', 's_used',
             'force_x', 'force_y', 'Edot', 'Udot',
             'pressure_residual', 'work_velocity_residual']
with (OUT / 'selected_faces.tsv').open('w') as fh:
    w = csv.writer(fh, delimiter='\t'); w.writerow(face_keys)
    for k in np.flatnonzero(np.isin(ids[i], sorted(selected_ids))):
        w.writerow([int(ids[i[k]]), int(ids[j[k]]), P[i[k]], P[j[k]], f[k, 7], pair[k], dist[k],
            *S[k], geom[k], work_velocity[k], s_geom[k], s_used[k],
            *(-f[k, 7]*S[k]), f[k, 9], f[k, 9]+f[k, 7]*viS[k],
            (P[i[k]]-f[k, 7])*s_geom[k], -f[k, 7]*(work_velocity[k]-geom[k])])

# Independent Qhull recomputation under x +/- delta*v. Only interior cells
# satisfying edge_distance > 2*circumradius are compared, proving that omitted
# periodic images cannot affect these cells. The periodic lattice is fixed.
def independent_volume(points):
    vor = Voronoi(points)
    areas, radii = [], []
    for k in selected:
        reg = vor.regions[vor.point_region[k]]
        assert reg and -1 not in reg
        a = vor.vertices[reg] - points[k]
        b = np.roll(a, -1, axis=0)
        areas.append(.5*abs(np.sum(a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0])))
        radii.append(np.max(np.linalg.norm(a, axis=1)))
    edge = np.min(np.column_stack((points[selected], box-points[selected])), axis=1)
    assert np.all(edge > 2*np.array(radii))
    return np.array(areas)

start = time.monotonic()
independent_V = independent_volume(xy)
fd_rows = []
for delta in (1e-5, 5e-6, 2.5e-6, 1.25e-6):
    finite_diff = (independent_volume(xy+delta*vel)-independent_volume(xy-delta*vel))/(2*delta)
    err = finite_diff-Vdot[selected]
    fd_rows.append(dict(delta=delta, max_absolute_error=float(abs(err).max()),
        l2_relative_error=float(np.linalg.norm(err)/max(np.linalg.norm(Vdot[selected]), 1e-300)),
        focus=[dict(id=int(ids[k]), analytic=float(Vdot[k]), finite_difference=float(finite_diff[a]),
                    absolute_error=float(err[a])) for a, k in enumerate(selected) if int(ids[k]) in focus_ids]))

# Passing quarter path: retain stage values, RK combination, and actual entropy
# from IC (an adiabatic reference, not a forced solution under Riemann exchange).
with (BASE / 'quarters/rk4_stage_r0.tsv').open() as fh:
    quarter_rows = list(csv.DictReader(fh, delimiter='\t'))
quarter = {}
for ident in focus_ids:
    records = [r for r in quarter_rows if int(r['id']) == ident]
    states = []
    for r in records:
        if r['stage'] not in ('k1_rhs', 'k2_rhs', 'k3_rhs', 'k4_rhs', 'final_geometry'):
            continue
        a = iclookup[ident]; vol, internal = float(r['V']), float(r['U_raw'])
        Uad = float(ic['ie'][a])*(float(ic['volume'][a])/vol)**(gamma-1)
        states.append(dict(step=int(r['step']), stage=r['stage'], t_start=float(r['t_start']),
            h=float(r['h']), U=internal, V=vol, U_adiabatic_from_IC=Uad,
            entropy_ratio=internal/Uad, Udot=float(r['Udot']) if r['rhs_valid']=='1' else None))
    quarter[str(ident)] = states

result = dict(scope='same captured failure interval; no long-run stability claim or adopted new operator',
    source=manifest['source'], replay=manifest, native_first_RHS_time=meta['t'],
    gamma=gamma, all_cells=N, faces=len(f), selected_ids=sorted(selected_ids), focus_ids=focus_ids,
    reconstruction=dict(force_max_absolute_error=float(abs(recon_F-native_F).max()),
        Edot_max_absolute_error=float(abs(recon_Edot-native_Edot).max()),
        Udot_max_absolute_error=float(abs(recon_Udot-native_Udot).max()),
        residual_decomposition_max_absolute_error=float(abs(residual-pressure_residual-velocity_residual-closure_residual).max()),
        maximum_surface_closure=float(np.linalg.norm(closure, axis=1).max()),
        selected_independent_volume_relative_error_upper_bound=float(abs(independent_V-V[selected]).max()/independent_V.min()),
        pressure_EOS_relative_error_max=float(np.max(abs(P-(gamma-1)*U/V)/P))),
    adjoint_reference=dict(description='F_i=sum_j P_j dV_j/dx_i; Udot_i=-P_i Vdot_i, fixed weights; no evolution performed',
        net_force=adj_F.sum(axis=0).tolist(), sum_absolute_force=adj_force_scale,
        total_power_residual=adj_power,
        sum_absolute_pressure_work=float(np.sum(abs(P*Vdot))),
        native_net_force=native_F.sum(axis=0).tolist(),
        native_total_hydro_power=float(np.sum(native_Edot)),
        finite_difference=fd_rows, finite_difference_seconds=time.monotonic()-start,
        torque_caveat='Fixed periodic square is not globally rotation invariant. No unjustified zero global torque assertion.'),
    selected_cells=[cell_row(k) for k in selected], quarter_focus_stages=quarter,
    pair_caveat='Production helper is used. cs is reconstructed from stored float P/rho; exact pair attribution is asserted only where base pair pressure is zero or reproduction precision is validated.',
    constraints='No steering, no physical viscosity, no mass/weight change, no floor, no spin damping.')
with (OUT / 'results.json').open('w') as fh:
    json.dump(result, fh, indent=2, allow_nan=False); fh.write('\n')
print(json.dumps({k: result[k] for k in ('reconstruction', 'adjoint_reference')}, indent=2), flush=True)
print('focus', json.dumps([cell_row(lookup[k]) for k in focus_ids], indent=2), flush=True)

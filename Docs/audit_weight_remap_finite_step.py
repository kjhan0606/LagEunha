"""Frozen-face RK4 diagnostic using the experimental C/CUDA remap kernel.

Run on lageunha. This is not a Kepler evolution: positions and face rates are
frozen, only the conservative cell state evolves. No gravity force is applied.
"""
import argparse
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output', type=Path)
    args = ap.parse_args()
    args.output.mkdir(exist_ok=False)
    base = Path('/gpfs/kjhan/LagForce')
    kernel = base / 'Codex_weight_remap_flux_contact_20261017'
    geom = base / 'Codex_weighted_compression_geometry_20261017/audit.py'
    source = base / 'Codex_gradient_operator_20261017/operator/rhs.dat'
    spec = importlib.util.spec_from_file_location('geometry', geom)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    a = np.loadtxt(source)
    a = a[np.argsort(a[:, 0])]
    x, v = a[:, 1:3]-3, a[:, 3:5]
    r = np.linalg.norm(x, axis=1)
    omega = (x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])/r**2
    omega[r >= 2.5] = 0
    c, s = np.cos(omega*.001), np.sin(omega*.001)
    x = np.column_stack((c*x[:, 0]-s*x[:, 1], s*x[:, 0]+c*x[:, 1]))
    v = np.column_stack((-omega*x[:, 1], omega*x[:, 0]))
    target = base / 'Codex_evolving_weight_transport_20261017'
    left, right = np.load(target/'weights_0.npz'), np.load(target/'weights_0.001.npz')
    ids, w = right['ids'], right['weights']
    wd = (w-left['weights'])/.001
    volume, _, lap, _ = mod.geometry(x, w, v, ids)
    lap = lap.tocoo()
    ii, jj = ids[lap.row], lap.col
    take = (ii < jj) & (wd[ii] != wd[jj])
    ii, jj = ii[take], jj[take]
    rate = .01*(-lap.data[take])*(wd[ii]-wd[jj])
    nodes = np.unique(np.r_[ii, jj])
    assert np.isin(nodes, ids).all(), 'All exchanging cells must have measured volumes'
    index = np.full(len(a), -1, dtype=int)
    index[nodes] = np.arange(len(nodes))
    i, j = index[ii], index[jj]
    volumes = dict(zip(ids, volume))
    mass = .01*np.array([volumes[k] for k in nodes])
    pos, vel = x[nodes], v[nodes]
    gamma = 1.6666669845581055
    e = 2.7777796321446218e-5/((gamma-1)*.01)
    initial = np.column_stack((mass, mass[:, None]*vel,
                              mass*(e+.5*np.sum(vel**2, axis=1))))
    cpu = ctypes.CDLL(str(kernel/'cpu.so'))
    gpu = ctypes.CDLL(str(kernel/'cuda.so'))
    arr = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
    for fun in (cpu.wr_flux_batch, gpu.wr_cuda_batch):
        fun.argtypes = [ctypes.c_int, ctypes.c_size_t, arr, arr]
        fun.restype = ctypes.c_int
    cpu.wr_set_threads(4)

    def rhs(q, fun):
        m = q[:, 0]
        assert np.all(m > 0)
        u = q[:, 1:3]/m[:, None]
        thermal = q[:, 3]/m-.5*np.sum(u*u, axis=1)
        data = np.ascontiguousarray(np.column_stack((rate, u[i], u[j], thermal[i], thermal[j])))
        flux = np.empty((len(i), 4))
        assert fun(2, len(i), data, flux) == 0
        out = np.zeros_like(q)
        np.add.at(out, i, flux)
        np.add.at(out, j, -flux)
        return out

    def energy(q):
        k = .5*np.sum(q[:, 1:3]**2, axis=1)/q[:, 0]
        return k, q[:, 3]-k

    def angular(q):
        return float(np.sum(pos[:, 0]*q[:, 2]-pos[:, 1]*q[:, 1]))

    # Match the canonical spline potential, RC=.5, about (3,3).
    u = np.linalg.norm(pos, axis=1)/.5
    phi = np.empty(len(u))
    low, mid, high = u < .5, (u >= .5) & (u < 1), u >= 1
    z = u[low]
    phi[low] = (16/3*z*z-48/5*z**4+32/5*z**5-14/5)/.5
    z = u[mid]
    phi[mid] = (1/(15*z)+32/3*z*z-16*z**3+48/5*z**4-32/15*z**5-16/5)/.5
    phi[high] = -1/(.5*u[high])
    derivative = rhs(initial, cpu.wr_flux_batch)
    mdot = derivative[:, 0]
    depletion = float(np.min(mass[mdot < 0]/-mdot[mdot < 0]))
    horizon = min(.001, .2*depletion)
    k0, u0 = energy(initial)
    rows = []
    for steps in (1, 2, 4, 8, 16):
        finals = []
        for name, fun in (('openmp4', cpu.wr_flux_batch), ('cuda', gpu.wr_cuda_batch)):
            q = initial.copy()
            dt = horizon/steps
            minimum_stage_u = float(u0.min())
            for _ in range(steps):
                k1 = rhs(q, fun)
                q2 = q+.5*dt*k1
                k2 = rhs(q2, fun)
                q3 = q+.5*dt*k2
                k3 = rhs(q3, fun)
                q4 = q+dt*k3
                k4 = rhs(q4, fun)
                q = q+dt/6*(k1+2*k2+2*k3+k4)
                minimum_stage_u = min(minimum_stage_u, *(float(energy(z)[1].min()) for z in (q2, q3, q4, q)))
            assert np.isfinite(q).all()
            k, thermal = energy(q)
            change = q-initial
            rows.append({'backend': name, 'steps': steps, 'dt': dt,
                         'global_m_p_E_change': change.sum(axis=0).tolist(),
                         'kinetic_change_over_initial_thermal': float((k-k0).sum()/u0.sum()),
                         'thermal_change_over_initial_thermal': float((thermal-u0).sum()/u0.sum()),
                         'minimum_stage_thermal_energy': minimum_stage_u,
                         'minimum_mass_over_initial': float(np.min(q[:, 0]/mass)),
                         'point_angular_change': angular(q)-angular(initial),
                         'point_potential_energy_change': float(phi@change[:, 0])})
            finals.append(q)
        assert np.array_equal(finals[0].view('u8'), finals[1].view('u8'))
    # A boost changes point torque by the first mass-moment transfer x mdot.
    boost = np.array([17., -9.])
    boosted = initial.copy()
    boosted[:, 1:3] += mass[:, None]*boost
    boosted[:, 3] += initial[:, 1:3]@boost+.5*(boost@boost)*mass
    boosted_rhs = rhs(boosted, cpu.wr_flux_batch)
    first_moment_rate = pos.T@mdot
    expected = first_moment_rate[0]*boost[1]-first_moment_rate[1]*boost[0]
    residual = angular(boosted_rhs)-angular(derivative)-expected
    assert abs(residual) < 1e-13
    result = {'scope': 'Frozen actual face graph and constant weight-induced mass rates; uniform-density manufactured state at canonical positions. No generator evolution, face rebuilding, hydro forces or MPI. Not a Kepler stability run.',
              'kernel_mode': 2, 'cells': len(nodes), 'faces': len(i), 'horizon': horizon,
              'mass_depletion_time': depletion, 'all_cpu_cuda_final_states_bitwise_equal': True,
              'rows': rows, 'point_mass_first_moment_rate': first_moment_rate.tolist(),
              'boosted_point_torque_identity_residual': float(residual),
              'point_potential_rate': float(phi@mdot),
              'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, geom, kernel/'weight_remap_flux.h', kernel/'cpu.so', kernel/'cuda.so', Path(__file__))},
              'interpretation': 'Time-step convergence can reduce kinetic integration error, but cannot repair the semidiscrete point torque or point potential transfer. Full cell moments and compatible potential transport remain required; adding spin without kinetic realizability is insufficient.'}
    (args.output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

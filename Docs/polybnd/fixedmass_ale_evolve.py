"""Four-cell compressible ALE prototype with fixed cell masses, no gravity.

Closure: piecewise-constant Euler states, so generator velocity equals p_i/m_i.
No centroid correction. Weight rates enforce zero NET mass flux per cell;
momentum/total-energy fluxes remain. Shared faces use one Rusanov ALE flux.
Rigid stationary walls; their impulse is recorded, energy flux is zero.
Research prototype, not a Kepler test or the enriched-subcell method.
"""
import json
import numpy as np
from shapely.geometry.polygon import orient
from weighted_mass_patch import power_cells

GAMMA = 5/3
DOMAIN = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])


def unpack(y):
    return y[:8].reshape(4, 2), y[8:12], y[12:20].reshape(4, 2), y[20:24]


def faces(sites, weights, cells):
    result = []
    for i, c in enumerate(cells):
        p = np.array(orient(c, sign=1).exterior.coords[:-1])
        for a, b in zip(p, np.roll(p, -1, axis=0)):
            edge = b-a; length = np.linalg.norm(edge)
            if length == 0: continue
            normal = np.array([edge[1], -edge[0]])/length; mid = (a+b)/2
            neighbor = -1
            for j in range(4):
                if j == i: continue
                aa = 2*(sites[j]-sites[i]); bb = np.dot(sites[j], sites[j])-np.dot(sites[i], sites[i])+weights[i]-weights[j]
                if abs(np.dot(aa, mid)-bb) < 1e-10:
                    neighbor = j; break
            if neighbor < 0 or i < neighbor: result.append((i, neighbor, mid, normal, length))
    return result


def initial(kind):
    sites = np.array([[-.55, -.48], [.48, -.55], [.53, .49], [-.46, .53]])
    weights = np.array([.03, -.01, .02, 0.])
    area = np.array([p.area for p in power_cells(sites, weights, DOMAIN)])
    rho = np.ones(4); pressure = np.ones(4); velocity = np.zeros((4, 2))
    if kind == 'acoustic': pressure = np.array([1.05, .95, 1.02, .98])
    elif kind == 'contact': rho = np.array([1., .7, .7, 1.])
    elif kind != 'uniform': raise ValueError(kind)
    mass = rho*area; momentum = mass[:, None]*velocity
    energy = pressure*area/(GAMMA-1)+.5*mass*np.sum(velocity*velocity, axis=1)
    return mass, np.r_[sites.ravel(), weights, momentum.ravel(), energy, np.zeros(2)]


def rhs(y, mass):
    sites, weights, momentum, energy = unpack(y)
    cells = power_cells(sites, weights, DOMAIN); area = np.array([p.area for p in cells])
    if (area <= 0).any(): raise RuntimeError('empty power cell')
    velocity = momentum/mass[:, None]; rho = mass/area
    internal = energy-.5*np.sum(momentum*momentum, axis=1)/mass
    pressure = (GAMMA-1)*internal/area
    if (pressure <= 0).any() or not np.isfinite(pressure).all(): raise RuntimeError('nonpositive pressure')
    sound = np.sqrt(GAMMA*pressure/rho)
    U = np.column_stack((rho, momentum/area[:, None], energy/area))
    geometry = faces(sites, weights, cells); records = []
    laplacian = np.zeros((4, 4))
    for i, j, mid, normal, length in geometry:
        vi = velocity[i]; vni = np.dot(vi, normal)
        if j >= 0:
            vj = velocity[j]; vnj = np.dot(vj, normal); distance = np.linalg.norm(sites[j]-sites[i])
            base = (np.dot(sites[j], vj)-np.dot(sites[i], vi)-np.dot(vj-vi, mid))/distance
            density = (rho[i]+rho[j])/2; coefficient = length*density/(2*distance)
            laplacian[i, i] += coefficient; laplacian[j, j] += coefficient
            laplacian[i, j] -= coefficient; laplacian[j, i] -= coefficient
            wavespeed = max(sound[i]+abs(vni-base), sound[j]+abs(vnj-base))
        else:
            vnj = -vni; base = 0.; distance = 0.; density = rho[i]
            wavespeed = sound[i]+abs(vni)
        records.append([i, j, mid, normal, length, base, distance, density, vni, vnj, wavespeed])
    # The mass constraint is linear for fixed Rusanov wave speeds. Increase
    # them until they also bound the waves in the resulting moving-face frame.
    for iteration in range(100):
        residual = np.zeros(4)
        for i, j, mid, n, length, base, distance, density, vni, vnj, speed in records:
            if j < 0: continue  # reflective wall: exactly zero mass flux
            fm = .5*(rho[i]*vni+rho[j]*vnj)-density*base-.5*speed*(rho[j]-rho[i])
            residual[i] += length*fm; residual[j] -= length*fm
        wdot = np.r_[np.linalg.solve(laplacian[:3, :3], residual[:3]), 0.]
        bounded = True
        for r in records:
            i, j, _, _, _, base, distance, _, vni, vnj, speed = r
            if j < 0: continue
            wn = base+(wdot[i]-wdot[j])/(2*distance)
            required = max(sound[i]+abs(vni-wn), sound[j]+abs(vnj-wn))
            if required > speed*(1+1e-12): r[-1] = required; bounded = False
        if bounded: break
    else: raise RuntimeError('ALE wave-speed bound failed to converge')
    change = np.zeros((4, 4)); impulse = np.zeros(2); max_face_speed = 0.
    for i, j, mid, normal, length, base, distance, density, vni, vnj, speed in records:
        ui = U[i]; vi = velocity[i]
        if j >= 0:
            uj = U[j]; vj = velocity[j]; pj = pressure[j]
            wn = base+(wdot[i]-wdot[j])/(2*distance)
        else:
            vj = vi-2*vni*normal; pj = pressure[i]
            uj = np.r_[rho[i], rho[i]*vj, energy[i]/area[i]]; wn = 0.
        fi = np.r_[rho[i]*vni, rho[i]*vi*vni+pressure[i]*normal, (ui[3]+pressure[i])*vni]
        fj = np.r_[uj[0]*vnj, uj[1:3]*vnj+pj*normal, (uj[3]+pj)*vnj]
        flux = (.5*(fi+fj)-wn*.5*(ui+uj)-.5*speed*(uj-ui))*length
        change[i] -= flux
        if j >= 0: change[j] += flux
        else: impulse -= flux[1:3]
        max_face_speed = max(max_face_speed, abs(wn))
    if np.max(abs(change[:, 0])) > 1e-9*mass.sum(): raise RuntimeError('fixed-mass flux constraint failed')
    out = np.r_[velocity.ravel(), wdot, change[:, 1:3].ravel(), change[:, 3], impulse]
    return out, dict(mass_flux=float(np.max(abs(change[:, 0]))), pressure_min=float(pressure.min()),
                     pressure_max=float(pressure.max()), area_min=float(area.min()),
                     wave_iterations=iteration+1, face_speed=max_face_speed)


def evolve(kind, dt, end=.2, rhs_fn=rhs, initial_fn=initial):
    mass, y = initial_fn(kind); y0 = y.copy(); initial_energy=y[20:24].sum()
    diagnostics = dict(mass_flux=0., pressure_min=float('inf'), area_min=float('inf'), wave_iterations=0)
    steps = round(end/dt)
    if abs(steps*dt-end) > 1e-12: raise ValueError('end must be multiple of dt')
    for _ in range(steps):
        k1, d1 = rhs_fn(y, mass); k2, d2 = rhs_fn(y+dt*k1/2, mass)
        k3, d3 = rhs_fn(y+dt*k2/2, mass); k4, d4 = rhs_fn(y+dt*k3, mass)
        y += dt*(k1+2*k2+2*k3+k4)/6
        for d in (d1,d2,d3,d4):
            diagnostics['mass_flux'] = max(diagnostics['mass_flux'], d['mass_flux'])
            diagnostics['pressure_min'] = min(diagnostics['pressure_min'], d['pressure_min'])
            diagnostics['area_min'] = min(diagnostics['area_min'], d['area_min'])
            diagnostics['wave_iterations'] = max(diagnostics['wave_iterations'], d['wave_iterations'])
    _, dfinal = rhs_fn(y, mass)
    sites, w, momentum, energy = unpack(y)
    return dict(kind=kind, dt=dt, end=end, **diagnostics,
                energy_relative=float(energy.sum()/initial_energy-1),
                momentum_ledger_error=float(np.linalg.norm(momentum.sum(0)-y0[12:20].reshape(4,2).sum(0)-y[24:26])),
                final_pressure_range=[dfinal['pressure_min'],dfinal['pressure_max']],
                generator_displacement=float(np.linalg.norm(sites-y0[:8].reshape(4,2))),
                weight_change=float(np.linalg.norm(w-y0[8:12])), final_state=y.tolist())


if __name__ == '__main__':
    for kind in ('uniform', 'contact', 'acoustic'):
        for dt in (.01, .005, .0025, .00125):
            print(json.dumps(evolve(kind, dt)), flush=True)

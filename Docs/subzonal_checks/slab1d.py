"""1D slab (mirror wall at x=0, free surface at x0) expanding into vacuum with a rigid cap
of half-width L/2 = 2.5 Delta attached to the end particle.  Models: 'cell' (m e(m/V), V of the
capped cell) and 'sub' (staggered subzonal, cap sub-zone rigid).  Exact: centred rarefaction,
free-surface speed 2c0/(gamma-1) = 3 c0.  Checks: edge-particle velocity, density L1, dE/E.
"""
import numpy as np
g = 5.0/3.0

def minmod(a, b):
    return np.where(a*b > 0, np.sign(a)*np.minimum(np.abs(a), np.abs(b)), 0.0)

def run(model, N=100, x0=1.0, tend=0.5, Lfac=5.0, useQ=True, cfl=0.2):
    D = x0/N; L = Lfac*D
    x = (np.arange(N) + 0.5)*D; v = np.zeros(N); m = np.full(N, D)   # rho0 = 1
    K = np.ones(N); rho0 = 1.0; P0 = 1.0; c0 = np.sqrt(g*P0/rho0)

    def gaps(x):
        dL = np.empty(N); dR = np.empty(N)
        dL[0] = 2*x[0]; dL[1:] = x[1:] - x[:-1]
        dR[:-1] = x[1:] - x[:-1]; dR[-1] = L
        return dL, dR

    def U_G(x, K):
        dL, dR = gaps(x)
        if model == 'cell':
            V = 0.5*(dL + dR); G = m*(m/V)**(g-1)/(g-1)
        else:
            G = 0.5*m*((m/dL)**(g-1) + (m/dR)**(g-1))/(g-1)
        return K*G, G

    def rhs(x, v, K):
        dL, dR = gaps(x)
        if model == 'cell':
            V = 0.5*(dL + dR); P = K*(m/V)**g
            F = np.zeros(N)
            # F_i = sum_j P_j dV_j/dx_i : dV_{i-1}/dx_i = 1/2, dV_{i+1}/dx_i = -1/2, dV_i/dx_i: 0 interior
            F[1:] += 0.5*P[:-1]; F[:-1] -= 0.5*P[1:]
            F[0] += 0.5*P[0]      # mirror: V_0 = (x_0 + x_1)/2
            F[-1] += 0.5*P[-1]    # cap: V_{N-1} = (x_{N-1}-x_{N-2})/2 + L/2
        else:
            V = 0.5*(dL + dR)
            # gap pressures: gap i+1/2 between i and i+1, rho = m/d
            d = x[1:] - x[:-1]; rg = m[:-1]/d
            Phat = 0.5*(K[:-1]*rg**g + K[1:]*rg**g)
            F = np.zeros(N)
            F[:-1] -= Phat; F[1:] += Phat
            F[0] += K[0]*(m[0]/(2*x[0]))**g    # mirror gap, d = 2 x_0, dd/dx_0 = 2, times 1/2
            # cap gap: d_R = L fixed -> no force
        rho = m/V; Pc = K*rho**g; c = np.sqrt(g*Pc/rho)
        HQ = np.zeros(N)
        if useQ:
            du_cell = v[:-1] - v[1:]; dxf = x[1:] - x[:-1]
            s_f = -du_cell/dxf
            slope = np.zeros(N); slope[1:-1] = minmod(s_f[:-1], s_f[1:])
            xf = 0.5*(x[:-1] + x[1:])
            uL = v[:-1] + slope[:-1]*(xf - x[:-1]); uR = v[1:] + slope[1:]*(xf - x[1:])
            mm = minmod(uL - uR, du_cell)
            rhob = 2*rho[:-1]*rho[1:]/(rho[:-1] + rho[1:]); cb = 0.5*(c[:-1] + c[1:])
            a = 0.125*(g+1)*np.abs(mm); Z = 0.5*rhob*(a + np.sqrt(a*a + cb*cb))
            Q = Z*mm
            F[:-1] -= Q; F[1:] += Q
            HQ[:-1] += 0.5*Q*du_cell; HQ[1:] += 0.5*Q*du_cell
        U, G = U_G(x, K)
        return v, F/m, HQ/G, (rho, c)

    t = 0.0; E0 = np.sum(0.5*m*v**2) + U_G(x, K)[0].sum()
    nstep = 0
    while t < tend - 1e-14:
        _, _, _, (rho, c) = rhs(x, v, K)
        dL, dR = gaps(x)
        dt = cfl*np.min(np.minimum(dL, dR)/(c + np.abs(np.gradient(v)) + 1e-30)); dt = min(dt, tend - t)
        def f(x, v, K):
            a, b, cc, _ = rhs(x, v, K); return a, b, cc
        k1 = f(x, v, K); k2 = f(x+0.5*dt*k1[0], v+0.5*dt*k1[1], K+0.5*dt*k1[2])
        k3 = f(x+0.5*dt*k2[0], v+0.5*dt*k2[1], K+0.5*dt*k2[2]); k4 = f(x+dt*k3[0], v+dt*k3[1], K+dt*k3[2])
        x = x + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); v = v + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        K = K + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]); t += dt; nstep += 1
        if np.any(np.diff(x) <= 0): print("collision"); break
    E = np.sum(0.5*m*v**2) + U_G(x, K)[0].sum()
    dL, dR = gaps(x); V = 0.5*(dL + dR); rho = m/V
    # exact rarefaction (t < x0/c0): head at x0 - c0 t, tail (free surface) at x0 + 3 c0 t
    xh = x0 - c0*t
    u_ex = np.where(x < xh, 0.0, 2/(g+1)*((x - x0)/t + c0))
    c_ex = c0 - 0.5*(g-1)*u_ex; rho_ex = rho0*(np.maximum(c_ex, 0)/c0)**(2/(g-1))
    # Lagrangian mass coordinate of the particle: exact density at its position compared on interior
    sel = x < x0 + 3*c0*t
    l1 = np.sum(np.abs(rho[sel] - rho_ex[sel]))/np.sum(np.abs(rho_ex[sel]))
    l1v = np.sum(np.abs(v[sel] - u_ex[sel]))/np.sum(np.abs(u_ex[sel]) + 1e-30)
    # exact velocity of the edge mass element: the last particle holds mass m centred at x0 - D/2;
    # the exact mass element at Lagrangian depth D/2 ... compare to the free-surface speed instead
    print("%4s Q=%d N=%d: steps=%d  v_edge/3c0=%.3f  x_edge=%.3f (fs exact %.3f)  rho_edge/rho0=%.3f  L1rho=%.3e L1v=%.3e  dE/E=%.1e  Kmax-1=%.1e"
          % (model, useQ, N, nstep, v[-1]/(3*c0), x[-1], x0 + 3*c0*t, rho[-1], l1, l1v, (E-E0)/E0, K.max()-1))
    return x, v, rho, K

if __name__ == '__main__':
    for N in (50, 100, 200):
        for model in ('cell', 'sub'):
            run(model, N=N, useQ=False)
    for model in ('cell', 'sub'):
        run(model, N=100, useQ=True)
    # longer time, does the edge run away?
    for model in ('cell', 'sub'):
        run(model, N=100, useQ=False, tend=0.7)

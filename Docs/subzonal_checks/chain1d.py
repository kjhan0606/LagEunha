"""1D Lagrangian Voronoi chain, gamma=5/3.
Variants: 'cell' (a), 'sub' (b: (나) staggered vN form), 'lin' (c: (나') linearized).
Dissipation: Q = Z*minmod(du_rec, du_cell), Z = rhobar*(cbar + beta*|mm|) (beta=0 linear, 1 quadratic).
State: x, v, K (entropy variable); U_i = K_i G_i(x). Semi-discrete total energy exact.
Boundaries: ghost copies of end cells (velocity mirrored for 'wall' at x=0 in Noh via symmetric domain).
"""
import numpy as np
g = 5.0/3.0

def minmod(a, b):
    return np.where(a*b > 0, np.sign(a)*np.minimum(np.abs(a), np.abs(b)), 0.0)

class Chain:
    def __init__(self, x, v, rho, P, m, variant='cell', beta=1.0, kappa=1.0, bc='ghost'):
        self.x = x.astype(float).copy(); self.v = v.astype(float).copy(); self.m = m.astype(float).copy()
        self.K = P/rho**g
        self.variant = variant; self.beta = beta; self.kappa = kappa; self.bc = bc
        self.N = len(x)

    # geometry with ghost ends: pad x by reflection of spacing
    def geom(self, x):
        N = self.N
        xp = np.empty(N+2); xp[1:-1] = x
        xp[0] = 2*x[0]-x[1]; xp[-1] = 2*x[-1]-x[-2]
        d = np.diff(xp)             # N+1 gaps: d[i] = gap between cell i-1 and i (ghost at i=0), d[i+1] right gap of i
        dL = d[:-1]; dR = d[1:]
        V = 0.5*(dL+dR)
        return dL, dR, V

    def G_and_grad(self, x):
        """Return G_i (so U_i = K_i G_i) and dG_i/d(dL_i), dG_i/d(dR_i)."""
        m = self.m; dL, dR, V = self.geom(x)
        gm1 = g-1
        if self.variant == 'cell':
            G = m*(m/V)**gm1/gm1
            dGdV = -m*(m/V)**gm1/V      # = -P/K * ... : d/dV [m rho^{g-1}/(g-1)] = -(m/V)^g
            return G, 0.5*dGdV, 0.5*dGdV
        if self.variant == 'sub':
            # U_i = K [ (m/2) rho_L^{g-1} + (m/2) rho_R^{g-1} ]/(g-1), rho_f = m/d_f
            GL = 0.5*m*(m/dL)**gm1/gm1; GR = 0.5*m*(m/dR)**gm1/gm1
            dGL = -0.5*(m/dL)**g; dGR = -0.5*(m/dR)**g
            return GL+GR, dGL, dGR
        if self.variant == 'lin':
            # U_i = K m rho^{g-1}/(g-1) + K * (kappa/2) g rho^{g-1} (m/2) [ (V/dL-1)^2 + (V/dR-1)^2 ]
            rho = m/V; k = self.kappa
            G0 = m*rho**gm1/gm1
            aL = V/dL-1; aR = V/dR-1
            S = aL**2 + aR**2
            G1 = 0.5*k*g*rho**gm1*0.5*m*S
            G = G0+G1
            # derivatives wrt dL, dR (V depends on both with 1/2)
            dG0dV = -(m/V)**g
            dG1dV = 0.5*k*g*0.5*m*( -gm1*rho**gm1/V*S + rho**gm1*2*(aL/dL + aR/dR) )
            dG1_dL_expl = 0.5*k*g*rho**gm1*0.5*m*2*aL*(-V/dL**2)
            dG1_dR_expl = 0.5*k*g*rho**gm1*0.5*m*2*aR*(-V/dR**2)
            dGdL = 0.5*(dG0dV+dG1dV) + dG1_dL_expl
            dGdR = 0.5*(dG0dV+dG1dV) + dG1_dR_expl
            return G, dGdL, dGdR
        raise ValueError

    def rhs(self, x, v, K):
        N = self.N; m = self.m
        dL, dR, V = self.geom(x)
        G, dGdL, dGdR = self.G_and_grad(x)
        # dU_i/d dL_i = K dGdL; dL_i = x_i - x_{i-1}; dR_i = x_{i+1}-x_i
        # F_a = -sum_i K_i [dGdL_i d(dL_i)/dx_a + dGdR_i d(dR_i)/dx_a]
        fL = K*dGdL; fR = K*dGdR
        F = -(fL - fR)                      # own-cell terms: d dL_i/dx_i = +1, d dR_i/dx_i = -1
        F[:-1] += fL[1:]                    # dL_{i+1}/dx_i = -1 -> -(-1) fL_{i+1}
        F[1:]  -= fR[:-1]                   # dR_{i-1}/dx_i = +1 -> -(+1) fR_{i-1}
        # mirror ghosts: dL_0 = dR_0 = x1-x0 ; dR_{N-1} = dL_{N-1} = x_{N-1}-x_{N-2}
        F[0] += 2*fL[0]; F[1] -= fL[0]
        F[-1] -= 2*fR[-1]; F[-2] += fR[-1]
        # cell state for Q
        rho = m/V; P0 = K*rho**g; c = np.sqrt(g*P0/rho)
        # reconstruction (minmod slopes) at faces i+1/2, i=0..N-2
        du_f = np.diff(v); dx_f = np.diff(x)
        s_f = du_f/dx_f
        slope = np.zeros(N); slope[1:-1] = minmod(s_f[:-1], s_f[1:])
        xf = 0.5*(x[:-1]+x[1:])
        uLf = v[:-1] + slope[:-1]*(xf-x[:-1]); uRf = v[1:] + slope[1:]*(xf-x[1:])
        du_rec = uLf - uRf; du_cell = v[:-1]-v[1:]
        mm = minmod(du_rec, du_cell)
        rhob = 2*rho[:-1]*rho[1:]/(rho[:-1]+rho[1:]); cb = 0.5*(c[:-1]+c[1:])
        if self.beta == 'W':   # Wilkins / exact two-shock Hugoniot impedance
            a = 0.25*(g+1)*np.abs(mm)
            Z = rhob*(a + np.sqrt(a*a + cb*cb))
        else:
            Z = rhob*(cb + self.beta*np.abs(mm))
        Q = Z*mm
        F[:-1] -= Q; F[1:] += Q
        HQ = np.zeros(N); HQ[:-1] += 0.5*Q*du_cell; HQ[1:] += 0.5*Q*du_cell
        Kdot = HQ/G
        if self.bc == 'fixed_ends':
            F[0] = 0; F[-1] = 0
        return v, F/m, Kdot, (V, rho, P0, c, dL, dR, Q, du_cell)

    def energy(self, x, v, K):
        G, _, _ = self.G_and_grad(x)
        return np.sum(0.5*self.m*v**2) + np.sum(K*G)

    def dt(self, x, v, K, cfl):
        _, _, _, aux = self.rhs(x, v, K)
        V, rho, P0, c, dL, dR, Q, du = aux
        dmin = np.minimum(dL, dR)
        # local omega_max ~ 2c/dmin (Nyquist) ; include |du|
        return cfl*np.min(dmin/(c + np.abs(v - np.roll(v, 1)) + 1e-30))

    def run(self, tend, cfl=0.25, record=None):
        x, v, K = self.x, self.v, self.K
        t = 0.0; E0 = self.energy(x, v, K); nstep = 0
        mingap_hist = []
        self.dtmin = np.inf; self.aborted = None
        while t < tend - 1e-14:
            h = min(self.dt(x, v, K, cfl), tend-t)
            self.dtmin = min(self.dtmin, h)
            if h < 1e-8 or nstep > 200000:
                self.aborted = "dt=%.1e at t=%.4f step %d" % (h, t, nstep); print("  ABORT", self.aborted); break
            def f(x, v, K):
                dx, dv, dK, _ = self.rhs(x, v, K); return dx, dv, dK
            k1 = f(x, v, K)
            k2 = f(x+0.5*h*k1[0], v+0.5*h*k1[1], K+0.5*h*k1[2])
            k3 = f(x+0.5*h*k2[0], v+0.5*h*k2[1], K+0.5*h*k2[2])
            k4 = f(x+h*k3[0], v+h*k3[1], K+h*k3[2])
            x = x + h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
            v = v + h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
            K = K + h/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
            t += h; nstep += 1
            if np.any(np.diff(x) <= 0):
                print("  COLLISION at t=%.4f step %d" % (t, nstep)); break
            dL, dR, V = self.geom(x)
            mingap_hist.append((t, np.min(np.diff(x)/(0.5*(V[:-1]+V[1:])))))
        self.x, self.v, self.K = x, v, K
        self.t = t; self.nstep = nstep; self.dE = (self.energy(x, v, K)-E0)/abs(E0)
        self.mingap_hist = np.array(mingap_hist)
        return x, v, K

    def state(self):
        dL, dR, V = self.geom(self.x)
        rho = self.m/V; P0 = self.K*rho**g
        return V, rho, P0, dL, dR

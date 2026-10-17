"""Audit of Docs/design_astra.md in 1D, implemented exactly as specified:
  reversible: cell-volume Lagrangian F_a = sum_i P_i dV_i/dx_a (no Phi)
  Q        : QLAW=2, Q = Z m^2/d_u, m = minmod(d_r, d_u), Z = max(Z_eff(Davis, analytic), Z_H), face subzonal states
  hourglass: xi_if = ln(B_if/W_if) (1D: ln(V_i/d_f)), r_if = xidot (2-hop affine filter is the identity in 1D, verified),
             C(q) = rho_if h W [b q + sqrt(b^2 q^2 + c_if^2)], b = (g+1) h/4, F = -sum C r g^T, H_i = sum C r^2
  dt       : min( 0.25 min h_if/(c_if+|d_u|), 0.25/(Omega + Lambda_D) )
State (x, v, K), U_i = K_i G_i, RK4. bc in {'periodic', 'ghost', 'fixed_ends'}.
"""
import numpy as np
g = 5.0/3.0


def minmod(a, b):
    return np.where(a*b > 0, np.sign(a)*np.minimum(np.abs(a), np.abs(b)), 0.0)


class Design1D:
    def __init__(self, x, v, rho, P, m, bc='periodic', Lbox=None, hg='hug', hgscale=1.0, hgcoef=1.0,
                 qlaw=2, zface=True, zh=True, qon=True, dtmode='design', cfl=0.25, zmode='design'):
        self.x = x.astype(float).copy(); self.v = v.astype(float).copy(); self.m = m.astype(float).copy()
        self.K = P/rho**g; self.N = len(x); self.bc = bc; self.L = Lbox
        self.hg = hg; self.hgscale = hgscale; self.hgcoef = hgcoef
        self.qlaw = qlaw; self.zface = zface; self.zh = zh; self.qon = qon
        self.dtmode = dtmode; self.cfl = cfl; self.zmode = zmode
        self.cum = dict(HQ=0.0, Hhg=0.0)

    # ---------------- geometry ----------------
    def faces(self, x):
        N = self.N
        if self.bc == 'periodic':
            d = np.roll(x, -1) - x; d[-1] += self.L          # face f: (f, f+1 mod N)
            dR = d; dL = np.roll(d, 1)
            cR = np.full(N, 0.5); cL = np.full(N, 0.5)
        else:
            d = np.diff(x)                                   # N-1 faces
            dR = np.empty(N); dL = np.empty(N)
            dR[:-1] = d; dL[1:] = d
            dR[-1] = d[-1]; dL[0] = d[0]                     # mirrored ghost gaps
            cR = np.full(N, 0.5); cL = np.full(N, 0.5)
            cR[-1] = 0.0; cL[0] = 0.0                        # end cells: V = d (own gap counted once)
            cR[0] = 1.0; cL[-1] = 1.0
        V = 0.5*(dL + dR) if self.bc == 'periodic' else np.where(np.arange(N) == 0, dR, np.where(np.arange(N) == N-1, dL, 0.5*(dL+dR)))
        return d, dL, dR, V, cL, cR

    def G(self, x):
        _, _, _, V, _, _ = self.faces(x)
        return self.m*(self.m/V)**(g-1)/(g-1)

    def energy(self, x, v, K):
        return np.sum(0.5*self.m*v**2) + np.sum(K*self.G(x))

    # ---------------- rhs ----------------
    def rhs(self, x, v, K, want_dt=False):
        N = self.N; m = self.m
        d, dL, dR, V, cL, cR = self.faces(x)
        per = self.bc == 'periodic'
        Nf = N if per else N-1
        fi = np.arange(Nf); fj = (fi + 1) % N
        rho = m/V; P = K*rho**g; c = np.sqrt(g*P/rho)
        # reversible: dV_i/dd_f = cR_i (face is right face of i), dV_j/dd_f = cL_j
        pf = P[fi]*cR[fi] + P[fj]*cL[fj]
        F = np.zeros(N)
        np.add.at(F, fj, pf); np.add.at(F, fi, -pf)
        # ---- face subzonal states
        rho_iR = m[fi]/d; rho_jL = m[fj]/d
        P_iR = K[fi]*rho_iR**g; P_jL = K[fj]*rho_jL**g
        c_iR = np.sqrt(g*P_iR/rho_iR); c_jL = np.sqrt(g*P_jL/rho_jL)
        # ---- Q
        du = v[fi] - v[fj]
        # minmod slopes per cell from adjacent face slopes
        sf = -du/d                                            # (v_j - v_i)/d
        slope = np.zeros(N)
        if per:
            slope = minmod(np.roll(sf, 1), sf)
        else:
            slope[1:-1] = minmod(sf[:-1], sf[1:])
        uL = v[fi] + slope[fi]*0.5*d; uR = v[fj] - slope[fj]*0.5*d
        dr = uL - uR
        mm = minmod(dr, du)
        if self.zface:
            rL, pL, cl = rho_iR, P_iR, c_iR; rR, pR, cr = rho_jL, P_jL, c_jL
        else:
            rL, pL, cl = rho[fi], P[fi], c[fi]; rR, pR, cr = rho[fj], P[fj], c[fj]
        mu = rL*rR/(rL + rR); Pc = (rR*pL + rL*pR)/(rL + rR); cM = np.maximum(cl, cr)
        # analytic Davis Z_eff with delta = d_r
        Zeff = np.where(dr >= 0, mu*(cM + dr),
                        np.where(Pc + mu*cM*dr > 0, mu*cM, -Pc/np.where(dr == 0, 1.0, dr)))
        aH = (g+1)/8*np.abs(mm)
        ZH = mu*(aH + np.sqrt(aH*aH + cM*cM))
        Z = np.maximum(Zeff, ZH) if self.zh else Zeff
        if self.zmode == 'ZF':      # Fable's full-amplitude two-shock Hugoniot: rhobar [a + sqrt(a^2+c^2)], a = (g+1)|m|/4
            aF = (g+1)/4*np.abs(mm); Z = 2*mu*(aF + np.sqrt(aF*aF + cM*cM))
        elif self.zmode == 'ZH2':   # Z_H with rho_H instead of rho_H/2
            Z = np.maximum(Zeff, 2*ZH)
        elif self.zmode == 'davis2':
            Z = 2*Zeff
        if self.qlaw == 2:
            Q = np.where(du != 0, Z*mm*mm/np.where(du == 0, 1.0, du), 0.0)
        else:
            Q = Z*mm
        if not self.qon:
            Q = 0*Q
        np.add.at(F, fi, -Q); np.add.at(F, fj, Q)
        HQ = np.zeros(N); np.add.at(HQ, fi, 0.5*Q*du); np.add.at(HQ, fj, 0.5*Q*du)
        # ---- hourglass (cells with two real faces)
        H = np.zeros(N); Fh = np.zeros(N)
        Dhg_row = np.zeros(N)
        if self.hg is not None:
            if per:
                cells = np.arange(N)
            else:
                cells = np.arange(1, N-1)
            ip = (cells + 1) % N; im = (cells - 1) % N
            dRc = dR[cells]; dLc = dL[cells]; Vc = V[cells]
            dRdot = v[ip] - v[cells]; dLdot = v[cells] - v[im]
            Vdot = 0.5*(dLdot + dRdot)
            rR_ = Vdot/Vc - dRdot/dRc; rL_ = Vdot/Vc - dLdot/dLc
            # g vectors: (coef on i-1, i, i+1)
            gR = (-0.5/Vc, 1.0/dRc, 0.5/Vc - 1.0/dRc)
            gL = (-0.5/Vc + 1.0/dLc, -1.0/dLc, 0.5/Vc)
            out = []
            for dd, rr, gg in ((dRc, rR_, gR), (dLc, rL_, gL)):
                h = 0.5*dd; W = 0.5*dd; rhof = m[cells]/dd; cf = np.sqrt(g*K[cells]*rhof**(g-1))
                q = np.abs(rr)
                if self.hg == 'hug':
                    b = self.hgscale*(g+1)/4*h
                    s = np.sqrt(b*b*q*q + cf*cf)
                    C = self.hgcoef*rhof*h*W*(b*q + s)
                    Cp = self.hgcoef*rhof*h*W*(2*b*q + s + b*b*q*q/s)     # C + q C'
                else:
                    C = self.hgcoef*rhof*cf*h*W; Cp = C
                qq = C*rr
                np.add.at(Fh, im, -qq*gg[0]); np.add.at(Fh, cells, -qq*gg[1]); np.add.at(Fh, ip, -qq*gg[2])
                H[cells] += C*rr*rr
                # row sums of D_hg = sum (C+qC') g g^T  (|g_a| sum_b |g_b|)
                gsum = np.abs(gg[0]) + np.abs(gg[1]) + np.abs(gg[2])
                np.add.at(Dhg_row, im, Cp*np.abs(gg[0])*gsum); np.add.at(Dhg_row, cells, Cp*np.abs(gg[1])*gsum)
                np.add.at(Dhg_row, ip, Cp*np.abs(gg[2])*gsum)
            F += Fh
        G = self.m*(self.m/V)**(g-1)/(g-1)
        Kdot = (HQ + H)/G
        if self.bc == 'fixed_ends':
            F[0] = 0.0; F[-1] = 0.0
        dtinfo = None
        if want_dt:
            # geometric/acoustic limit: 0.25 min h_if/(c_if + |d_u|) over both sides of each face
            hf = 0.5*d
            dt_geom = 0.25*np.min(np.minimum(hf/(c_iR + np.abs(du)), hf/(c_jL + np.abs(du))))
            # reversible stiffness row sums: K_ab = -dF_a/dx_b ; |K| row sum / sqrt(m_a m_b)
            # dP_i/dd_f = -g P_i/V_i * dV_i/dd_f ; F_j += pf, F_i -= pf  => entries +-(dpf/dd) with dd/dx = +-1
            dpf = g*(P[fi]/V[fi]*cR[fi]**2 + P[fj]/V[fj]*cL[fj]**2)   # |d pf/d d_f| dominant own-gap part
            # cross terms: pf also depends on neighbouring gaps through V_i (dL_i) and V_j (dR_j): same magnitude scale
            Krow = np.zeros(N)
            np.add.at(Krow, fi, 2*dpf/np.sqrt(m[fi]*m[fi])); np.add.at(Krow, fj, 2*dpf/np.sqrt(m[fj]*m[fj]))
            Omega = np.sqrt(Krow.max())
            # damping: hg rows + Q rows
            dQdu = np.where(np.abs(du) <= np.abs(dr), Z + np.abs(du)*mu*(g+1)/4, Z)   # |dQ/d d_u|
            dQdr = np.where(np.abs(du) <= np.abs(dr), 0.0, 2*Z + mu*np.abs(dr))        # |dQ/d d_r|
            DQ = (dQdu + 1.5*dQdr)*2.0                                                  # row sum of 2x2 block (+ d_r spread)
            Drow = Dhg_row.copy()
            np.add.at(Drow, fi, DQ); np.add.at(Drow, fj, DQ)
            Lam = (Drow/m).max()
            dt_spec = 0.25/(Omega + Lam)
            dtinfo = (dt_geom, dt_spec, Omega, Lam)
        return v, F/m, Kdot, dict(d=d, V=V, rho=rho, P=P, c=c, Q=Q, du=du, HQ=HQ, H=H, Fh=Fh, dt=dtinfo, mm=mm, dr=dr, Z=Z)

    # ---------------- integrator ----------------
    def run(self, tend, record_every=0, maxsteps=400000):
        x, v, K = self.x, self.v, self.K
        t = 0.0; E0 = self.energy(x, v, K); nstep = 0
        self.hist = []; self.aborted = None; self.collided = False
        self.dt_hist = []; self.limiter = {'geom': 0, 'spec': 0}
        while t < tend - 1e-14:
            _, _, _, aux = self.rhs(x, v, K, want_dt=True)
            dtg, dts, Om, Lam = aux['dt']
            if self.dtmode == 'design':
                dt = min(dtg, dts); self.limiter['geom' if dtg <= dts else 'spec'] += 1
            elif self.dtmode == 'geom':
                dt = dtg; self.limiter['geom'] += 1
            elif isinstance(self.dtmode, float):      # spectral coefficient instead of 0.25
                dts2 = dts*self.dtmode/0.25
                dt = min(dtg, dts2); self.limiter['geom' if dtg <= dts2 else 'spec'] += 1
            else:   # legacy chain1d-like: cfl * min gap/(c+|du|)
                dt = self.cfl*np.min(aux['d']/(aux['c'][:len(aux['d'])] + np.abs(aux['du'])))
            dt = min(dt, tend - t)
            self.dt_hist.append((t, dt))
            if dt < 1e-12 or nstep > maxsteps:
                self.aborted = "dt=%.1e at t=%.4f step %d" % (dt, t, nstep); break
            def f(x, v, K):
                dx, dv, dK, a = self.rhs(x, v, K); return dx, dv, dK, a
            k1 = f(x, v, K)
            k2 = f(x+0.5*dt*k1[0], v+0.5*dt*k1[1], K+0.5*dt*k1[2])
            k3 = f(x+0.5*dt*k2[0], v+0.5*dt*k2[1], K+0.5*dt*k2[2])
            k4 = f(x+dt*k3[0], v+dt*k3[1], K+dt*k3[2])
            x = x + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
            v = v + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
            K = K + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
            for kk in (k1, k2, k3, k4):
                w = dt/6*(2 if kk is k2 or kk is k3 else 1)
                self.cum['HQ'] += w*kk[3]['HQ'].sum(); self.cum['Hhg'] += w*kk[3]['H'].sum()
            t += dt; nstep += 1
            d, dL, dR, V, _, _ = self.faces(x)
            if np.any(d <= 0) or np.any(K <= 0):
                self.collided = True; self.aborted = "COLLISION at t=%.4f step %d" % (t, nstep); break
            if self.bc == 'periodic':
                hloc = 0.5*(V + np.roll(V, -1))
            else:
                hloc = 0.5*(V[:-1] + V[1:])
            self.hist.append((t, (d/hloc).min(), dt, np.sum(0.5*self.m*v**2)))
        self.x, self.v, self.K = x, v, K; self.t = t; self.nstep = nstep
        self.dE = (self.energy(x, v, K) - E0)/abs(E0)
        self.hist = np.array(self.hist)
        return x, v, K

    def summary(self):
        h = self.hist
        n = len(h); q = max(1, n//10)
        dt0 = np.median(h[:q, 2]); dt1 = np.median(h[-q:, 2])
        return dict(steps=self.nstep, t=self.t, mingap=h[:, 1].min(), dt_first=dt0, dt_last=dt1,
                    dtmin=h[:, 2].min(), dE=self.dE, HQ=self.cum['HQ'], Hhg=self.cum['Hhg'],
                    lim=self.limiter, aborted=self.aborted)


def filter_check():
    """Explicit 2-hop weighted LSQ affine filter in 1D: r_filtered == r_raw, L^T == g^T (P^T g = 0)."""
    rng = np.random.default_rng(1)
    N = 41; x = np.sort(np.arange(N) + rng.uniform(-0.3, 0.3, N)); v = rng.normal(size=N)
    out = []
    for i in range(2, N-2):
        d = np.diff(x); dL = d[i-1]; dR = d[i]; V = 0.5*(dL+dR)
        g_ = {i-1: -0.5/V, i: 1.0/dR, i+1: 0.5/V - 1.0/dR}          # xi_R
        w1m = 1/dL**2; w1p = 1/dR**2
        st = {i-2: V*w1m/d[i-2]**2, i-1: w1m, i: w1m+w1p, i+1: w1p, i+2: V*w1p/d[i+1]**2}
        ks = sorted(st); W = np.array([st[k] for k in ks]); X = np.array([[1.0, x[k]-x[i]] for k in ks])
        A = X.T @ (W[:, None]*X); Pm = X @ np.linalg.solve(A, X.T*W)   # fit operator on stencil
        gv = np.array([g_.get(k, 0.0) for k in ks]); vv = v[ks]
        r_raw = gv @ vv; r_filt = gv @ (vv - Pm @ vv)
        LT = gv - Pm.T @ gv
        out.append((abs(r_filt - r_raw), np.abs(LT - gv).max(), abs(r_raw)))
    out = np.array(out)
    print("1D 2-hop affine filter: max|r_filt - r_raw| = %.1e, max|L^T - g^T| = %.1e  (max|r| = %.2f)" %
          (out[:, 0].max(), out[:, 1].max(), out[:, 2].max()))


if __name__ == '__main__':
    filter_check()

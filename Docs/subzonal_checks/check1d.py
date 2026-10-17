import numpy as np
g = 5.0/3.0
K = 1.0
def e(rho): return K*rho**(g-1)/(g-1)
def P(rho): return K*rho**g

N = 16; L = 16.0; D = L/N; m = 1.0  # rho=1
def gaps(x):
    xs = np.append(x, x[0]+L)
    return np.diff(xs)
def U_sub(x):   # option (나): 1D centroid = cell midpoint -> m_{i,f} = m/2 exactly
    d = gaps(x); dl = np.roll(d, 1)
    return np.sum(0.5*m*e(m/dl) + 0.5*m*e(m/d))
def U_cell(x):
    d = gaps(x); V = 0.5*(d+np.roll(d,1))
    return np.sum(m*e(m/V))
def Phi(x):
    d = gaps(x); dl = np.roll(d,1); V = 0.5*(d+dl)
    return 0.5*m*e(m/dl)+0.5*m*e(m/d) - m*e(m/V)
def hess(U, x, h=1e-4):
    n = len(x); H = np.zeros((n,n))
    for a in range(n):
        for b in range(n):
            xpp = x.copy(); xpp[a]+=h; xpp[b]+=h
            xpm = x.copy(); xpm[a]+=h; xpm[b]-=h
            xmp = x.copy(); xmp[a]-=h; xmp[b]+=h
            xmm = x.copy(); xmm[a]-=h; xmm[b]-=h
            H[a,b] = (U(xpp)-U(xpm)-U(xmp)+U(xmm))/(4*h*h)
    return H
x0 = np.arange(N)*D
for name, U in [("cell-volume", U_cell), ("subzonal (나)", U_sub)]:
    H = hess(U, x0)
    w, v = np.linalg.eigh(H)
    nyq = np.array([(-1)**i for i in range(N)])/np.sqrt(N)
    k_nyq = nyq @ H @ nyq
    print(f"{name}: Hessian eigen min/max = {w[0]:.4e} {w[-1]:.4e}; Nyquist stiffness = {k_nyq:.4e}")
rho = m/D; c2 = g*P(rho)/rho
print("prediction 4 rho c^2 / D =", 4*rho*c2/D, "  (omega_max = 2c/D =", 2*np.sqrt(c2)/D, ")")
# dispersion check of (나): k_mode vs 4 rho c^2/D sin^2(kD/2)
H = hess(U_sub, x0)
for kk in [1,2,4,8]:
    q = 2*np.pi*kk/L
    vec = np.cos(q*x0); vec/=np.linalg.norm(vec)
    print(f"k={kk}: stiffness {vec@H@vec:.4e}  pred {4*rho*c2/D*np.sin(q*D/2)**2:.4e}")
# Phi >= 0 and =0 for uniform, random configuration
rng = np.random.default_rng(1)
xr = np.sort(x0 + rng.uniform(-0.3,0.3,N)*D)
print("Phi uniform:", Phi(x0).min(), " Phi random min:", Phi(xr).min(), "max:", Phi(xr).max(), " Phi/U:", Phi(xr).sum()/U_cell(xr))
# force check: F_i = Phat_{i-1/2} - Phat_{i+1/2}
def F_fd(U,x,h=1e-6):
    return np.array([-(U(np.where(np.arange(N)==a, x+h, x))-U(np.where(np.arange(N)==a, x-h, x)))/(2*h) for a in range(N)])
d = gaps(xr); Pz = P(m/d)   # uniform m,s: staggered zone pressure
print("max |F_fd - (P_{i-1/2}-P_{i+1/2})| =", np.abs(F_fd(U_sub,xr) - (np.roll(Pz,1)-Pz)).max())

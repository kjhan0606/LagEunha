import numpy as np, importlib.util, sys
spec = importlib.util.spec_from_file_location("c2", "check2d.py")
src = open("check2d.py").read().split("rng = np.random.default_rng(3)")[0]
exec(src)
rng = np.random.default_rng(3)
x0 = (np.array([[i,j] for i in range(n) for j in range(n)]) + 0.5)*L/n
xr = (x0 + rng.uniform(-0.25,0.25,x0.shape)*L/n) % L
nb0 = neighbors(xr)[0]
def U_at(s):
    xt = xr.copy(); xt[0,0] = (xr[0,0]+s) % L; return Utot(xt)
def nb_at(s):
    xt = xr.copy(); xt[0,0] = (xr[0,0]+s) % L; return neighbors(xt)[0]
lo, hi = 0.080, 0.0805
for _ in range(60):
    mid=0.5*(lo+hi)
    if nb_at(mid)==nb0: lo=mid
    else: hi=mid
sf=0.5*(lo+hi)
for eps in [1e-4,1e-5,1e-6,1e-7]:
    dl=(U_at(sf-eps)-U_at(sf-2*eps))/eps; dr=(U_at(sf+2*eps)-U_at(sf+eps))/eps
    print(f"eps={eps:.0e}: U(sf+eps)-U(sf-eps)={U_at(sf+eps)-U_at(sf-eps):.3e}  dU left={dl:.6f} right={dr:.6f} jump={dr-dl:.4e}")

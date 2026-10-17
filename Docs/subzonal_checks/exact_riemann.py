import numpy as np
from scipy.optimize import brentq

def exact_riemann(rhoL, uL, pL, rhoR, uR, pR, g, xs, t, x0=0.0):
    cL = np.sqrt(g*pL/rhoL); cR = np.sqrt(g*pR/rhoR)
    def fK(p, rhoK, pK, cK):
        if p > pK:
            A = 2/((g+1)*rhoK); B = (g-1)/(g+1)*pK
            return (p-pK)*np.sqrt(A/(p+B))
        return 2*cK/(g-1)*((p/pK)**((g-1)/(2*g))-1)
    f = lambda p: fK(p, rhoL, pL, cL)+fK(p, rhoR, pR, cR)+uR-uL
    ps = brentq(f, 1e-12, 1e6*max(pL, pR))
    us = 0.5*(uL+uR)+0.5*(fK(ps, rhoR, pR, cR)-fK(ps, rhoL, pL, cL))
    rho = np.zeros_like(xs); u = np.zeros_like(xs); p = np.zeros_like(xs)
    xi = (xs-x0)/t
    for k, s in enumerate(xi):
        if s < us:  # left
            if ps > pL:  # shock
                rsL = rhoL*((ps/pL+(g-1)/(g+1))/((g-1)/(g+1)*ps/pL+1))
                SL = uL-cL*np.sqrt((g+1)/(2*g)*ps/pL+(g-1)/(2*g))
                if s < SL: rho[k], u[k], p[k] = rhoL, uL, pL
                else: rho[k], u[k], p[k] = rsL, us, ps
            else:
                rsL = rhoL*(ps/pL)**(1/g); csL = cL*(ps/pL)**((g-1)/(2*g))
                SH = uL-cL; ST = us-csL
                if s < SH: rho[k], u[k], p[k] = rhoL, uL, pL
                elif s > ST: rho[k], u[k], p[k] = rsL, us, ps
                else:
                    u[k] = 2/(g+1)*(cL+(g-1)/2*uL+s)
                    cc = 2/(g+1)*(cL+(g-1)/2*(uL-s))
                    rho[k] = rhoL*(cc/cL)**(2/(g-1)); p[k] = pL*(cc/cL)**(2*g/(g-1))
        else:
            if ps > pR:
                rsR = rhoR*((ps/pR+(g-1)/(g+1))/((g-1)/(g+1)*ps/pR+1))
                SR = uR+cR*np.sqrt((g+1)/(2*g)*ps/pR+(g-1)/(2*g))
                if s > SR: rho[k], u[k], p[k] = rhoR, uR, pR
                else: rho[k], u[k], p[k] = rsR, us, ps
            else:
                rsR = rhoR*(ps/pR)**(1/g); csR = cR*(ps/pR)**((g-1)/(2*g))
                SH = uR+cR; ST = us+csR
                if s > SH: rho[k], u[k], p[k] = rhoR, uR, pR
                elif s < ST: rho[k], u[k], p[k] = rsR, us, ps
                else:
                    u[k] = 2/(g+1)*(-cR+(g-1)/2*uR+s)
                    cc = 2/(g+1)*(cR-(g-1)/2*(uR-s))
                    rho[k] = rhoR*(cc/cR)**(2/(g-1)); p[k] = pR*(cc/cR)**(2*g/(g-1))
    return rho, u, p, ps, us

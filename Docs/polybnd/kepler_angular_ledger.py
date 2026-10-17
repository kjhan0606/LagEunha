"""Fluid angular momentum from mixed moments, including boundary exchange.

These diagnostics do not replace physical field or orbital accuracy tests.
The caller supplies the active-cell mask used by the evolution.
"""
import numpy as np


def components(state,mass,center,active):
    mass=np.asarray(mass,dtype=float);center=np.asarray(center,dtype=float)
    active=np.asarray(active,dtype=bool);state=np.asarray(state,dtype=float)
    if mass.ndim!=1 or state.shape!=(20*len(mass)+4,) or center.shape!=(2,) or active.shape!=mass.shape:
        raise ValueError('incompatible angular diagnostic shapes')
    if np.any(mass<=0) or not all(np.isfinite(x).all() for x in (state,mass,center)):
        raise ValueError('finite state and positive masses required')
    rows=state[:-4].reshape(-1,20);q=rows[:,3:5];p=rows[:,9:11]
    origin=rows[:,12]-rows[:,13]
    centroid=(q[:,0]*p[:,1]-q[:,1]*p[:,0])/mass
    shift=center[0]*p[:,1]-center[1]*p[:,0]
    return dict(total=float(np.sum((origin-shift)[active])),
                orbital_centroid=float(np.sum((centroid-shift)[active])),
                intrinsic=float(np.sum((origin-centroid)[active])))


def balance(initial,final,mass,center,active):
    before=components(initial,mass,center,active);after=components(final,mass,center,active)
    ledger=np.asarray(final)[-4:]-np.asarray(initial)[-4:]
    exchange=float(ledger[3]-center[0]*ledger[1]+center[1]*ledger[0])
    change=after['total']-before['total'];residual=change-exchange
    return dict(initial=before,final=after,angular_change=change,
                hydro_boundary_angular_exchange=exchange,ledger_residual=residual,
                relative_ledger_residual=abs(residual)/abs(before['total']) if before['total'] else None)

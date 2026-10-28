"""Newtonian stress quadrature contributions; not yet coupled to evolution.

B[a,b]=d v_a/d x_b; D[a,b]=integral rho x_a v_b.
Planar inputs use the 3-D Newtonian constitutive law with v_z=0.
No advective, pressure, gravity or heat-conduction terms are included.
"""
import numpy as np

KEPLER_REYNOLDS = np.inf


def viscosity(reynolds, reference_speed, reference_length):
    if not np.isfinite(reference_speed) or reference_speed <= 0:
        raise ValueError('positive finite reference speed required')
    if not np.isfinite(reference_length) or reference_length <= 0:
        raise ValueError('positive finite reference length required')
    if np.isnan(reynolds) or reynolds <= 0:
        raise ValueError('positive Reynolds number required')
    return reference_speed * reference_length / reynolds


def stress(gradient, density, nu, bulk_viscosity=0.):
    B=np.asarray(gradient,dtype=float)
    if B.shape not in ((2,2),(3,3)) or not np.isfinite(B).all():
        raise ValueError('finite 2-D or 3-D velocity gradient required')
    if not all(np.isfinite(x) for x in (density,nu,bulk_viscosity)) or density<=0 or nu<0 or bulk_viscosity<0:
        raise ValueError('positive density and nonnegative viscosities required')
    theta=np.trace(B);mu=density*nu
    return mu*(B+B.T-2*theta*np.eye(len(B))/3)+bulk_viscosity*theta*np.eye(len(B))


def face_rates(point, velocity, normal, tau, area_weight):
    """Outward physical traction contribution, with no mass transport.

    Use the SAME point, velocity and stress and opposite normals on a pair.
    Pi denotes integral rho v outer v, not total-energy density.
    """
    x=np.asarray(point);v=np.asarray(velocity)
    force=area_weight*(np.asarray(tau)@np.asarray(normal))
    return dict(p=force,D=np.outer(x,force),
                Pi=np.outer(v,force)+np.outer(force,v),E=float(v@force))


def volume_rates(gradient, tau, volume_weight):
    """Integration-by-parts volume terms accompanying face_rates."""
    B=np.asarray(gradient);t=np.asarray(tau)
    return dict(p=np.zeros(len(B)),D=-volume_weight*t.T,
                Pi=-volume_weight*(B@t.T+t@B.T),E=0.)

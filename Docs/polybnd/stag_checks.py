"""Auditable geometry/remap primitives; no floors or boundary mass sources.

The 1e-9 relative closure tolerance is inherited from stagP's mass check.
It detects failed partitions, never repairs geometry or changes physical state.
"""
import numpy as np
from shapely.geometry import Polygon
from shapely.strtree import STRtree


def signed_area(p):
    a, b, c = np.asarray(p)
    return 0.5*((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]))


def check_orientation(p, reference):
    area = signed_area(p)
    if not np.isfinite(area) or not np.isfinite(reference) or area*reference <= 0:
        raise RuntimeError('triangle orientation changed or degenerate')
    return area


def overlap_remap(old, new, mass, internal, rtol=1e-9):
    """Remap extensive m,U and check BOTH donor and receiver area closure.

    Inputs contain only the replaced patch. Each mesh must itself be a valid
    nonoverlapping partition (orientation checks occur before this routine).
    Return mass, internal energy, and independent closure residuals.
    """
    old = [Polygon(p) for p in old]
    new = [Polygon(p) for p in new]
    if not old or not new or any(not p.is_valid or p.area <= 0 for p in old+new):
        raise RuntimeError('invalid remap polygon/empty patch')
    ao = np.array([p.area for p in old]); an = np.array([p.area for p in new])
    mass = np.asarray(mass); internal = np.asarray(internal)
    if (not np.isfinite(mass).all() or not np.isfinite(internal).all()
            or (mass <= 0).any() or (internal <= 0).any()):
        raise RuntimeError('nonpositive/nonfinite donor state')
    # Marginal area sums alone do not detect balanced overlaps and holes.
    # Check each partition independently before constructing the transfer.
    for polygons in (old, new):
        partition_tree = STRtree(polygons)
        for j, p in enumerate(polygons):
            for i in partition_tree.query(p):
                if i < j and p.intersection(polygons[i]).area > rtol*min(p.area, polygons[i].area):
                    raise RuntimeError('remap closure: overlapping partition polygons')
    tree = STRtree(old); covered_old = np.zeros(len(old)); covered_new = np.zeros(len(new))
    mn = np.zeros(len(new)); un = np.zeros(len(new))
    for j, p in enumerate(new):
        for i in tree.query(p):
            overlap = p.intersection(old[i]).area
            covered_old[i] += overlap; covered_new[j] += overlap
            mn[j] += overlap*mass[i]/ao[i]; un[j] += overlap*internal[i]/ao[i]
    donor_error = np.max(np.abs(covered_old/ao-1))
    receiver_error = np.max(np.abs(covered_new/an-1))
    mass_error = abs(mn.sum()/mass.sum()-1)
    energy_error = abs(un.sum()/internal.sum()-1)
    residual = np.array([donor_error, receiver_error, mass_error, energy_error])
    if not np.isfinite(residual).all() or residual.max() > rtol:
        raise RuntimeError('remap closure donor/receiver/mass/U: '+repr(residual))
    if (mn <= 0).any() or (un <= 0).any():
        raise RuntimeError('nonpositive receiver state')
    return mn, un, residual


def atomic_substeps(state, dt, advance, max_substeps=16):
    """advance must not mutate its input or external ledgers.

    Its returned records become visible only after the entire macrostep passes.
    A rejected attempt's records are discarded with its tentative state.
    """
    nsub = 1; rejected = []
    while True:
        current = state; records = []
        try:
            for _ in range(nsub):
                current, record = advance(current, dt/nsub)
                records.append(record)
            return current, records, rejected
        except RuntimeError as exc:
            rejected.append((nsub, str(exc)))
            if nsub >= max_substeps:
                raise RuntimeError('macrostep failed after %d substeps: %s' % (nsub, exc)) from exc
            nsub *= 2

"""Radial binning helpers shared by Sedov3D, Noh3D and Evrard3D."""
import numpy as np


def radii(d, center, box=None, periodic=(0, 0, 0)):
    dx = np.stack([d["x"] - center[0], d["y"] - center[1], d["z"] - center[2]], axis=1)
    if box is not None:
        L = np.array([box[1] - box[0], box[3] - box[2], box[5] - box[4]])
        for k in range(3):
            if periodic[k]:
                dx[:, k] -= L[k] * np.round(dx[:, k] / L[k])
    r = np.linalg.norm(dx, axis=1)
    vr = (dx[:, 0] * d["vx"] + dx[:, 1] * d["vy"] + dx[:, 2] * d["vz"]) / np.maximum(r, 1e-300)
    return r, vr, dx


def binned(r, q, w, edges):
    """Weighted mean of q in radial bins; returns centres, mean, count."""
    idx = np.digitize(r, edges) - 1
    nb = len(edges) - 1
    ok = (idx >= 0) & (idx < nb)
    sw = np.bincount(idx[ok], weights=w[ok], minlength=nb)
    sq = np.bincount(idx[ok], weights=(w * q)[ok], minlength=nb)
    cnt = np.bincount(idx[ok], minlength=nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sq / sw
    return 0.5 * (edges[1:] + edges[:-1]), mean, cnt


def peak_radius(rc, q):
    """Radius of the maximum of a binned profile, parabolic refinement."""
    good = np.isfinite(q)
    i = int(np.nanargmax(np.where(good, q, -np.inf)))
    if 0 < i < len(q) - 1 and np.all(np.isfinite(q[i - 1:i + 2])):
        y0, y1, y2 = q[i - 1], q[i], q[i + 1]
        den = y0 - 2 * y1 + y2
        if den < 0:
            dxb = rc[i + 1] - rc[i]
            return rc[i] + 0.5 * (y0 - y2) / den * dxb, y1
    return rc[i], q[i]

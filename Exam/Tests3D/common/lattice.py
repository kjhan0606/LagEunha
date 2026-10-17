"""Particle-position generators for the 3D tests.

cubic_lattice : cell-centred simple cubic lattice (the 2D tests use the
                same cell-centred grid, x=(i+1/2)dx, when no glass file).
glass         : periodic "glass" by short-range repulsive relaxation of a
                random distribution (cKDTree, k nearest neighbours). This is
                the standard particle-code recipe (White 1994 style); it is
                not a GADGET glass file, and quality is reported (min/mean
                nearest-neighbour distance) so the README can say which was
                used.
"""
import numpy as np


def cubic_lattice(nx, ny, nz, box=(0, 1, 0, 1, 0, 1)):
    x0, x1, y0, y1, z0, z1 = box
    dx, dy, dz = (x1 - x0) / nx, (y1 - y0) / ny, (z1 - z0) / nz
    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                          indexing="ij")
    x = x0 + (i.ravel() + 0.5) * dx
    y = y0 + (j.ravel() + 0.5) * dy
    z = z0 + (k.ravel() + 0.5) * dz
    return x, y, z


def glass(n, box=(0, 1, 0, 1, 0, 1), iters=60, k=16, seed=1, verbose=False):
    """Relax n random points in a periodic box to a glass-like state."""
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    L = np.array([box[1] - box[0], box[3] - box[2], box[5] - box[4]])
    lo = np.array([box[0], box[2], box[4]])
    p = rng.random((n, 3)) * L
    h = (L.prod() / n) ** (1.0 / 3.0)       # mean spacing
    for it in range(iters):
        tree = cKDTree(p, boxsize=L)
        d, idx = tree.query(p, k=k + 1)
        d, idx = d[:, 1:], idx[:, 1:]
        dr = p[:, None, :] - p[idx]
        dr -= L * np.round(dr / L)
        w = np.clip(1.5 * h - d, 0, None) / (d + 1e-12 * h)   # linear repulsion
        f = (dr * w[..., None]).sum(axis=1)
        step = 0.15 * h * (1 - it / iters) + 0.02 * h
        fn = np.linalg.norm(f, axis=1) + 1e-30
        disp = f / fn[:, None] * np.minimum(step, fn)[:, None]
        p = np.mod(p + disp, L)
        if verbose and (it % 10 == 0 or it == iters - 1):
            print("  glass iter %3d: dmin/h=%.3f dmean/h=%.3f" %
                  (it, d[:, 0].min() / h, d[:, 0].mean() / h))
    q = glass_quality(p, L)
    p = p + lo
    return p[:, 0], p[:, 1], p[:, 2], q


def glass_quality(p, L):
    from scipy.spatial import cKDTree
    n = len(p)
    h = (np.prod(L) / n) ** (1.0 / 3.0)
    d, _ = cKDTree(p, boxsize=L).query(p, k=2)
    nn = d[:, 1] / h
    return dict(nn_min=float(nn.min()), nn_mean=float(nn.mean()),
                nn_std=float(nn.std()))


def positions(kind, nx, ny, nz, box, seed=1, verbose=False):
    """kind = 'cubic' or 'glass'. Returns x, y, z, info-dict."""
    if kind == "cubic":
        x, y, z = cubic_lattice(nx, ny, nz, box)
        return x, y, z, dict(kind="cubic", nn_min=1.0, nn_mean=1.0, nn_std=0.0)
    if kind == "glass":
        x, y, z, q = glass(nx * ny * nz, box, seed=seed, verbose=verbose)
        q["kind"] = "glass"
        return x, y, z, q
    raise ValueError(kind)

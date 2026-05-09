#!/usr/bin/env python3
"""
Generate a 2D glass particle distribution on a periodic box [0,Lx]×[0,Ly]
using Lloyd's algorithm (iterative centroidal Voronoi).

Usage:  python3 generate_glass_2d.py [nx] [ny] [Lx] [Ly] [n_iter] [output]

Binary output format:
    int32  np
    float64 Lx
    float64 Ly
    float64[np] x positions
    float64[np] y positions

For periodic BC in Voronoi, we use the 3×3 tiling trick: tile particles
into the 8 neighboring copies of the central box, compute Voronoi on the
tiled set, extract cells for the central tile only, then wrap centroids
back into the central box.
"""
import sys, struct, time
import numpy as np
from scipy.spatial import Voronoi


def tile_periodic(pts, Lx, Ly):
    """Tile points into 3×3 periodic grid. Returns (tiled_pts, center_idx_of_each_tiled)."""
    offs = [(dx*Lx, dy*Ly) for dx in (-1, 0, 1) for dy in (-1, 0, 1)]
    tiles = []
    owner = []
    for ox, oy in offs:
        tiles.append(pts + np.array([ox, oy]))
        owner.append(np.arange(len(pts)))
    return np.concatenate(tiles, axis=0), np.concatenate(owner, axis=0)


def lloyd_periodic(pts, Lx, Ly, n_iter=50, damping=1.0):
    """Lloyd's iteration with periodic BC. Returns relaxed positions."""
    pts = pts.copy()
    center_ids_in_tiled = None
    for it in range(n_iter):
        tiled, owner = tile_periodic(pts, Lx, Ly)
        vor = Voronoi(tiled)
        # For each central-box particle, find its Voronoi cell and compute centroid
        if center_ids_in_tiled is None:
            # Central tile is the 5th copy (0-indexed 4), offset [0,0]
            n_per_tile = len(pts)
            center_ids_in_tiled = np.arange(4*n_per_tile, 5*n_per_tile)

        new_pts = np.empty_like(pts)
        # point_region[i] gives the Voronoi region index for point i
        for i_center, tiled_idx in enumerate(center_ids_in_tiled):
            region_idx = vor.point_region[tiled_idx]
            vertex_ids = vor.regions[region_idx]
            if -1 in vertex_ids or len(vertex_ids) < 3:
                # Unbounded region — shouldn't happen with 3×3 tiling, but safeguard
                new_pts[i_center] = pts[i_center]
                continue
            verts = vor.vertices[vertex_ids]  # (n_verts, 2)
            # Polygon centroid using the shoelace formula
            x = verts[:, 0]; y = verts[:, 1]
            x_next = np.roll(x, -1); y_next = np.roll(y, -1)
            cross = x*y_next - x_next*y
            A = 0.5*np.sum(cross)
            if abs(A) < 1e-30:
                new_pts[i_center] = pts[i_center]
                continue
            cx = np.sum((x + x_next)*cross) / (6.0*A)
            cy = np.sum((y + y_next)*cross) / (6.0*A)
            # Damping: mix old + new
            new_pts[i_center, 0] = (1-damping)*pts[i_center, 0] + damping*cx
            new_pts[i_center, 1] = (1-damping)*pts[i_center, 1] + damping*cy

        # Wrap back to [0,Lx)×[0,Ly)
        new_pts[:, 0] = np.mod(new_pts[:, 0], Lx)
        new_pts[:, 1] = np.mod(new_pts[:, 1], Ly)

        # Report
        displ = np.linalg.norm(new_pts - pts, axis=1)
        mean_d = displ.mean(); max_d = displ.max()
        pts = new_pts
        print(f"  iter {it+1:3d}/{n_iter}: mean_disp={mean_d:.5f}  max_disp={max_d:.5f}")
    return pts


def save_binary(pts, Lx, Ly, path):
    n = len(pts)
    with open(path, 'wb') as f:
        f.write(struct.pack('i', n))
        f.write(struct.pack('dd', Lx, Ly))
        f.write(pts[:, 0].astype(np.float64).tobytes())
        f.write(pts[:, 1].astype(np.float64).tobytes())
    print(f"Saved {n} particles to {path} ({8*2*n + 4 + 16} bytes)")


def save_image(pts, Lx, Ly, path):
    """Save a PNG visualization of the glass distribution.
    For large N: 2D histogram density map + scatter overlay of small patch.
    For small N: scatter of all particles."""
    from PIL import Image, ImageDraw
    n = len(pts)
    # Full-box 2D histogram at 1 particle/pixel average
    nbin_x = max(200, int(np.sqrt(n/(Lx*Ly)*Lx*Lx)))
    nbin_y = int(nbin_x * Ly/Lx)
    H, _, _ = np.histogram2d(pts[:,0], pts[:,1], bins=[nbin_x, nbin_y],
                              range=[[0,Lx],[0,Ly]])
    # Density normalized
    H = H.T  # (ny, nx)
    mean_H = n / (nbin_x*nbin_y)
    vmin, vmax = 0, 2*mean_H
    norm = np.clip((H - vmin)/(vmax - vmin + 1e-30), 0, 1)
    # Grayscale
    arr = (norm * 255).astype(np.uint8)
    full_img = Image.fromarray(arr, mode='L').resize((nbin_x*2, nbin_y*2), Image.NEAREST)

    # Scatter plot of small patch (center 10% × 10%)
    PAD_W = 20
    patch_sz = 600
    patch_img = Image.new('RGB', (patch_sz, patch_sz), 'white')
    pd = ImageDraw.Draw(patch_img)
    cx, cy = Lx/2, Ly/2
    half_w = Lx * 0.05
    half_h = Lx * 0.05  # square patch in world units
    xlo, xhi = cx-half_w, cx+half_w
    ylo, yhi = cy-half_h, cy+half_h
    mask = (pts[:,0]>=xlo)&(pts[:,0]<xhi)&(pts[:,1]>=ylo)&(pts[:,1]<yhi)
    sel = pts[mask]
    for x, y in sel:
        px = int((x-xlo)/(xhi-xlo) * patch_sz)
        py = patch_sz - int((y-ylo)/(yhi-ylo) * patch_sz)
        pd.ellipse((px-2, py-2, px+2, py+2), fill='black')
    pd.text((10, 10), f"center patch {2*half_w:.3f}x{2*half_h:.3f} ({len(sel)} pts)", fill='blue')

    # Combine: density map (left) + scatter patch (right)
    W = full_img.width + PAD_W + patch_sz
    H_out = max(full_img.height, patch_sz)
    canvas = Image.new('RGB', (W, H_out), 'white')
    canvas.paste(full_img.convert('RGB'), (0, 0))
    canvas.paste(patch_img, (full_img.width + PAD_W, 0))
    canvas.save(path)
    print(f"Saved glass image to {path} ({W}x{H_out})")


def main():
    nx    = int(sys.argv[1]) if len(sys.argv) > 1 else 512
    ny    = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
    Lx    = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    Ly    = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
    n_iter= int(sys.argv[5]) if len(sys.argv) > 5 else 50
    out   = sys.argv[6] if len(sys.argv) > 6 else f"glass_{nx}x{ny}.bin"

    np_tot = nx*ny
    # Initial: grid + small random jitter (~10% of cell size)
    dx, dy = Lx/nx, Ly/ny
    xs = (np.arange(nx) + 0.5)*dx
    ys = (np.arange(ny) + 0.5)*dy
    xv, yv = np.meshgrid(xs, ys, indexing='xy')
    rng = np.random.default_rng(42)
    # Large initial jitter (50% of cell size) breaks square-lattice symmetry
    # so Lloyd's can converge to hexagonal CVT (the 2D optimum).
    jitter_x = (rng.random(xv.shape) - 0.5)*0.5*dx
    jitter_y = (rng.random(yv.shape) - 0.5)*0.5*dy
    pts = np.column_stack([(xv+jitter_x).ravel(), (yv+jitter_y).ravel()])
    # Wrap
    pts[:, 0] = np.mod(pts[:, 0], Lx)
    pts[:, 1] = np.mod(pts[:, 1], Ly)

    print(f"Generating glass: nx={nx}, ny={ny}, np={np_tot}, box=[0,{Lx}]x[0,{Ly}], n_iter={n_iter}")
    t0 = time.time()
    pts = lloyd_periodic(pts, Lx, Ly, n_iter=n_iter, damping=1.0)
    dt = time.time() - t0
    print(f"Lloyd's iteration done in {dt:.1f}s ({dt/n_iter:.2f}s/iter)")

    # Stats: nearest-neighbor distance distribution (diagnostic)
    tiled, _ = tile_periodic(pts, Lx, Ly)
    from scipy.spatial import cKDTree
    tree = cKDTree(tiled)
    # central tile particles
    n_per_tile = len(pts)
    center_ids = np.arange(4*n_per_tile, 5*n_per_tile)
    # 2nd-nearest (1st is self)
    d, _ = tree.query(tiled[center_ids], k=2)
    nn = d[:, 1]
    target_nn = np.sqrt(Lx*Ly/np_tot)
    print(f"Nearest-neighbor stats: target~{target_nn:.5f}, mean={nn.mean():.5f}, "
          f"std={nn.std():.5f} (std/mean={nn.std()/nn.mean():.3f})")

    save_binary(pts, Lx, Ly, out)
    png_path = out.rsplit('.', 1)[0] + '.png'
    save_image(pts, Lx, Ly, png_path)


if __name__ == '__main__':
    main()

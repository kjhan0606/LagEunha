"""Cross-check of C arc-mode g (GFS_HG_ARC_DUMP) against the python reference (arcgeom.tangent)."""
import sys, numpy as np, arcgeom as g
f = sys.argv[1]; L = open(f).read().split('\n'); i = 0; worst = 0; worstxi = 0; ncell = 0; skipped = 0
while i < len(L):
    if not L[i].startswith('CELL'): i += 1; continue
    t = L[i].split(); R = float(t[3]); n = int(t[5]); xi_c = float(t[7]); i += 1
    E = [L[i+k].split() for k in range(n)]; i += n
    Gc = np.array([[float(x) for x in L[i+k].split()[2:4]] for k in range(n+1)]); i += n+1
    # neighbour offsets for fluid edges; support index q = k+1 for edge k
    D = []; qmap = []
    for k, e in enumerate(E):
        if int(e[2]) == 1: D.append([float(e[3]), float(e[4])]); qmap.append(k+1)
    D = np.array(D)
    pcs = g.clipped_cell(D, R)
    if pcs is None: skipped += 1; continue
    V, c, arcs = g.arc_xi(pcs, R)
    if len(arcs) != 1: skipped += 1; continue
    xi_p = arcs[0][2]
    # python gradient wrt x_q: neighbour j -> d/dD_j ; i -> -sum_j d/dD_j
    Gp = np.zeros((n+1, 2))
    for jj in range(len(D)):
        for comp in (0, 1):
            dD = np.zeros_like(D); dD[jj, comp] = 1
            _, _, _, _, xs = g.tangent(pcs, D, dD, R)
            Gp[qmap[jj], comp] += xs[0][1]; Gp[0, comp] -= xs[0][1]
    sc = np.abs(Gp).max()
    worst = max(worst, np.abs(Gc - Gp).max()/sc); worstxi = max(worstxi, abs(xi_c - xi_p)/max(abs(xi_p), 1e-12)); ncell += 1
print('cells compared %d (skipped %d): max rel |xi_C - xi_py| = %.2e, max |g_C - g_py|/max|g| = %.2e' % (ncell, skipped, worstxi, worst))

import json
import numpy as np
from weighted_mass_patch import CapacityPatch, extended_support
from test_weighted_mass_patch import SQUARE, IRREGULAR

for name, x in [('square', SQUARE), ('irregular', IRREGULAR)]:
    for rho in ((1., 1.), (1., .01), (.01, 1.)):
        patch = CapacityPatch(x, rho); w, iterations, error = patch.solve(x)
        mass, volume, cells = patch.evaluate(x, w)
        support = extended_support(x, rho)
        centers = []
        for c, m in zip(cells, mass):
            first = np.zeros(2)
            for p, r in zip(patch.donors, rho):
                ov = c.intersection(p)
                if ov.area > 0: first += ov.area*r*np.array(ov.centroid.coords[0])
            centers.append(first/m)
        centers = np.array(centers)
        violations = [max(w[j]-w[i]-np.dot(x[j]-xi, x[j]-xi) for j in range(4) if j != i) for i, xi in enumerate(x)]
        print(json.dumps(dict(shape=name, rho=rho, iterations=iterations,
              mass_relative_error=error, weights=w.tolist(), volumes=volume.tolist(),
              nonincident_mass_fraction=float(support['nonincident']/patch.total),
              extended_support_row_error=float(np.max(abs(support['allocation'].sum(1)-patch.fixed))),
              mass_centroid_generator_distance=float(np.max(np.linalg.norm(centers-x, axis=1))),
              site_own_cell_power_violation=violations)), flush=True)

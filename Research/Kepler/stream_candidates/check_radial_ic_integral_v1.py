"""Independent original-density integral checks for the IC boundary formula."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from prepare_radial_resolution_v1 import density, integrate_density, resolve_near_origin_edges

rows = []
for shift in (0., .046876):
    b = .046875
    vertices = np.array([[-b+shift, -b], [b+shift, -b], [b+shift, b], [-b+shift, b]])
    origin = np.repeat([[shift, 0.]], 4, axis=0)
    fan = (np.array([4*b*b]), np.array([np.sqrt(2)*b]), origin, vertices,
           np.roll(vertices, -1, axis=0), np.zeros(4, dtype=int), np.full(4, shift == 0))
    fan = resolve_near_origin_edges(fan)
    mass = integrate_density(32, fan)[0]
    def section(x):
        return quad(lambda y: float(density(np.hypot(x, y))), -b, b,
                    points=[0.], epsabs=1e-17, epsrel=2e-11, limit=300)[0]
    reference, error = quad(section, shift-b, shift+b,
        points=[0.] if shift < b else None, epsabs=1e-17, epsrel=2e-11, limit=300)
    relative = float(abs(mass-reference)/reference)
    rows.append(dict(square_shift=shift, original_nested_quad=reference,
                     boundary_formula=mass, relative_difference=relative,
                     reference_outer_error_estimate=error,
                     pass_1e_minus8=relative <= 1e-8))
print(json.dumps(dict(scope='two original-density cells; no hydrodynamic evolution', cases=rows), indent=2))
if not all(r['pass_1e_minus8'] for r in rows):
    raise SystemExit(2)

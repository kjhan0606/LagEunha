"""Spatial and temporal checks of the hot pressure-response patch."""
import json
from feedback_patch import run, phase_norm

for h in (.2, .1, .05):
    for dt in (.00125, .000625):
        for remap in (False, True):
            base, _, _ = run(h, .278, dt, 1, remap)
            plus, yp, _ = run(h, .278, dt, 1, remap, 1e-6)
            minus, ym, _ = run(h, .278, dt, 1, remap, -1e-6)
            gain = phase_norm(yp-ym)/phase_norm(yp[0]-ym[0])
            print(json.dumps(dict(base, gain_final=float(gain[-1]), gain_max=float(gain.max()),
                                  pair_flips=[plus['flips'], minus['flips']])), flush=True)

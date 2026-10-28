"""ID-matched radial comparison of accepted aligned W=0 t20 and t30 states."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


BASE = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
AUDITS = {20.0: BASE/"fast_to20_a10_v1/audit_t20.json",
          30.0: BASE/"wallselect_to30_a10_v1/audit_t30.json"}
OUTPUT = Path(__file__).with_name("KEPLER_Q2_ALIGNED_RADIAL_T20_T30_2026-10-27.json")
BANDS = ((.5,.65),(.65,.75),(.75,1.),(1.,1.25),(1.25,1.5),(1.5,1.75),(1.75,2.))
GAMMA_MINUS_ONE = float(np.float32(1.666667))-1.


def weighted_abs(value, mass):
    return float(np.sum(np.abs(value)*mass,dtype=np.longdouble)/
                 np.sum(mass,dtype=np.longdouble))


def main():
    rows=[]
    initial_path=None
    initial_hash=None
    for t, audit_path in AUDITS.items():
        audit=json.loads(audit_path.read_text())
        if audit["time"] != t or not all(audit["checks"].values()):
            raise ValueError(f"Not an accepted t={t} endpoint")
        if initial_path is None:
            initial_path=Path(audit["initial"])
            initial_hash=audit["initial_sha256"]
        if Path(audit["initial"]) != initial_path or audit["initial_sha256"] != initial_hash:
            raise ValueError("Not the same initial state")
        if digest(initial_path) != initial_hash:
            raise ValueError("Changed initial or audit")
        endpoint=Path(audit["endpoint"])
        if digest(endpoint) != audit["endpoint_sha256"]:
            raise ValueError("Changed endpoint")
        p0,m0=ordered(initial_path,55351)
        p,m=ordered(endpoint,55351)
        if m0["t"] != 0 or abs(m["t"]-t)>1e-10 or not np.array_equal(p0["mass"],p["mass"]):
            raise ValueError("Wrong time or changed masses")
        x0=p0["x"].astype(float)-3.
        y0=p0["y"].astype(float)-3.
        r0=np.hypot(x0,y0)
        angle=t*(r0*r0+.01175**2)**(-.75)
        xb=3.+x0*np.cos(angle)-y0*np.sin(angle)
        yb=3.+x0*np.sin(angle)+y0*np.cos(angle)
        position=np.hypot(p["x"].astype(float)-xb,p["y"].astype(float)-yb)
        density=p["den"].astype(float)/p0["den"].astype(float)-1.
        entropy=(np.log(p["ie"].astype(float)/p0["ie"].astype(float))+
                 GAMMA_MINUS_ONE*np.log(p["volume"].astype(float)/p0["volume"].astype(float)))
        mass=p0["mass"].astype(float)
        for lo,hi in BANDS:
            mask=(r0>=lo)&(r0<hi)
            rows.append({"time":t,"r_min":lo,"r_max":hi,"sites":int(mask.sum()),
                         "density_abs_mass_weighted":weighted_abs(density[mask],mass[mask]),
                         "orbit_abs_mass_weighted":weighted_abs(position[mask],mass[mask]),
                         "entropy_abs_mass_weighted":weighted_abs(entropy[mask],mass[mask])})
    result={"audits":{str(t):{"path":str(path),"sha256":digest(path)}
                      for t,path in AUDITS.items()},
            "initial":str(initial_path),"initial_sha256":initial_hash,
            "rows":rows,"scope":"Same IDs, initial-radius bands, fixed masses and softened analytic circular orbit as accepted native endpoint audits. Radial breakdown is diagnostic, not a new global physics gate."}
    with OUTPUT.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    for row in rows:
        print(row)


if __name__=="__main__":
    main()

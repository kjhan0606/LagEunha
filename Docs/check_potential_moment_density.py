"""Actual continuous initial moments: positive fit including Epot, no evolution."""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_split_face_modules_20261018'))
from potential_moment_density import PotentialDensity
from run_canonical_av_branches_20261018 import read
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
source=base/'Codex_continuous_kepler_interface_split_20261018/moments_split_q16.npz';data=np.load(source);ids=data['ids'];q=data['moments']
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');x=np.column_stack([ic['x']-3,ic['y']-3])
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000');order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0]
root=base/'Codex_potential_moment_density_20261018';root.mkdir(exist_ok=False)
rows=[];failures=[];coefficients=[]
for i,a in zip(ids,q):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]];vv=ff[:,6:8]
    poly=vv[np.argsort(np.arctan2(vv[:,1],vv[:,0]))]+x[i]
    try:
        d=PotentialDensity(poly,a[0],a[1:3],np.array([[a[3],a[4]],[a[4],a[5]]]),a[16],geom.potential,geom.quadrature,24)
        pts,w=geom.quadrature(poly,48);rho=d.value(pts);dm=w*rho
        assert np.isfinite(rho).all() and np.all(rho>0)
        F=d.polynomial(pts);target=d.weighted_mass@d.polynomial(d.points)
        spatial_error=float(np.max(abs(dm@F-target))/a[0]);potential_error=float(dm@geom.potential(pts)-a[16])
        rows.append(dict(id=int(i),iterations=d.iterations,condition=d.condition,potential_scale=d.potential_scale,
            fit_residual=d.residual,independent_spatial_relative_error=spatial_error,independent_potential_error=potential_error,
            minimum_vertex_density=float(d.value(poly).min())))
        coefficients.append(np.r_[i,d.origin,d.scale,d.beta,d.potential_scale,d.theta,d.logZ])
    except Exception as exc:failures.append(dict(id=int(i),error=repr(exc)))
out=dict(scope=__doc__,source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cells=len(ids),fitted=len(rows),failures=failures,
    max_spatial_relative_error=max(r['independent_spatial_relative_error'] for r in rows),
    max_potential_error=max(abs(r['independent_potential_error']) for r in rows),
    sum_potential_error=sum(r['independent_potential_error'] for r in rows),
    maximum_condition=max(r['condition'] for r in rows),minimum_potential_scale=min(r['potential_scale'] for r in rows),
    minimum_vertex_density=min(r['minimum_vertex_density'] for r in rows),maximum_iterations=max(r['iterations'] for r in rows),
    limits='Initial-state fit only with one added independent potential moment. No guarantee during evolution or of face density accuracy. Orthogonalization changes basis, not target conserved moments; no clipping or thermal adjustment.')
np.savez(root/'coefficients.npz',coefficients=coefficients)
(root/'cells.json').write_text(json.dumps(rows,indent=2));(root/'report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

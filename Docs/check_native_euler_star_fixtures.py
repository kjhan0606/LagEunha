"""Scalar/native exact Euler star fixtures and root error envelopes."""
import argparse,hashlib,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_face_record_native_modules_v2_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_coupled_relaxation_evolution_20261018')]
os.environ['LAGFORCE_NATIVE_EULER_STAR_LIB']=str(modules/'libnative_euler_star.so')
from native_euler_star import star_batch,provenance
from euler_star_pressure import star
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
g=float(np.float32(1.666667));critical=4*np.sqrt(g)/(g-1)
fixtures=[('equal',[1,0,1,1,0,1,g]),('Sod',[1,0,1,.125,0,.1,g]),('strong_shock',[1,10,.01,1,-10,.01,g]),('rarefaction',[1,-1,1,1,1,1,g]),('near_vacuum',[1,-.4995*critical,1,1,.4995*critical,1,g]),('vacuum',[1,-.51*critical,1,1,.51*critical,1,g]),('negative_pressure',[1,0,-1,1,0,1,g]),('gamma_one',[1,0,1,1,0,1,1]),('nonfinite_velocity',[1,np.inf,1,1,0,1,g]),('nonfinite_pressure',[1,0,np.nan,1,0,1,g]),('bracket_overflow',[1,1e200,1,1,-1e200,1,g])]
inputs=np.array([v for _,v in fixtures]);out=star_batch(inputs,return_failures=True);rows=[];passed=True
for (name,a),v in zip(fixtures,out):
 try:
  old=np.array(star(a[0],[a[1],0],a[2],a[3],[a[4],0],a[5],[1,0],a[6]))
  assert v[11]==0
  delta=abs(v[:2]-old);envelope=v[7:9]+8*np.finfo(float).eps*abs(old)
  gate=1e-14+1e-10*abs(old)
  assert np.all(delta<=envelope) and np.all(delta<=gate)
  row=dict(name=name,status='passed',native_status=0,scalar=old.tolist(),native=v[:2].tolist(),difference=delta.tolist(),error_envelope=envelope.tolist(),bracket=v[2:6].tolist(),iterations=int(v[9]))
 except (ValueError,RuntimeError) as exc:
  expected=2 if 'vacuum' in str(exc) else 3 if 'overflow' in str(exc) else 1
  okay=int(v[11])==expected;passed &= okay
  row=dict(name=name,status='passed' if okay else 'failed',native_status=int(v[11]),scalar_failure=str(exc),expected_status=expected)
 except Exception as exc:passed=False;row=dict(name=name,status='failed',error=repr(exc),native=v.tolist())
 rows.append(row)
report=dict(status='passed' if passed else 'fixture_gate_failed',fixtures=rows,provenance=provenance(),
 source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),modules/'native_euler_star.py',modules/'native_euler_star.cpp',modules/'libnative_euler_star.so')},
 limits='Same installed SciPy Brent kernel; wave/bracket roundoff envelopes assume usual libm accuracy and are not interval proofs. Extreme overflow is explicitly rejected, no alternate solver.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));assert passed

"""Explicit arc nodes and Cartesian roundoff envelope for cold tangent flow."""
import argparse,hashlib,json,os,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_face_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ['LAGFORCE_NATIVE_EULER_STAR_LIB']=str(modules/'libnative_euler_star.so')
from curved_disk_cell import DiskCell
from curved_face_nodes import build_geometry,prepare_states,arc_prior_cuts
from face_record_batch import build
from closed_weight_root import solve
from kinetic_factor_rhs import ALPHA
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
R=2.5;rho=.01;cell=DiskCell([[-R,-R],[R,-R],[R,R],[-R,R]],[-1,-2,-3,-4],R,17)
class ConstantDensity:
 origin=np.zeros(2)
 def value_trace_batch(self,points,inward_direction,return_branch=False):
  values=np.full(len(points),rho);branch=np.zeros(len(points),np.int8);return (values,branch) if return_branch else values
mass=np.array([rho*cell.area]);factors=np.array([np.eye(2)*np.sqrt(mass[0])*1e-9]);states=[dict(c=np.zeros(2),v=np.zeros(2),B=np.array([[0.,-.4],[.4,0.]]),P=.001,density=ConstantDensity())];positions=np.array([[.125,-.25]]);gen=np.array([[.1,.05]]);priors=[dict(flat=True,value=rho)];rows={}
for order in (8,16,32):
 geometry=build_geometry([cell],np.array([17]),positions,mass,factors,states,gen,np.zeros(1),priors,arc_order=order)
 assert np.max(abs(np.linalg.norm(geometry['points'],axis=1)-R))<2e-15
 assert np.max(abs(geometry['normals']-geometry['points']/R))<4*np.finfo(float).eps
 cache=prepare_states(geometry,states,factors,mass);records,stats,audit=build(cache,[],states,positions,gen,np.zeros(1,bool),np.zeros((1,2)),np.zeros((1,2)),lambda p:np.zeros(2),5/3,ALPHA)
 z,values,norm,diag=solve(records,1,return_diagnostics=True);assert norm<2e-13 and z[0]==0.
 actual=[];oldbound=[]
 for rec,(wl,wr,wn) in zip(records,values):
  actual.append(abs(rec[2]*(sum(wl)+sum(wr))));oldbound.append(128*np.finfo(float).eps*rec[2]*(sum(abs(wl))+sum(abs(wr))))
 rows[str(order)]=dict(status='passed',nodes=len(records),norm=norm,old_normal_speed_only_bound_violation_max=float(max(a-b for a,b in zip(actual,oldbound))),maximum_wall_flux=max(actual),new_cartesian_bound=diag['history'][-1]['wall_mass_flux_bound_max'],compatibility=diag['history'][-1]['compatibility_residual'],compatibility_bound=diag['history'][-1]['compatibility_bound'])
# True prior intersections of x=0 partition the full wall at pi/2,3pi/2.
cut=arc_prior_cuts(cell.arcs[0],dict(flat=False,normal=[1.,0.],offset=0.));assert np.allclose(cut,[0.,np.pi/2,3*np.pi/2,2*np.pi],rtol=0,atol=1e-15)
report=dict(status='passed',scope=__doc__,results=rows,arc_cut_angles=cut,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'curved_face_nodes.py',modules/'closed_weight_root.py']},limits='Prescribed synthetic sample width1e-9, no production covariance floor/noise. Static arc/cold mirror tests only; full dynamic domain root is unverified.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

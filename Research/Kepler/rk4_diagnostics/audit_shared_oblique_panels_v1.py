"""Exact same-domain candidate panel coverage, linear chart and prefix audit."""
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,add,mul,div,neg,encloses,verify
from audit_shared_oblique_slabs_v1 import coordinates
from audit_shared_regular_q_transfer_v1 import rows,sha
ROOT=Path(__file__).resolve().parent
def main():
 run=ROOT/'shared_oblique_panels_v1_mpi1';c=verify(run)
 for name in ('adoption_input_integrity.txt','adoption_output_integrity.txt'):
  v=(run/name).read_text().splitlines();assert v and all(s.endswith(': OK') for s in v)
 _,frames=wire(run/'unchanged_oblique_frames_rank0.bin');_,original=wire(ROOT/'shared_pending29_rebind_v1_mpi1/oblique_coordinate_frames_rank0.bin');assert frames==original
 _,slabs=wire(run/'unchanged_actual_slabs_rank0.bin');_,original=wire(ROOT/'shared_oblique_slabs_v1_mpi1/actual_oblique_slabs_rank0.bin');assert slabs==original
 fm={struct.unpack_from('<Q',b)[0]:b for b in frames};sm={(struct.unpack_from('<Q',b)[0],struct.unpack_from('<Q',b,32)[0]):b for b in slabs};expected={k for k,b in sm.items() if struct.unpack_from('<Q',b,40)[0]==0};assert len(expected)==40
 hp,panels=wire(run/'candidate_oblique_panels_rank0.bin');ha,axes=wire(run/'candidate_axis_proofs_rank0.bin');hl,linear=wire(run/'same_chart_linear_coordinates_rank0.bin');assert (hp[3],ha[3],hl[3])==(192,656,144) and len(linear)==len(frames)==14
 state=(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes();assert sha(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin')=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
 for f,b in zip(frames,linear):
  source=struct.unpack_from('<Q',f,8)[0];apex,part=struct.unpack_from('<2i',f,48);v=[tuple(F(x) for x in struct.unpack_from('<2d',state,64+1064*source+16*k)) for k in range(3)]
  ax,ay=v[0];bx,by=v[1];cx,cy=v[2];det=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax);g1=((cy-ay)/det,-(cx-ax)/det);g2=(-(by-ay)/det,(bx-ax)/det);gr=[tuple(-g1[k]-g2[k] for k in range(2)),g1,g2];i,j=(apex+1)%3,(apex+2)%3
  B=[F(1,3)]*3 if part else [F(int(k in (apex,(apex+1)%3)),2) for k in range(3)];C=[F(int(k in (apex,(apex+2)%3)),2) for k in range(3)] if part else [F(1,3)]*3;d0=B[i]*C[j]-B[j]*C[i]
  for axis,at in ((0,512),(1,448)):
   D=[ball(f,at+32*k)[0] for k in range(2)];li=sum(gr[i][k]*D[k] for k in range(2));lj=sum(gr[j][k]*D[k] for k in range(2));q=(B[i]*lj-B[j]*li)/d0;s=(li*C[j]-lj*C[i])/d0+q;encloses(ball(b,32*axis),(s,F(0)));encloses(ball(b,64+32*axis),(q,F(0)))
  pivot=struct.unpack_from('<i',b,128)[0];assert pivot==int(abs(quad(b,32))>abs(quad(b,0))) and abs(ball(b,32*pivot)[0])>ball(b,32*pivot)[1]
 ledger={int(r['slot_id']):r for r in rows(run/'oblique_panel_cost_clearance.csv')};pb=defaultdict(list);groups=defaultdict(list)
 for b in panels:
  slot,source,target,chart,index,slab,child=struct.unpack_from('<7Q',b);j,n=struct.unpack_from('<2i',b,56);assert (source,target,chart)==struct.unpack_from('<3Q',fm[slot],8) and index==len(pb[slot]) and (slot,slab) in expected;pb[slot].append(b);groups[slot,slab,child].append(b);assert 0<=j<n and n>=16
  encloses(ball(b,128),(F(j,n),F(0)));encloses(ball(b,160),(F(j+1,n),F(0)));assert ball(b,64)[0]+ball(b,64)[1]<ball(b,96)[0]-ball(b,96)[1]
 covered=set()
 for slot,slab in expected:
  children=sorted((key[2],g) for key,g in groups.items() if key[:2]==(slot,slab));assert [k for k,g in children]==list(range(len(children)))
  if ledger[slot]['plan_complete']!='1':continue
  assert children;parent=sm[slot,slab];assert children[0][1][0][64:96]==parent[80:112] and children[-1][1][0][96:128]==parent[112:144]
  previous=None
  for child,g in children:
   n=struct.unpack_from('<i',g[0],60)[0];assert len(g)==n and [struct.unpack_from('<i',b,56)[0] for b in g]==list(range(n))
   assert all(b[64:128]==g[0][64:128] for b in g)
   if previous is not None:assert previous==g[0][64:96]
   previous=g[0][96:128]
  covered.add((slot,slab))
 ab=defaultdict(list);passes=Counter()
 for b in axes:
  slot,index,axis,passed=struct.unpack_from('<4Q',b);assert index<len(pb[slot]) and axis<2;panel=pb[slot][index];ab[slot].append(b)
  U=(F(1,2),F(1,2) if axis else F(5,4));E=(F(1,2),F(5,4) if axis else F(1,2));u=add(ball(panel,64),mul(add(ball(panel,96),neg(ball(panel,64))),U));e=add(ball(panel,128),mul(add(ball(panel,160),neg(ball(panel,128))),E));encloses(ball(b,32),u);encloses(ball(b,64),e)
  if passed:
   passes[slot]+=1;J=ball(b,224);assert J[0]>J[1]>=0
   q,inv,R,prior=quad(b,352),quad(b,368),quad(b,400),quad(b,416);assert 0<=q<1 and 0<=inv<1 and 0<=R<=prior
   measure=mul(mul(mul(ball(fm[slot],576),ball(b,192)),add(ball(panel,96),neg(ball(panel,64)))),add(ball(panel,160),neg(ball(panel,128))));encloses(ball(b,288),measure)
 for slot,r in ledger.items():
  a=ab[slot];expected_order=[(j,axis) for j in range(len(pb[slot])) for axis in range(2)];assert [(struct.unpack_from('<Q',b,8)[0],struct.unpack_from('<Q',b,16)[0]) for b in a]==expected_order[:len(a)]
  assert len(pb[slot])==int(r['candidate_panels']) and len(a)==int(r['checked_axes'])
  if r['clearance_complete']=='1':assert len(a)==passes[slot]==2*len(pb[slot]) and r['plan_complete']=='1'
  elif a:assert passes[slot]==len(a)-1 and not struct.unpack_from('<Q',a[-1],24)[0]
 assert len(panels)==c['candidate_panels'] and len(axes)==c['new_complex_axis_evaluations'] and sum(passes.values())==c['new_complex_axis_passes']
 assert c['completed_candidate_plans']==sum(r['plan_complete']=='1' for r in ledger.values()) and c['completed_oblique_chart_geometry_density_clearance']==sum(r['clearance_complete']=='1' for r in ledger.values())
 assert c['new_physical_quadrature_queries']==c['new_inverse_solves']==c['new_canonical_root_jobs']==c['normal_field_or_panel_replays']==0 and not c['full_quadrature_allowed']
 print(json.dumps(dict(scope='exact original same-chart linear coefficients, actual shared-event endpoints and complete40-strip/eta coverage where planned; native complex-axis prefix/contraction/J/measure accounting, retaining native analytic field/libm premises',run=str(run),independent_oblique_candidate_plan_accounting_pass=True,covered_actual_nonempty_slabs=len(covered),candidate_panels=len(panels),checked_complex_axes=len(axes),passed_complex_axes=sum(passes.values()),completed_oblique_chart_geometry_density_clearance=c['completed_oblique_chart_geometry_density_clearance'],remaining_incomplete_atlas_slots=29,all12_analytic_remainders_complete=False,full_quadrature_allowed=False),indent=2))
if __name__=='__main__':main()

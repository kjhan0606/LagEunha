"""Original first failures only, reversible observer and exact guard predicates."""
from collections import Counter
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import verify,string,ball
ROOT=Path(__file__).resolve().parent
def main():
 run=ROOT/'shared_oblique_guard_witness_v1_mpi1';c=verify(run);old=ROOT/'shared_oblique_panels_v1_mpi1'
 original=(old/'material_shared_oblique_point.h').read_text();observer=(run/'material_shared_oblique_point.h').read_text();expected=original
 changes=[('need(sf::lower(out.s)>0,"collected oblique regular-Duffy radial pole excluded");','if(!(sf::lower(out.s)>0))throw sd::DomainFailure("collected oblique regular-Duffy radial pole excluded",out.s,Q(0));')]
 for phase,expression in [(70,'os::graph(m,p,f,u0,lo,hi,real)'),(71,'face(p,f,l,u0,S(root.n))'),(72,'face(p,f,l,u,S(root.n))'),(73,'face(p,f,l,u,S(predictor.n,sf::add(predictor.e,prior)),&local_s)'),(74,'sd::evaluate(m,f.original,c.s,c.t,&c.angular)')]:
  changes.append((expression,'[&]{sd::TracePhase observer_phase('+str(phase)+');return '+expression+';}()'))
 for a,b in changes:assert expected.count(a)==1;expected=expected.replace(a,b)
 assert observer==expected
 for a,b in reversed(changes):assert expected.count(b)==1;expected=expected.replace(b,a)
 assert expected==original
 _,previous=wire(old/'candidate_axis_proofs_rank0.bin');failed={struct.unpack_from('<3Q',b):b for b in previous if not struct.unpack_from('<Q',b,24)[0]};assert len(failed)==14
 h,observed=wire(run/'first_guard_arguments_rank0.bin');assert h[3]==656 and len(observed)==14;phases=Counter();reasons=Counter();seen=set()
 for b in observed:
  key=struct.unpack_from('<3Q',b);assert key in failed and key not in seen;seen.add(key);p=failed[key];assert b[32:96]==p[32:96] and not struct.unpack_from('<Q',b,24)[0]
  reason=string(b,456,192);assert reason==string(p,456,192);phase,edge=struct.unpack_from('<2i',b,448);phases[phase]+=1;reasons[reason]+=1
  n,e=ball(b,320);limit=quad(b,432);assert e>=0
  if 'binomial domain' in reason:assert abs(n)+e>=limit and limit>0 and b[320:352]==p[320:352]
  elif reason=='regular-q canonical edge denominator excludes zero':assert n-e<=0 and limit==0 and edge in range(3) and b[320:352]==p[320:352]
  else:assert reason=='collected oblique regular-Duffy radial pole excluded' and n-e<=0 and limit==0
 assert seen==set(failed) and phases=={73:12,72:2}
 assert c['matched_original_reasons']==14 and c['passing_complex_axis_replays']==c['new_canonical_root_jobs']==c['new_inverse_solves']==c['new_physical_quadrature_queries']==0 and not c['full_quadrature_allowed']
 print(json.dumps(dict(scope='exact original14first-failed inputs/reasons, byte-reversible diagnostic-only patch and actual guard inequality/phase accounting; no physical instability inference',run=str(run),independent_oblique_guard_witness_pass=True,first_failures=14,phase_counts=dict(phases),reason_counts=dict(reasons),hessian_neighborhood_failures=12,constant_v_forcing_neighborhood_failures=2,passing_axis_replays=0,existing800passing_axes_preserved=True,all_target_actual_atlas_complete=False,full_quadrature_allowed=False),indent=2))
if __name__=='__main__':main()

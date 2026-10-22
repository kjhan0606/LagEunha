"""Exact saved-event ownership, reference maps, inverse inequalities and frames."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json
import struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import Q,rows,sha
ROOT=Path(__file__).resolve().parent
def ball(b,a):return quad(b,a),quad(b,a+16)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[0],a[1]
def mul(a,b):return a[0]*b[0],abs(a[0])*b[1]+abs(b[0])*a[1]+a[1]*b[1]
def div(a,b):
    assert abs(b[0])>b[1];n=a[0]/b[0]
    return n,(a[1]+abs(n)*b[1])/(abs(b[0])-b[1])
def encloses(outer,inner):
    assert outer[1]>=0 and inner[1]>=0 and abs(outer[0]-inner[0])+inner[1]<=outer[1]
def string(b,a,n):return b[a:a+n].split(b'\0',1)[0].decode()
def verify(run):
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    assert (run/'native_exit_code.txt').read_text().strip()=='0'
    for name in ('previous_source_integrity.txt','previous_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        v=(run/name).read_text().splitlines();assert v and all(x.endswith(': OK') for x in v)
    return json.loads((run/'contract.json').read_text())
def main():
    first=ROOT/'shared_pending29_projection_v1_mpi1';run=ROOT/'shared_pending29_rebind_v1_mpi1';c0=verify(first);c=verify(run)
    pending={int(r['slot_id']):r for r in rows(ROOT/'shared_atlas_adoption_v5/remaining29_atlas_slots.csv')};assert len(pending)==29
    node={};edges={}
    for k in range(4):
        for r in rows(ROOT/f'shared_registry_domains_v1_mpi4/domain_nodes_rank{k}.csv'):
            slot=int(r['slot_id'])
            if slot in pending:node[slot,int(r['node'])]=r
        for r in rows(ROOT/f'shared_registry_domains_v1_mpi4/domain_boundaries_rank{k}.csv'):
            slot=int(r['slot_id'])
            if slot in pending:edges.setdefault(slot,[]).append(r)
    used={(slot,int(r[k])) for slot,ev in edges.items() for r in ev for k in ('first_node','last_node')}
    assert len(used)==109
    _,old=wire(first/'actual_reference_events_rank0.bin');h,new=wire(run/'actual_reference_events_rank0.bin');assert h[3]==640 and len(new)==len(old)==109
    _,oldinv=wire(first/'actual_power_inverse_certificates_rank0.bin');h,inverses=wire(run/'actual_power_inverse_certificates_rank0.bin');assert h[3]==608 and len(oldinv)==8 and len(inverses)==12 and inverses[:8]==oldinv
    state=(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes();assert sha(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin')=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    tri={s:64+s*1064 for s in {int(r['source']) for r in pending.values()}}
    vertices={s:[tuple(F(v) for v in struct.unpack_from('<2d',state,a+16*j)) for j in range(3)] for s,a in tri.items()}
    counts=Counter();seen=set();by_id={};changed=[]
    for before,b in zip(old,new):
        slot,source,target,j,kind,passed,reused,index=struct.unpack_from('<8Q',b);key=slot,j;assert key in used and key not in seen;seen.add(key);by_id[key]=b
        assert (source,target)==(int(pending[slot]['source']),int(pending[slot]['target'])) and passed
        r=node[key];assert string(b,288,160)==r['event_id']
        assert ball(b,64)==(Q(r['x']),Q(r['x_radius'])) and ball(b,96)==(Q(r['y']),Q(r['y_radius']))
        if struct.unpack_from('<Q',before,40)[0]:assert b==before
        else:changed.append(key)
        counts[kind]+=1;a=tri[source];v=vertices[source]
        if kind==0:
            ids=struct.unpack_from('<3i',state,a+968);local=ids.index(int(r['global_vertex']))
            for k in range(2):encloses(ball(b,128+32*k),(v[local][k],F(0)))
        elif kind==1:
            e=next(e for e in edges[slot] if int(e['canonical_edge'])>=0 and j in (int(e['first_node']),int(e['last_node'])))
            local=next(k for k in range(3) if struct.unpack_from('<i',state,a+344+200*k+184)[0]==int(e['canonical_edge']))
            start,end=struct.unpack_from('<2i',state,a+944+8*local)
            p=(Q(e['parameter_lo']),Q(e['lo_radius'])) if j==int(e['first_node']) else (Q(e['parameter_hi']),Q(e['hi_radius']))
            for k in range(2):encloses(ball(b,128+32*k),add((v[start][k],F(0)),mul(p,(v[end][k]-v[start][k],F(0)))))
        else:
            assert kind==2 and index<len(inverses);iv=inverses[index];assert struct.unpack_from('<Q',iv,8)[0]==source
            for k in range(2):assert ball(b,64+32*k)==ball(iv,48+32*k) and ball(b,128+32*k)==ball(iv,112+32*k)
    assert seen==used and len(changed)==4 and counts=={0:35,1:60,2:14},counts
    for iv in inverses:
        slot,source,target,chart,inside=struct.unpack_from('<5Q',iv);assert inside==1
        L,R,res=quad(iv,560),quad(iv,576),quad(iv,592);assert 0<=L<1 and R>0 and res>=0
        for k in range(2):
            pos,pe=ball(iv,176+32*k);pn,pb=ball(iv,48+32*k);assert res>=abs(pos-pn)+pe+pb
            assert ball(iv,112+32*k)==(quad(iv,528+16*k),R)
        assert R>=2*res/(1-L)
        X,Y=quad(iv,528),quad(iv,544);v=vertices[source];ax,ay=v[0];bx,by=v[1];cx,cy=v[2];det=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
        g1=((cy-ay)/det,-(cx-ax)/det);g2=(-(by-ay)/det,(bx-ax)/det);g0=tuple(-g1[k]-g2[k] for k in range(2))
        l1=g1[0]*(X-ax)+g1[1]*(Y-ay);l2=g2[0]*(X-ax)+g2[1]*(Y-ay)
        for k,(ln,g) in enumerate(zip((1-l1-l2,l1,l2),(g0,g1,g2))):
            given=ball(iv,400+32*k);assert given[0]>given[1]>=0;encloses(given,(ln,(abs(g[0])+abs(g[1]))*R))
    rebound=rows(run/'rebound_input_certificates.csv');assert len(rebound)==4
    for r in rebound:
        oldindex,newindex=int(r['original_inverse_index']),int(r['new_inverse_index']);a=oldinv[oldindex];b=inverses[newindex]
        assert a[176:400]==b[176:400] and a[496:576]==b[496:576] # same evaluated field, center, L
        assert a[48:112]!=b[48:112] and r['inside']=='1'
    hf,frames=wire(run/'oblique_coordinate_frames_rank0.bin');_,oldframes=wire(first/'oblique_coordinate_frames_rank0.bin');assert hf[3]==2064 and len(frames)==14 and frames==oldframes
    apices={r['slot_id']:r for r in json.loads((ROOT/'shared_oblique_apex_v1/contract.json').read_text())['chart_results']}
    directions=rows(ROOT/'shared_atlas_repairs_v1_mpi1/oblique_direction_proofs.csv');proj_by_frame={}
    for b in frames:
        slot,source,target,chart,passed,mask=struct.unpack_from('<6Q',b);assert passed and (source,target)==(int(pending[slot]['source']),int(pending[slot]['target']))
        d=[ball(b,448+32*k) for k in range(2)];e=[ball(b,512+32*k) for k in range(2)];assert all(x[1]==0 for x in d) and e==[d[1],neg(d[0])]
        determinant=ball(b,576);encloses(determinant,(d[0][0]**2+d[1][0]**2,F(0)));assert determinant[0]>determinant[1]
        cut=quad(b,1856);assert 0<cut<F(apices[slot]['inner_radial_cut_exact'])
        for r in directions:
            if int(r['slot_id'])!=slot:continue
            assert [x[0] for x in d]==[Q(r['direction_x']),Q(r['direction_y'])]
            competitors=sorted({int(r0['id']) for r0 in rows(ROOT/'curved_map_algebra_v2/advected_sites.csv') if int(r0['id'])!=target});index=competitors.index(int(r['competitor']));assert mask&(1<<index)
            gr=[(Q(r[f'cone_g{k}']),Q(r[f'cone_g{k}_radius'])) for k in ('x','y')]
            dv=add(mul(gr[0],d[0]),mul(gr[1],d[1]));actual=ball(b,1120+32*index);encloses(actual,dv);assert abs(actual[0])>actual[1]
        proj_by_frame[slot,chart]=b
    hp,projections=wire(run/'oblique_actual_event_projections_rank0.bin');assert hp[3]==320 and len(projections)==56
    seenp=set()
    for b in projections:
        slot,chart,j,kind=struct.unpack_from('<4Q',b);frame=proj_by_frame[slot,chart];key=slot,j;assert key in by_id and key not in seenp;seenp.add(key)
        event=by_id[key];assert b[32:96]==event[128:192] and string(b,160,160)==string(event,288,160)
        delta=[add(ball(b,32+32*k),neg(ball(frame,384+32*k))) for k in range(2)]
        for axis,at in ((0,512),(1,448)):
            expected=div(add(mul(ball(frame,at),delta[0]),mul(ball(frame,at+32),delta[1])),ball(frame,576));encloses(ball(b,96+32*axis),expected)
    assert seenp=={key for key in used if key[0] in apices}
    assert c['all29_used_reference_event_packets_complete'] and c['new_field_evaluations']==c['new_inverse_Newton_iterations']==c['new_physical_quadrature_queries']==0 and not c['full_quadrature_allowed']
    result=dict(scope='exact saved used-node ownership/reference coordinates, same-source inverse residual/inclusion inequalities,14 actual transverse frames and56 event projection bounds; retains original native field/L engineering premises',
        first_run=str(first),completion_run=str(run),independent_pending29_reference_event_audit_pass=True,
        used_events=109,direct_reference_vertex_events=counts[0],direct_reference_edge_events=counts[1],power_event_bindings=counts[2],
        reused_byte_exact_event_records=105,new_distinct_input_residual_certificates=4,
        all29_used_reference_event_packets_complete=True,all14_oblique_frames_complete=True,
        actual_event_slabs_and_typed_atlases_still_pending=True,remaining_incomplete_atlas_slots=29,full_quadrature_allowed=False)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()

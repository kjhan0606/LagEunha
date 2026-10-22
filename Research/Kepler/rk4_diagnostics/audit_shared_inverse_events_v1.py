"""Saved inverse/event identity, exact ordering and strip packet accounting."""
from collections import Counter,defaultdict
from pathlib import Path
import json
import struct
import sys
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import Q,rows,sha
ROOT=Path(__file__).resolve().parent
def ball(b,at):return quad(b,at),quad(b,at+16)
def main():
    run=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_inverse_events_v1_mpi1');assert run.parent==ROOT
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    for name in ('original_witness_source_integrity.txt','original_witness_output_integrity.txt','certified_inverse_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    assert (run/'native_exit_code.txt').read_text().strip()=='0';c=json.loads((run/'contract.json').read_text())
    policy=ROOT/'SHARED_ATLAS_BUDGET_POLICY_20261021.json';p=json.loads(policy.read_text());assert c['user_authorized_target_query_limit']==p['new_per_target_query_limit']==13074432
    source={};roots=defaultdict(list);root_windows={}
    for rank in range(4):
        _,v=wire(ROOT/f'shared_atlas_witness_v1_mpi4/failed_charts_rank{rank}.bin')
        for b in v:
            slot,s,t,chart=struct.unpack_from('<4Q',b)
            if s==133:source[slot,chart]=b
        _,v=wire(ROOT/f'shared_atlas_witness_v1_mpi4/chart_events_rank{rank}.bin')
        for b in v:
            slot,s,target,chart,face,side,complete,windows,index,count=struct.unpack_from('<10Q',b)
            if s!=133:continue
            assert complete==1;root_windows[slot,chart,face,side]=windows
            if count:roots[slot,chart].append((face,side,ball(b,80)))
    assert len(source)==8 and len(root_windows)==48
    _,v=wire(run/'bound_chart_events_rank0.bin');events=defaultdict(list)
    for b in v:
        slot,s,target,chart,index,kind=struct.unpack_from('<6Q',b);face,side,node=struct.unpack_from('<3q',b,48)
        assert s==133 and (slot,chart) in source
        e={'index':index,'kind':kind,'face':face,'side':side,'node':node,'free':ball(b,80),'s':ball(b,112),'t':ball(b,144)};events[slot,chart].append(e)
    bindings=0
    for key,v in events.items():
        assert [e['index'] for e in v]==list(range(len(v)))
        assert v[0]['kind']==v[-1]['kind']==0 and v[0]['free']==(0,0) and v[-1]['free']==(1,0)
        for a,b in zip(v,v[1:]):assert a['free'][0]+a['free'][1]<b['free'][0]-b['free'][1]
        actual=[(e['face'],e['side'],e['free']) for e in v if e['kind']==1]
        assert sorted(actual)==sorted(roots[key])
        power=[e for e in v if e['kind']==2];assert len(power)<=1
        for e in power:
            bindings+=1;assert e['node'] in (9,11) and e['free']==e['s']
            assert 0<e['free'][0]-e['free'][1]<e['free'][0]+e['free'][1]<1
            assert 0<e['t'][0]-e['t'][1]<e['t'][0]+e['t'][1]<1
    assert bindings==4
    _,new=wire(run/'new_failed_chart_arrangements_rank0.bin');assert len(new)==8
    ledger={(int(r['slot_id']),int(r['chart'])):r for r in rows(run/'inverse_chart_arrangements.csv')};assert set(ledger)==set(source)
    new_counts={}
    for b in new:
        slot,s,target,chart,complete,empty,whole,mask,strips,windows=struct.unpack_from('<10Q',b);key=slot,chart;old=source[key];r=ledger[key]
        assert complete==1 and empty==whole==0 and r['complete']=='1' and int(r['strips'])==strips
        assert b[80:88]==old[80:88] # selected axis and trim face unchanged
        assert b[96:464]==old[96:464] # reference/radial proof unchanged
        # unresolved/selected graph-face diagnostics are newly populated;
        # their old failure occurred before branch selection was reached.
        assert mask==struct.unpack_from('<Q',old,56)[0]
        assert b[592:1328]==old[592:1328] # original selected derivative intervals
        assert windows==int(r['saved_chart_root_windows'])==sum(w for k,w in root_windows.items() if k[:2]==key)
        assert r['reused_rootset_rows']=='6' and r['reused_inverse_count']=='1'
        new_counts[key]=strips
    _,strips=wire(run/'new_partial_strips_rank0.bin');assert len(strips)==22
    counts=Counter()
    for b in strips:
        slot,s,target,chart,index=struct.unpack_from('<5Q',b);key=slot,chart
        assert s==133 and index==counts[key];counts[key]+=1
        ch,dep,lower,upper=struct.unpack_from('<4i',b,48)
        assert ch==chart and dep==struct.unpack_from('<i',source[key],80)[0] and dep in (0,1) and -1<=lower<23 and -1<=upper<23
        lo,hi=ball(b,64),ball(b,96);assert lo[0]+lo[1]<hi[0]-hi[1]
        ev=events[key];assert any(a['free']==lo and z['free']==hi for a,z in zip(ev,ev[1:]))
    assert all(counts[k]==v for k,v in new_counts.items())
    assert c['completed_chart_arrangements']==8 and c['saved_partial_strips']==22 and c['reused_original_rootset_rows']==48 and c['reused_inverse_bindings']==8
    for key in ('new_inverse_solves','new_scalar_root_jobs','new_canonical_root_jobs','new_physical_quadrature_queries','normal_chart_reconstructions'):assert c[key]==0
    assert c['state_unchanged'] and not c['atlas_complete'] and not c['full_quadrature_allowed']
    print(json.dumps({'scope':'exact saved root/input identity, strict event ordering, unchanged chart/derivative witnesses and strip ABI accounting; relies on retained complete rootsets, gradient and shared inverse proofs plus native branch-selection algorithm','run':str(run),'contract':c,'independent_saved_inverse_event_accounting_pass':True,'all8original_inverse_failed_arrangements_complete':True,'strictly_ordered_event_sequences':8,'interior_shared_power_node_projected_into4_charts':4,'retained_original_rootset_rows':48,'actual_strip_count':22,'remaining_incomplete_atlas_slots':33,'analytic_remainders_complete':False,'full_quadrature_allowed':False},indent=2))
if __name__=='__main__':main()

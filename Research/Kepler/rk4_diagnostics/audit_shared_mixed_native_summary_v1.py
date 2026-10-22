"""Plan/analytic/oracle/target accounting after the native integer point audit.

Uses saved native packets. No fields, roots or physical queries are replayed.
Full source closure and transactional state update remain separate gates.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,verify
from audit_shared_regular_q_transfer_v1 import rows,Q,extended,sha

ROOT=Path(__file__).resolve().parent
NAMES=('mass','px','py','angular','kinetic','physical_U','volume','qx','qy','polar','physical_Udot','potential')

def main():
    run=ROOT/'shared_mixed_native_v2_mpi1';c=verify(run)
    point_run=ROOT/'shared_mixed_native_points_audit_v2'
    native=json.loads((point_run/'contract.json').read_text())
    assert native['independent_integer_source_enclosures_pass'] and native['failed_accounting_slots']==0
    assert (point_run/'native_exit_code.txt').read_text().strip()=='0'
    for name in ('source.sha256','binary.sha256','output.sha256'):
        for line in (point_run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    registry=[r for r in rows(ROOT/'shared_atlas_adoption_v7/panel_evidence_registry.csv') if r['panel_kind']!='retained_group']
    assert len(registry)==90121
    prior=json.loads((ROOT/'shared_native_transfer_adoption_v1/contract.json').read_text())
    remaining=set(prior['remaining_target_ids']);assert len(remaining)==19
    expected=[(i,r) for i,r in enumerate(registry) if int(r['target']) in remaining]
    expected.sort(key=lambda pair:int(pair[1]['slot_id']))  # Stable order within each source overlap.
    index=rows(run/'runtime_panel_index.csv');assert len(index)==len(expected)==84489
    proofs=(ROOT/'shared_saved_all12_v1_mpi1/saved_field_all12_panel_proofs_rank0.bin').read_bytes()
    header=struct.unpack_from('<8Q',proofs);assert header[3:5]==(1664,90121)
    sums=defaultdict(lambda:[F(0) for _ in NAMES]);next_panel=defaultdict(int);planned={}
    for r,(original,old) in zip(index,expected):
        slot=int(r['slot_id']);sequence=int(r['runtime_sequence']);kind=int(r['panel_kind'])
        assert slot==int(old['slot_id']) and sequence==next_panel[slot];next_panel[slot]+=1
        assert int(r['native_panel'])==int(old['native_panel']) and int(r['typed_packet_id'])==int(old['typed_packet_id']) and int(r['typed_record'])==int(old['typed_record']) and int(r['analytic_record'])==original
        proof=memoryview(proofs)[64+1664*original:64+1664*(original+1)]
        assert struct.unpack_from('<4Q',proof)==tuple(int(old[k]) for k in ('slot_id','source','target','native_panel'))
        whole,passed=struct.unpack_from('<2Q',proof,32);assert passed and whole==(kind==0)
        assert (kind==2)==(old['panel_kind']=='oblique')
        for k in range(12):sums[slot][k]+=quad(proof,1280+16*k)
        planned[slot,sequence]=r
    _,packet=wire(run/'native_source_packets_rank0.bin');assert len(packet)==465==len(next_panel)
    result={};total_points=complete=0
    point_rows={int(r['slot_id']):r for r in rows(point_run/'exact_native_prefix_accounting.csv')}
    for p in packet:
        slot,source,target,number,count,done,noracle=struct.unpack_from('<7Q',p)
        assert len(p)==832 and slot in next_panel and target in remaining and number==256*next_panel[slot]
        assert count<=number and (not done or count==number)
        r=point_rows[slot];assert int(r['points'])==count and int(r['complete'])==done and r['exact_source_enclosure_pass']=='1'
        for k in range(12):assert quad(p,448+16*k)>=sums[slot][k]
        result[slot]=p;total_points+=count;complete+=done
    assert total_points==native['accounted_actual_points']==c['actual_saved_queries']
    assert complete==native['completed_overlap_slots']==c['completed_source_overlap_slots']
    sites={(int(r['slot_id']),int(r['runtime_sequence']),int(r['iu']),int(r['iv'])) for r in rows(run/'fixed_oracle_sites.csv')}
    assert len(sites)==c['fixed_oracle_sites']
    _,oracles=wire(run/'native_oracle_comparisons_rank0.bin');seen=set();failed=[];by_slot=defaultdict(int)
    for o in oracles:
        slot,source,target,sequence,iu,iv,point,passed=struct.unpack_from('<8Q',o)
        key=slot,sequence,iu,iv;assert key in sites and key not in seen and (slot,sequence) in planned;seen.add(key)
        assert point==sequence*256+iu*16+iv and point<struct.unpack_from('<Q',result[slot],32)[0]
        assert (source,target)==struct.unpack_from('<2Q',result[slot],8)
        valid=True
        for n,at,bt in ((17,64,608),(12,1152,1536)):
            for j in range(n):
                a,b=ball(o,at+32*j),ball(o,bt+32*j)
                valid=valid and abs(a[0]-b[0])<=a[1]+b[1]
        assert bool(passed)==valid
        if not valid:failed.append(dict(slot=slot,runtime_sequence=sequence,iu=iu,iv=iv))
        by_slot[slot]+=1
    assert len(oracles)==c['saved_oracle_comparisons']
    for slot,p in result.items():assert struct.unpack_from('<Q',p,48)[0]==by_slot[slot]
    criteria={(int(r['target']),r['moment']):Q(r['physical_tau']) for r in rows(ROOT/'shared_normalizers_v1_mpi4/all12_physical_tau.csv')}
    budgets=rows(run/'target_native_error_budget.csv');assert len(budgets)==19*12
    target_pass={};worst=[];incomplete=[]
    for r in budgets:
        target=int(r['target']);k=NAMES.index(r['moment']);assert target in remaining
        source=[p for p in result.values() if struct.unpack_from('<Q',p,16)[0]==target]
        all_complete=all(struct.unpack_from('<Q',p,40)[0] for p in source)
        assert bool(int(r['all_selected_sources_complete']))==all_complete
        passed=bool(int(r['readonly_candidate_pass']));assert target not in target_pass or target_pass[target]==passed;target_pass[target]=passed
        tau=Q(r['ORIGINAL_tau']);assert tau==criteria[target,r['moment']]
        if not all_complete:
            assert not passed and r['total_reported_error']=='unproved'
            if target not in incomplete:incomplete.append(target)
            continue
        n,e=Q(r['value']),Q(r['native_input_cache_geometry_tree_error'])
        exact=sum(ball(p,64+32*k)[0] for p in source);radius=sum(ball(p,64+32*k)[1] for p in source)
        assert abs(n-exact)+radius<=e
        analytic=Q(r['analytic_remainder']);assert analytic>=sum(quad(p,448+16*k) for p in source)
        conversion=Q(r['final_R80_conversion_error']);reported=Q(r['total_reported_error'])
        assert reported>=e+analytic+conversion
        if passed:
            out=(run/f'readonly_target_packet_{target}.bin').read_bytes()
            assert len(out)==400 and struct.unpack_from('<2Q',out)==(1,target)
            value=extended(out[16+16*k:32+16*k]);error=extended(out[208+16*k:224+16*k])
            assert conversion==abs(value-n) and error==reported<tau
            if k in (0,5):assert value>error
        worst.append(dict(target=target,moment=r['moment'],native_error_over_tau=float(e/tau),analytic_over_tau=float(analytic/tau),total_error_over_tau=float(reported/tau)))
    assert sum(target_pass.values())==c['readonly_target_candidates_passing_original_total_budget_and_positivity']
    worst.sort(key=lambda r:r['total_error_over_tau'],reverse=True)
    print(json.dumps(dict(scope='fixed full plan identity including panel kind/typed record, saved integer point/tree audit, same-point oracle17/all12 and actual original19target native/analytic/final conversion budgets; no source closure or physical state commit',
        saved_mixed_runtime_accounting_pass=True,completed_overlap_slots=complete,accounted_actual_queries=total_points,
        remaining_unacquired_queries=21629184-total_points,all_saved_analytic_plan_bindings=84489,
        fixed_oracle_sites=len(sites),saved_oracle_comparisons=len(oracles),failed_oracle_comparisons=failed,
        remaining_targets_passing_original_budget=sum(target_pass.values()),passed_target_ids=sorted(t for t,p in target_pass.items() if p),
        incomplete_targets=sorted(incomplete),all_remaining19_total_budgets_pass=all(target_pass.values()) and len(target_pass)==19,
        worst_target_budgets=worst[:16],native_seconds=c['seconds'],integer_audit_seconds=native['seconds'],
        new_field_or_physical_queries_in_audit=0,completed_five_targets_replayed=False,
        full384source_numerical_closure_complete=False,full24target_transfer_complete=False,
        original_caller_state_committed=False,accepted_RK4_evolution=False,long_Kepler_validated=False),indent=2))

if __name__=='__main__':main()

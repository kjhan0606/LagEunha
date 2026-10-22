"""Independent saved ABI/failure accounting and exact binary input predicates."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_regular_q_transfer_v1 import decode_quad,rows,sha

ROOT=Path(__file__).resolve().parent

def wire(path):
    raw=path.read_bytes();h=struct.unpack_from('<8Q',raw)
    assert h[:3]==(0x4753415749544e31,1,1) and h[6:]==(24,384)
    assert len(raw)==64+h[3]*h[4]
    return h,[raw[64+j*h[3]:64+(j+1)*h[3]] for j in range(h[4])]

def quad(raw,at):return decode_quad(int.from_bytes(raw[at:at+16],'little'))

def main():
    run=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_atlas_witness_v1_mpi4')
    assert run.parent==ROOT
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    for name in ('old_source_integrity.txt','old_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    assert int((run/'native_exit_code.txt').read_text())==0
    c=json.loads((run/'contract.json').read_text())
    old=ROOT/'shared_all_atlas_v2_mpi4'
    failures={int(r['slot_id']):r for k in range(4) for r in rows(old/f'actual_atlas_cost_rank{k}.csv') if r['atlas_complete']=='0'}
    ex={(int(r['slot_id']),int(r['chart'])):r for k in range(4) for r in rows(old/f'actual_chart_status_rank{k}.csv') if int(r['slot_id']) in failures and r['complete']=='0'}
    assert len(failures)==35 and len(ex)==22
    original_prefix={}
    for k in range(4):
        raw=(old/f'typed_new_atlas_rank{k}.bin').read_bytes();h=struct.unpack_from('<8Q',raw)
        for j in range(h[4]):
            b=raw[64+j*704:64+(j+1)*704];slot=struct.unpack_from('<Q',b)[0]
            if slot in failures and failures[slot]['reason']=='fixed combined ordinary/clipped panel budget':
                original_prefix[slot,struct.unpack_from('<Q',b,16)[0]]=b
    sites={int(r['id']):[F(float(r[k])) for k in ('mapped_x','mapped_y')] for r in rows(ROOT/'curved_map_algebra_v2/advected_sites.csv')}
    weights={int(r['id']):F(float(r['weight'])) for r in rows(ROOT/'curved_events_v9/site_weights.csv')}
    # The exact vertex is the original binary64 coefficient, not a rounded decimal.
    state=(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes()
    assert sha(ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin')=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    # checkpoint: 8 uint64, immutable triangles, then (4 uint64, 4 binary64).
    sh=struct.unpack_from('<8Q',state)
    assert sh[6]==32 and len(state)==64+sh[4]*sh[5]+sh[3]*64
    coefficients={}
    for j in range(sh[3]):
        at=64+sh[4]*sh[5]+j*64
        gid=struct.unpack_from('<Q',state,at)[0]
        assert gid not in coefficients
        coefficients[gid]=[F(v) for v in struct.unpack_from('<2d',state,at+32)]
    x=coefficients[118]
    counts=Counter();seen=set();members=[];typed_counts=Counter();cost_seen={};inv=[]
    for rank in range(4):
        for r in rows(run/f'failure_witness_rank{rank}.csv'):
            slot=int(r['slot_id']);assert slot in failures and slot//24%4==rank
            assert r['source']==failures[slot]['source'] and r['target']==failures[slot]['target']
            assert r['reason_reproduced']==r['status_equal']=='1' and r['original_reason']==r['observed_reason']
            key=(slot,r['kind'],int(r['chart']));assert key not in seen;seen.add(key);counts[r['kind']]+=1
            if r['kind']=='failed_chart':assert r['observed_reason']==ex[slot,int(r['chart'])]['reason']
            else:assert r['observed_reason']==failures[slot]['reason']
        for stem in ('axis_trials','chart_events','inverse_iterations','inverse_final','radius_pieces','radius_final','failed_charts','partial_strips','saved_cost_prefix'):
            h,records=wire(run/f'{stem}_rank{rank}.bin');assert h[5]==rank;typed_counts[stem]+=len(records)
            for b in records:
                slot,source,target=struct.unpack_from('<3Q',b)
                if stem=='saved_cost_prefix':
                    slot,source,index,target=struct.unpack_from('<4Q',b);assert b==original_prefix[slot,index];cost_seen[slot]=index+1
                assert slot in failures and source==slot//24 and source%4==rank and target==int(failures[slot]['target'])
                if stem=='failed_charts':
                    chart,complete,empty,whole,active,strips,windows=struct.unpack_from('<7Q',b,24);dep=struct.unpack_from('<i',b,80)[0];e=ex[slot,chart]
                    assert (complete,empty,whole,dep,strips,windows)==tuple(int(e[k]) for k in ('complete','empty','whole','dependent','strips','root_windows'))
                    assert b[-192:].split(b'\0',1)[0].decode()==e['reason']
                if stem=='radius_final':
                    # padded 3Q to offset32, then lo2,hi2,r2,displacement2.
                    # Outward down(0) is a negative minimum subnormal; the
                    # original guard is r2>0, not a nonnegative assertion.
                    r2=quad(b,96)
                    assert (r2<=0)==(failures[slot]['reason']=='actual clipped reference region density-origin separation unproved')
        for r in rows(run/f'origin_membership_rank{rank}.csv'):
            target=int(r['target']);assert int(r['global_vertex'])==118 and int(r['source'])%4==rank
            assert [float(r['chi_x']),float(r['chi_y'])]==[float(a) for a in x]
            own=sites[target];positive=zero=0;owner=target
            for other,b in sites.items():
                if other==target:continue
                n=[b[k]-own[k] for k in range(2)]
                offset=(sum(a*a for a in b)-sum(a*a for a in own)+weights[target]-weights[other])/2
                value=sum(n[k]*x[k] for k in range(2))-offset
                positive+=value>0;zero+=value==0
                if value==0:owner=min(owner,other)
            assert (positive,zero)==(int(r['positive_planes']),int(r['zero_planes']))
            assert r['classification']==('outside' if positive else 'boundary' if zero else 'inside')
            assert r['half_open_owner']==('NOT_MEMBER' if positive else str(owner))
            members.append(r)
        inv.extend(rows(run/f'inverse_final_rank{rank}.csv'))
    assert set(int(k[0]) for k in seen)==set(failures)
    assert counts=={'failed_chart':22,'radius_primitive':15,'saved_cost_guard':2}
    assert cost_seen=={3248:3192,4499:3192} and len(original_prefix)==6384
    assert typed_counts['failed_charts']==22 and typed_counts['inverse_final']==8 and typed_counts['inverse_iterations']==128
    assert len(inv)==8 and all(r['inside']=='0' and float(r['input_x'])==float(r['input_y'])==0 for r in inv)
    assert len({tuple(r[k] for k in ('center_x','center_y','radius','residual','L','lambda0','lambda1','lambda2')) for r in inv})==1
    assert len(members)==19
    selected=[r for r in members if failures[int(r['slot_id'])]['reason']=='actual clipped reference region density-origin separation unproved']
    origin_counts=Counter(r['classification'] for r in selected);assert origin_counts=={'inside':7,'outside':8}
    assert c['witnesses_matching_original']==39 and c['failed_chart_reconstructions']==22 and c['normal_chart_reconstructions']==0 and c['compiled_selected_sources']==18
    assert c['new_canonical_root_jobs']==c['new_physical_quadrature_queries']==0 and c['state_unchanged']
    print(json.dumps({'scope':'independent witness ABI, exact original failure accounting, unchanged raw prefix and binary64 power predicates; no new atlas/field proof','run':str(run),'contract':c,'independent_saved_witness_accounting_pass':True,'typed_record_counts':dict(typed_counts),'radius_failure_exact_origin_memberships':dict(origin_counts),'inverse_same_physical_node':True,'inverse_all16_steps_exhausted':True,'atlas_complete':False,'full_quadrature_allowed':False},indent=2))

if __name__=='__main__':main()

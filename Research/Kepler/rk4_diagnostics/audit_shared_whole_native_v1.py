"""Exact saved native point/source/target/error accounting, no field replay.

Sparse integer superaccumulators sum the binary80 native values/radii exactly.
The native basis/elementary-function forward model remains a stated premise;
384 saved Q128 field/moment oracle comparisons are also checked exactly.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json,struct,sys
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,verify
from audit_shared_regular_q_transfer_v1 import rows,Q,extended

ROOT=Path(__file__).resolve().parent
MASK=(1<<48)-1
LOW49=(1<<49)-1
NAMES=('mass','px','py','angular','kinetic','physical_U','volume','qx','qy','polar','physical_Udot','potential')

def integer_binary80_q(data,offset):
    lo,hi=struct.unpack_from('<QQ',data,offset)
    exponent=(hi>>48)&0x7fff;mantissa=((hi&MASK)<<64)|lo
    assert exponent<0x7fff
    if not exponent:
        assert not mantissa  # Actual acquired values are normal or zero.
        return 0,0
    mantissa|=1<<112
    assert not mantissa&LOW49  # Q128 serialization of an exact native binary80.
    return (-1 if hi>>63 else 1)*(mantissa>>49),exponent-16446

def total(bank):
    if not bank:return F(0)
    minimum=min(bank);value=sum(m<<(unit-minimum) for unit,m in bank.items())
    return F(value*(1<<minimum)) if minimum>=0 else F(value,1<<(-minimum))

def main():
    run=ROOT/'shared_whole_native_v1_mpi1';c=verify(run)
    _,source_packets=wire(run/'native_source_packets_rank0.bin')
    _,oracle=wire(run/'native_oracle_comparisons_rank0.bin')
    jobs=rows(ROOT/'shared_whole_native_inputs_v1/slot_jobs.csv')
    job={int(r['slot_id']):r for r in jobs};assert len(job)==len(source_packets)==64
    source_rows={(int(r['slot_id']),r['moment']):r for r in rows(run/'source_native_integrals.csv')}
    result={};oracles=defaultdict(list)
    for o in oracle:
        slot,source,target,panel,iu,iv,index,passed=struct.unpack_from('<8Q',o)
        assert len(o)==1920 and passed==1 and iu==iv==0
        assert source==int(job[slot]['source']) and target==int(job[slot]['target']) and index==panel*256
        for j in range(17):
            a,b=ball(o,64+32*j),ball(o,608+32*j)
            assert abs(a[0]-b[0])<=a[1]+b[1]
        for j in range(12):
            a,b=ball(o,1152+32*j),ball(o,1536+32*j)
            assert abs(a[0]-b[0])<=a[1]+b[1]
        oracles[slot].append(o)
    assert len(oracle)==384
    for number,p in enumerate(source_packets,1):
        slot,source,target,planned,count,complete,noracle=struct.unpack_from('<7Q',p)
        assert len(p)==832 and slot in job and source==int(job[slot]['source']) and target==int(job[slot]['target'])
        assert planned==count==22528 and complete==1 and noracle==len(oracles[slot])==6
        sums=[defaultdict(int) for _ in range(12)];radii=[defaultdict(int) for _ in range(12)]
        path=run/f'positive_slot{slot}.bin';data=path.read_bytes();assert len(data)==1008*count
        for o in oracles[slot]:
            row=memoryview(data)[1008*struct.unpack_from('<Q',o,48)[0]:][:1008]
            for k in range(12):
                assert o[1152+32*k:1168+32*k]==row[432+16*k:448+16*k]
                assert o[1168+32*k:1184+32*k]==row[624+16*k:640+16*k]
            for j in range(17):assert o[80+32*j:96+32*j]==row[160+16*j:176+16*j]
        for i in range(count):
            at=1008*i;assert struct.unpack_from('<Q',data,at)[0]==source
            for k in range(12):
                n,unit=integer_binary80_q(data,at+432+16*k)
                if n:sums[k][unit]+=n
                er,eu=integer_binary80_q(data,at+624+16*k);assert er>=0
                if er:radii[k][eu]+=er
            # Positive mass and physical heat at every accepted actual node.
            for k in (0,5):
                nm=int.from_bytes(data[at+432+16*k:at+448+16*k],'little')
                er=int.from_bytes(data[at+624+16*k:at+640+16*k],'little')
                assert 0<er<nm<(1<<127)
        for k in range(12):
            n,e=ball(p,64+32*k);exact,er=total(sums[k]),total(radii[k])
            assert abs(n-exact)+er<=e
            r=source_rows[slot,NAMES[k]]
            assert r['complete']=='1' and int(r['points'])==count
            assert Q(r['value'])==n and Q(r['native_input_cache_geometry_tree_error'])==e
            assert Q(r['analytic_remainder'])==quad(p,448+16*k)
        result[slot]=p
        if number%8==0:print(f'exact_native_source_accounting={number}/64',file=sys.stderr,flush=True)
    criteria={(int(r['target']),r['moment']):Q(r['physical_tau']) for r in rows(ROOT/'shared_normalizers_v1_mpi4/all12_physical_tau.csv')}
    budgets=rows(run/'target_native_error_budget.csv');assert len(budgets)==48
    worst=[]
    for r in budgets:
        target=int(r['target']);k=NAMES.index(r['moment']);assert r['all_selected_sources_complete']==r['readonly_candidate_pass']=='1'
        related=[p for p in result.values() if struct.unpack_from('<Q',p,16)[0]==target];assert len(related)==16
        n,e=Q(r['value']),Q(r['native_input_cache_geometry_tree_error'])
        exact=sum(ball(p,64+32*k)[0] for p in related);er=sum(ball(p,64+32*k)[1] for p in related)
        assert abs(n-exact)+er<=e
        analytic=Q(r['analytic_remainder']);assert analytic>=sum(quad(p,448+16*k) for p in related)
        packet=(run/f'readonly_target_packet_{target}.bin').read_bytes();assert len(packet)==400 and struct.unpack_from('<2Q',packet)==(1,target)
        value=extended(packet[16+16*k:32+16*k]);reported=extended(packet[208+16*k:224+16*k])
        conversion=Q(r['final_R80_conversion_error']);assert conversion==abs(value-n)
        assert reported==Q(r['total_reported_error']) and reported>=e+analytic+conversion
        tau=Q(r['ORIGINAL_tau']);assert tau==criteria[target,r['moment']] and reported<tau
        if k in (0,5):assert value>reported
        worst.append(dict(target=target,moment=r['moment'],native_error_over_tau=float(e/tau),analytic_over_tau=float(analytic/tau),total_error_over_tau=float(reported/tau)))
    worst.sort(key=lambda r:r['total_error_over_tau'],reverse=True)
    assert c['completed_sources']==64 and c['actual_saved_queries']==1441792 and c['readonly_target_candidates_passing_original_total_budget_and_positivity']==4
    assert not c['original_caller_state_committed'] and not c['full24target_transfer_complete']
    print(json.dumps(dict(scope='all actual1441792native point values/radii, exact integer source sum/error enclosures,384saved canonical oracle comparisons and48original target post-R80 budgets; original native field/libm/cache model retained as premise, no query replay',
       independent_whole_native_saved_accounting_pass=True,selected_targets=[17894,20502,39878,44945],
       actual_native_queries_accounted=1441792,complete_source_overlap_packets=64,exact_source_sums=768,
       canonical_oracle_comparisons=384,readonly_target_candidates_pass=4,
       all48_original_target_total_error_and_positivity_criteria_pass=True,
       worst_target_budgets=worst[:12],native_seconds=c['seconds'],
       new_field_or_physical_queries_in_audit=0,original68526publication_replayed=False,
       original_caller_state_committed=False,full24target_transfer_complete=False,
       accepted_RK4_evolution=False,long_Kepler_validated=False),indent=2))

if __name__=='__main__':main()

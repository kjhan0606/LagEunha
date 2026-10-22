"""Split original grouped source accounting using stored positive rows only."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import csv,json,struct
from audit_shared_atlas_witness_v1 import quad
from audit_shared_regular_q_transfer_v1 import rows,Q,sha

ROOT=Path(__file__).resolve().parent
NAMES=('mass','px','py','angular','kinetic','physical_U','volume','qx','qy','polar','physical_Udot','potential')

def main():
    old=ROOT/'shared_atlas_v4_mpi4';dest=ROOT/'shared_original_whole4_source_accounting_v1';dest.mkdir()
    for line in (old/'output.sha256').read_text().splitlines():
        digest,path=line.split('  ',1);assert sha(Path(path))==digest
    declared=rows(old/'whole_inside_source_integrals.csv');ids=[int(r['source']) for r in declared]
    assert ids==[240,244,248,252] and all(r['all12_complete']=='1' for r in declared)
    sums=defaultdict(lambda:[F(0) for _ in NAMES]);errors=defaultdict(lambda:[F(0) for _ in NAMES]);analytic=defaultdict(lambda:[F(0) for _ in NAMES]);counts=defaultdict(int);inputs=[old/'output.sha256',old/'whole_inside_source_integrals.csv',old/'PARTIAL_four_source_all12_errors.csv']
    for rank,source in enumerate(ids):
        path=old/f'PARTIAL_whole_inside_queries_rank{rank}.bin';inputs.append(path)
        data=path.read_bytes();assert len(data)%448==0
        for at in range(0,len(data),448):
            r=memoryview(data)[at:at+448];assert struct.unpack_from('<Q',r)[0]==source
            assert quad(r,48)>0
            for k in range(12):
                n,e=quad(r,64+16*k),quad(r,256+16*k);assert e>=0
                if k in (0,5):assert n>e
                sums[source][k]+=n;errors[source][k]+=e
            counts[source]+=1
        assert counts[source]==int(declared[rank]['positive_points'])
        proof=old/f'actual_panel_certificates_rank{rank}.csv';inputs.append(proof)
        records=rows(proof);assert len(records)==int(declared[rank]['panels'])
        for r in records:
            assert int(r['source'])==source and r['pass']=='1'
            for k,name in enumerate(NAMES):analytic[source][k]+=Q(r[name+'_remainder'])
    grouped=rows(old/'PARTIAL_four_source_all12_errors.csv');assert len(grouped)==12
    group_checks=[]
    for k,r in enumerate(grouped):
        assert r['moment']==NAMES[k] and r['only_four_sources_complete']=='1'
        exact=sum(sums[s][k] for s in ids);point_error=sum(errors[s][k] for s in ids)
        n,e=Q(r['partial_value']),Q(r['partial_native_geometry_sum_error'])
        assert abs(n-exact)+point_error<=e
        remainder=Q(r['partial_analytic_remainder']);assert remainder>=sum(analytic[s][k] for s in ids)
        group_checks.append(dict(moment=NAMES[k],original_group_contains_exact_saved_source_sums=True,
            exact_group_center_discrepancy=str(abs(n-exact)),exact_total_point_radius=str(point_error),
            original_group_native_error=str(e),original_group_analytic_remainder=str(remainder)))
    output=[]
    for source in ids:
        for k,name in enumerate(NAMES):
            if k in (0,5):assert sums[source][k]>errors[source][k]+analytic[source][k]
            output.append(dict(slot_id=source*24+15,source=source,target=68526,moment=name,
                saved_positive_points=counts[source],exact_saved_native_center=str(sums[source][k]),
                exact_sum_point_radii=str(errors[source][k]),exact_sum_original_analytic_remainders=str(analytic[source][k]),
                native_point_arithmetic='original_binary64_fields_with_Q128_accumulation',
                original_group_packet_replaced=False,new_field_queries=0))
    # Derive target slot from the canonical original site registry, not an
    # assumed literal position (the original group uses the same24 sites).
    targets=sorted(int(r['id']) for r in rows(ROOT/'curved_map_algebra_v2/advected_sites.csv'))
    target_index=targets.index(68526)
    for r in output:r['slot_id']=24*r['source']+target_index
    with (dest/'original_whole4_source_accounting.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(output[0]),lineterminator='\n');writer.writeheader();writer.writerows(output)
    published=ROOT/'shared_missing_source_v1_mpi1/caller_target_packet.bin'
    assert sha(published)=='946b5639159c57f34d9104a7aac537f6dcfa6a0e232d40057ac5aa465df532f3';inputs.append(published)
    result=dict(scope='four original whole source contributions recovered by exact sums of all45056saved positive point rows and176original panel remainder records; original group and published target unchanged; native field/libm premises retained',
        source_ids=ids,actual_saved_points_accounted=sum(counts.values()),source_moment_records=48,
        all12_original_group_reconciliations_pass=True,individual_source_mass_and_U_positive=True,
        exact_group_reconciliations=group_checks,new_field_or_physical_queries=0,
        original68526publication_preserved=True,full384source_numerical_closure_complete=False,
        full24target_transfer_complete=False,accepted_RK4_evolution=False,long_Kepler_validated=False)
    (dest/'contract.json').write_text(json.dumps(result,indent=2)+'\n')
    (dest/'source.py').write_bytes(Path(__file__).read_bytes())
    (dest/'source.sha256').write_text(sha(dest/'source.py')+'  '+str(dest/'source.py')+'\n')
    (dest/'input.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in inputs))
    (dest/'output.sha256').write_text(''.join(sha(dest/n)+'  '+str(dest/n)+'\n' for n in ('contract.json','original_whole4_source_accounting.csv')))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()

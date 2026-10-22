"""Exact saved failed-input accounting; no native fields or root replays."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json,struct
from audit_shared_atlas_witness_v1 import wire
from audit_shared_pending29_projection_v1 import ball,verify,encloses
from audit_shared_regular_q_transfer_v1 import sha

ROOT=Path(__file__).resolve().parent
def main():
    run=ROOT/'shared_native_reference_inputs_v2_mpi1';c=verify(run)
    state=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin'
    assert sha(state)=='74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005'
    raw=state.read_bytes();h,records=wire(run/'native_barycentric_input_witnesses_rank0.bin')
    assert h[3]==624 and len(records)==c['saved_input_witnesses']==4272
    kinds=Counter();sources=Counter();seen=set();example=None
    from audit_shared_regular_q_transfer_v1 import rows
    targets=sorted(int(r['id']) for r in rows(ROOT/'curved_map_algebra_v2/advected_sites.csv'))
    for r in records:
        slot,source,target,kind,panel,iu,iv,oldmask,newmask=struct.unpack_from('<9Q',r)
        key=slot,kind,panel,iu,iv;assert key not in seen;seen.add(key)
        assert slot==24*source+targets.index(target) and iu<16 and iv<16
        original=[ball(r,144+32*j) for j in range(3)]
        old=[ball(r,240+32*j) for j in range(3)]
        direct=[ball(r,336+32*j) for j in range(3)]
        for a,b in zip(old,original):encloses(a,b)
        for a,b in zip(direct,original):encloses(a,b)
        assert abs(sum(n for n,e in original)-1)<=sum(e for n,e in original)
        for k in range(2):
            vertex=[F(struct.unpack_from('<d',raw,64+1064*source+16*j+8*k)[0]) for j in range(3)]
            center=sum(original[j][0]*vertex[j] for j in range(3))
            radius=sum(original[j][1]*abs(vertex[j]) for j in range(3))
            xn,xe=ball(r,80+32*k);assert abs(center-xn)<=radius+xe
        om=nm=0
        for e in range(3):
            a,b=struct.unpack_from('<2i',raw,64+1064*source+944+8*e)
            for values,offset in ((old,432),(direct,528)):
                expected=(values[a][0]+values[b][0],values[a][1]+values[b][1])
                encloses(ball(r,offset+32*e),expected)
            on,oe=ball(r,432+32*e);nn,ne=ball(r,528+32*e)
            if on<=oe:om|=1<<e
            if nn<=ne:nm|=1<<e
        assert oldmask==om and newmask==nm==0 and om
        kinds[kind]+=1;sources[source]+=1
        if example is None:example=dict(slot=slot,source=source,target=target,panel=panel,iu=iu,iv=iv,
            original_lambda=[dict(center=str(n),radius=str(e)) for n,e in original],old_bad_mask=om)
    assert c['geometry_input_failures']==c['direct_reference_lambda_inputs_with_nonseparated_edge_sum']==0
    assert c['legacy_Cartesian_inverse_inputs_with_nonseparated_edge_sum']==len(records)
    assert c['new_full_field_queries']==c['new_physical_quadrature_queries']==0
    print(json.dumps(dict(scope='exact accounting of4272stored input failures and corrected enclosure/edge sums; sample438120geometryinputs, not all remaining native queries or a field-error proof',
        saved_input_failure_accounting_pass=True,saved_lambda_partition_and_Cartesian_reconstruction_consistent=True,checked_actual_reference_inputs=c['checked_actual_reference_inputs'],
        legacy_nonseparated_inputs=len(records),direct_reference_nonseparated_inputs=0,
        failed_input_kind_counts=dict(kinds),failed_input_source_counts=dict(sources),example=example,
        native_seconds=c['seconds'],new_field_or_physical_queries=0,accepted_RK4_evolution=False,long_Kepler_validated=False),indent=2))

if __name__=='__main__':main()

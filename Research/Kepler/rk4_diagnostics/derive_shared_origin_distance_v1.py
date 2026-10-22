"""Exact saved-data source-vertex-to-target separating distance bounds."""
from fractions import Fraction as F
from pathlib import Path
import csv
import json
import struct
import sys
from audit_shared_regular_q_transfer_v1 import binary_round,scale_two,rows,sha

ROOT=Path(__file__).resolve().parent
def main():
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_origin_distance_v1')
    assert dest.parent==ROOT and not dest.exists();dest.mkdir()
    audit=ROOT/'SHARED_ATLAS_WITNESS_20261021.json'
    assert json.loads(audit.read_text())['independent_saved_witness_accounting_pass']
    state=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin';raw=state.read_bytes();h=struct.unpack_from('<8Q',raw)
    coords={}
    for j in range(h[3]):
        at=64+h[4]*h[5]+j*64;gid=struct.unpack_from('<Q',raw,at)[0]
        coords[gid]=tuple(F(v) for v in struct.unpack_from('<2d',raw,at+32))
    sitespath=ROOT/'curved_map_algebra_v2/advected_sites.csv';weightspath=ROOT/'curved_events_v9/site_weights.csv';boundspath=ROOT/'shared_real_bound_v1_mpi4/source_real_bounds.csv'
    sites={int(r['id']):tuple(F(float(r[k])) for k in ('mapped_x','mapped_y')) for r in rows(sitespath)}
    weights={int(r['id']):F(r['weight']) for r in rows(weightspath)}
    assert all(v==0 for v in weights.values()) # same original input, no rounding ambiguity
    bounds={int(r['source']):r for r in rows(boundspath)}
    failed={int(r['slot_id']):r for k in range(4) for r in rows(ROOT/f'shared_all_atlas_v2_mpi4/actual_atlas_cost_rank{k}.csv') if r['reason']=='actual clipped reference region density-origin separation unproved'}
    memberpaths=[ROOT/f'shared_atlas_witness_v1_mpi4/origin_membership_rank{k}.csv' for k in range(4)]
    results=[]
    for path in memberpaths:
        for r in rows(path):
            slot=int(r['slot_id'])
            if slot not in failed:continue
            source=int(r['source']);target=int(r['target']);x=coords[int(r['global_vertex'])];own=sites[target]
            rounded,_,unit=binary_round(F(bounds[source]['tight_L']),64,-16382)
            L=rounded+scale_two(F(1),unit) # exactly retained_sources nextafterl(+infinity)
            assert 0<L<1
            best=F(0);chosen=-1;besth=F(0);normal_norm=F(0)
            for other,p in sites.items():
                if other==target:continue
                n=tuple(p[k]-own[k] for k in range(2));offset=(sum(v*v for v in p)-sum(v*v for v in own)+weights[target]-weights[other])/2
                value=sum(n[k]*x[k] for k in range(2))-offset
                if value>0:
                    norm=sum(abs(v) for v in n);lower=value/(norm*(1+L))
                    if lower>best:best=lower;chosen=other;besth=value;normal_norm=norm
            assert (best>0)==(r['classification']=='outside')
            results.append({'slot_id':slot,'source':source,'target':target,'origin_membership':r['classification'],'separating_competitor':chosen,'exact_positive_plane_value':str(besth),'exact_normal_one_norm':str(normal_norm),'exact_retained_L':str(L),'exact_reference_radius_lower':str(best),'reference_radius_lower_display':float(best),'actual_origin_duffy_required':r['classification']=='inside','atlas_still_pending':True})
    assert len(results)==15 and sum(r['origin_membership']=='outside' for r in results)==8 and sum(r['origin_membership']=='inside' for r in results)==7
    paths=[Path(__file__).resolve(),audit,state,sitespath,weightspath,boundspath,*memberpaths]
    (dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in paths))
    contract={'scope':'exact saved-data separation, with retained ||Dchi-I||infinity<=L and segment inside convex reference triangle; no field/quadrature/root calls','origin_radius_cases':results,'outside8_positive_reference_distance_proved':True,'inside7_original_density_origin_present':True,'formula':'|X|2 >= |X|infinity >= h_i(chi(0))/(||n_i||1*(1+L)); maximize only among original positive separating planes','new_field_queries':0,'new_physical_quadrature_queries':0,'all15_atlases_still_pending':True}
    (dest/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    out=dest/'contract.json';(dest/'output.sha256').write_text(f'{sha(out)}  {out}\n')
    print(json.dumps({'scope':contract['scope'],'outside8_positive_reference_distance_proved':True,'inside7_actual_origin_present':True,'positive_lower_range':[min(r['reference_radius_lower_display'] for r in results if r['origin_membership']=='outside'),max(r['reference_radius_lower_display'] for r in results if r['origin_membership']=='outside')],'all15_atlases_still_pending':True},indent=2))

if __name__=='__main__':main()

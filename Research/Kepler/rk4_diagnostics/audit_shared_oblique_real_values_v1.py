"""Saved real oblique input/slab/branch/contraction accounting; no queries."""
from collections import defaultdict
from pathlib import Path
import json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball,verify
from audit_shared_regular_q_transfer_v1 import rows

ROOT=Path(__file__).resolve().parent

def main():
    run=ROOT/'shared_oblique_real_values_v1_mpi1';c=verify(run)
    inputs=ROOT/'shared_mixed_runtime_inputs_v1'
    groups={int(r['group']):r for r in rows(inputs/'real_graph_groups.csv') if r['kind']=='oblique_same_frame_real_branches'}
    ledger={int(r['group']):r for r in rows(run/'oblique_real_group_values.csv')}
    _,cache=wire(run/'real_oblique_boundary_values_rank0.bin');bygroup=defaultdict(list)
    for p in cache:
        group,slot,source,target,chart,panel,node,passed=struct.unpack_from('<8Q',p)
        assert len(p)==800 and group in groups and node<16 and passed==1
        g=groups[group];r=ledger[group]
        assert (slot,source,target,panel)==tuple(int(g[k]) for k in ('slot_id','source','target','first_native_panel'))
        assert chart==int(r['chart'])
        u=ball(p,64);assert u[0]-u[1]>0 or u[0]+u[1]<0
        for side,face in enumerate((int(r['selected_lower']),int(r['selected_upper']))):
            at=96+224*side
            if face>=0:
                q,inv,R,prior=quad(p,at+32),quad(p,at+48),quad(p,at+80),quad(p,at+96)
                assert 0<=q<1 and 0<=inv<1 and 0<R<=prior
                assert ball(p,at+192)==(0,R)
            else:
                assert all(quad(p,at+j)==0 for j in (0,16,32,48,64,80,96,112,192,208))
        bygroup[group].append(p)
    assert set(bygroup)==set(groups)==set(ledger) and len(groups)==274
    for group,p in bygroup.items():
        assert len(p)==16 and [struct.unpack_from('<Q',x,48)[0] for x in p]==list(range(16))
        assert ledger[group]['complete']=='1' and ledger[group]['stored_nodes']=='16'
    assert len(cache)==c['stored_real_node_records']==4384 and c['completed_groups']==274
    assert c['new_canonical_root_or_arrangement_jobs']==c['new_full_complex_field_queries']==c['new_physical_quadrature_queries']==0
    print(json.dumps(dict(scope='all274same slab/frame selected real branch/contraction/uncertain correction inputs; native analytic/real-root/Hessian/libm premises retained, no field or physical replay',
       independent_oblique_real_cache_accounting_pass=True,complete_oblique_groups=274,
       saved_real_node_records=4384,selected_implicit_branch_values=c['selected_branch_value_bound'],
       all_original_domain_and_work_guards_preserved=True,
       new_field_or_native_queries=0,mixed_native_transfer_complete=False,
       accepted_RK4_evolution=False,long_Kepler_validated=False),indent=2))

if __name__=='__main__':main()

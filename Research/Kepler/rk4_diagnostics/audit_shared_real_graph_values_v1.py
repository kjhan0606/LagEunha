"""Saved real-node/domain/coordinate accounting; no scalar face/value replay.

Retains the native interval-Newton inclusion algorithm and analytic face/libm
premises. An interval overlap alone is not used as a root-identity proof.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json,struct
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_pending29_projection_v1 import ball
from audit_shared_regular_q_transfer_v1 import rows,sha

ROOT=Path(__file__).resolve().parent

def verify(run):
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    assert (run/'native_exit_code.txt').read_text().strip()=='0'
    for name in ('previous_source_integrity.txt','registry_input_integrity.txt','registry_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        values=(run/name).read_text().splitlines();assert values and all(x.endswith(': OK') for x in values)
    return json.loads((run/'contract.json').read_text())

def main():
    run=ROOT/'shared_real_graph_values_v1_mpi1';c=verify(run)
    inputs=ROOT/'shared_mixed_runtime_inputs_v1'
    groups={int(r['group']):r for r in rows(inputs/'real_graph_groups.csv') if r['kind']=='regular_saved_axis1_domain'}
    packets={int(r['packet_id']):r for r in rows(ROOT/'shared_atlas_adoption_v7/typed_packet_registry.csv')}
    jobs={(int(r['slot_id']),int(r['native_panel'])):r for r in rows(inputs/'panel_jobs.csv') if r['panel_kind']=='regular'}
    data={};header={}
    def record(path,index):
        path=str(path)
        if path not in data:
            raw=Path(path).read_bytes();h=struct.unpack_from('<8Q',raw)
            assert h[0] in (0x4753415749544e31,0x475341544c415331) and len(raw)==64+h[3]*h[4]
            data[path]=raw;header[path]=h
        h=header[path];assert 0<=index<h[4]
        return memoryview(data[path])[64+h[3]*index:64+h[3]*(index+1)]
    _,cache=wire(run/'real_regular_boundary_values_rank0.bin')
    bygroup=defaultdict(list);contained=0;extensions=[]
    for r in cache:
        group,slot,source,target,chart,dependent,node,passed,il,ih=struct.unpack_from('<10Q',r)
        assert len(r)==1456 and group in groups and node<16
        g=groups[group];assert (slot,source,target)==tuple(int(g[k]) for k in ('slot_id','source','target'))
        j=jobs[slot,int(g['first_native_panel'])]
        p=record(packets[int(j['typed_packet_id'])]['path'],int(j['typed_record']))
        saved=record(g['prepared_packet'],int(g['prepared_record']))
        assert dependent==struct.unpack_from('<i',p,36)[0] and chart==struct.unpack_from('<i',p,48)[0]
        free=ball(r,80);oldfree=ball(saved,48)
        assert abs(free[0]-oldfree[0])+free[1]<=oldfree[1]
        bygroup[group].append(r)
        if not passed:continue
        for side,face in enumerate(struct.unpack_from('<2i',p,40)):
            if face<0:
                assert (il if side==0 else ih)==0
                continue
            assert 1<=(il if side==0 else ih)<=256
            derivative=ball(r,1136+32*side)
            assert derivative[0]-derivative[1]>0 or derivative[0]+derivative[1]<0
            if dependent==1:
                new=memoryview(r)[112+480*side:592+480*side]
                old=memoryview(saved)[80+608*side:560+608*side]
                for at in (0,32,96,128,160,192,240,272,304,336,368,400,432):
                    assert new[at:at+32]==old[at:at+32]
                assert new[224]==old[224] and new[464]==old[464]
                value,initial=ball(new,64),ball(old,64)
            else:
                value,initial=ball(r,1072+32*side),ball(saved,1296+128*side)
            # The native inclusion chain, same selected function and initial
            # existence certificate prove identity. This endpoint check only
            # accounts for possible final outward rounding beyond the initial.
            excess=abs(value[0]-initial[0])+value[1]-initial[1]
            if excess<=0:contained+=1
            else:extensions.append(dict(group=group,node=node,side=side,exact_excess=str(excess)))
    ledger={int(r['group']):r for r in rows(run/'real_graph_group_values.csv')}
    assert set(ledger)==set(groups)==set(bygroup)
    complete=0;failed=[]
    for group,g in groups.items():
        r=ledger[group];actual=bygroup[group]
        assert int(r['stored_nodes'])==len(actual)
        assert [struct.unpack_from('<Q',x,48)[0] for x in actual]==list(range(len(actual)))
        if r['complete']=='1':
            assert len(actual)==16 and all(struct.unpack_from('<Q',x,56)[0] for x in actual);complete+=1
        else:
            assert not struct.unpack_from('<Q',actual[-1],56)[0]
            b=actual[-1];failed.append(dict(group=group,slot=int(g['slot_id']),node=len(actual)-1,
                reason=r['reason'],argument_center=str(quad(b,1200)),argument_radius=str(quad(b,1216)),
                argument_limit=str(quad(b,1232)),phase=struct.unpack_from('<i',b,1248)[0],edge=struct.unpack_from('<i',b,1252)[0]))
    assert complete==c['completed_groups'] and len(cache)==c['stored_node_records']
    assert sum(int(r['face_jet_evaluations']) for r in ledger.values())==c['scalar_face_jet_evaluations']
    assert sum(int(r['real_branch_values']) for r in ledger.values())==c['actual_real_branch_values']
    assert c['new_canonical_root_or_arrangement_jobs']==c['new_full_complex_field_queries']==c['new_physical_quadrature_queries']==0
    print(json.dumps(dict(scope='saved actual Gauss free-domain/selected-branch/correction-only/refinement accounting; retains native interval-Newton inclusion and analytic face/libm premises; no root identity inferred from interval overlap alone',
        independent_real_graph_cache_accounting_pass=True,scheduled_regular_groups=len(groups),
        complete_groups=complete,stored_real_node_records=len(cache),
        initial_domains_containing_final_intervals=contained,final_outward_extensions=extensions,
        first_failed_groups=failed,all_regular_groups_complete=complete==len(groups),
        real_branch_values=c['actual_real_branch_values'],face_jet_evaluations=c['scalar_face_jet_evaluations'],
        family_floor_stops=c['parameter_family_floor_stops'],
        new_field_or_physical_queries_in_audit=0,mixed_native_transfer_complete=False,
        accepted_RK4_evolution=False,long_Kepler_validated=False),indent=2))

if __name__=='__main__':main()

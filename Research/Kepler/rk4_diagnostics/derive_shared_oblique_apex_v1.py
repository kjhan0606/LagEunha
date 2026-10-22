"""Exact saved apex membership and conservative inner-Duffy-neighborhood proof."""
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_atlas_witness_v1 import wire
from audit_shared_regular_q_transfer_v1 import Q,rows,sha
ROOT=Path(__file__).resolve().parent
def main():
    dest=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_oblique_apex_v1');assert dest.parent==ROOT and not dest.exists()
    state=ROOT/'shared_trace_rhs_v1_mpi4/shared_state.bin';data=state.read_bytes();h=struct.unpack_from('<8Q',data);assert h[5:7]==(1064,32)
    coefficients={}
    for j in range(h[3]):
        at=64+h[4]*h[5]+j*64;gid=struct.unpack_from('<Q',data,at)[0];coefficients[gid]=tuple(F(v) for v in struct.unpack_from('<2d',data,at+32))
    sitespath=ROOT/'curved_map_algebra_v2/advected_sites.csv';weightspath=ROOT/'curved_events_v9/site_weights.csv';boundspath=ROOT/'shared_registry_v2_mpi4/source_image_registry.csv'
    sites={int(r['id']):tuple(F(float(r[k])) for k in ('mapped_x','mapped_y')) for r in rows(sitespath)}
    weights={int(r['id']):F(r['weight']) for r in rows(weightspath)};assert all(v==0 for v in weights.values())
    bounds={int(r['source']):Q(r['L']) for r in rows(boundspath)}
    results=[];paths=[Path(__file__).resolve(),state,sitespath,weightspath,boundspath]
    for rank in range(4):
        path=ROOT/f'shared_atlas_witness_v1_mpi4/failed_charts_rank{rank}.bin';paths.append(path);_,v=wire(path)
        for b in v:
            slot,source,target,chart=struct.unpack_from('<4Q',b)
            reason=b[-192:].split(b'\0',1)[0].decode()
            if reason!='actual_gradient_cone_needs_oblique_event_projection':continue
            apex=struct.unpack_from('<i',b,96)[0];at=64+source*1064
            coords=struct.unpack_from('<6d',data,at);vertices=[tuple(F(x) for x in coords[2*j:2*j+2]) for j in range(3)]
            gid=struct.unpack_from('<i',data,at+968+apex*4)[0];chi=coefficients[gid];R=max(abs(vertices[j][k]-vertices[apex][k]) for j in range(3) for k in range(2));L=bounds[source];assert R>0 and 0<=L<1
            own=sites[target];faces=[]
            for other,p in sorted(sites.items()):
                if other==target:continue
                n=tuple(p[k]-own[k] for k in range(2));offset=(sum(x*x for x in p)-sum(x*x for x in own)+weights[target]-weights[other])/2
                value=sum(n[k]*chi[k] for k in range(2))-offset;M=sum(abs(x) for x in n)*(1+L)*R
                assert M>0;faces.append((other,value,M))
            positive=[(other,value,M) for other,value,M in faces if value>0];zero=[other for other,value,M in faces if value==0]
            if positive:
                chosen,value,M=max(positive,key=lambda x:x[1]/x[2]);cut=min(F(1,2),value/(2*M));classification='outside';proof='inner_original_apex_Duffy_neighborhood_excluded_by_one_exact_target_face';margin=value-M*cut;assert margin>0
            elif not zero:
                chosen,value,M=min(faces,key=lambda x:-x[1]/x[2]);cut=min(F(1,2),-value/(2*M));classification='inside';proof='inner_original_apex_Duffy_neighborhood_inside_all_original_target_faces';margin=min(-v-m*cut for _,v,m in faces);assert margin>0
            else:
                chosen=-1;cut=F(0);classification='boundary';proof='exact_boundary_apex_requires_shared_blowup_event';margin=F(0)
            results.append({'slot_id':slot,'source':source,'target':target,'chart':chart,'apex':apex,'canonical_vertex':gid,'exact_apex_membership':classification,'positive_plane_count':len(positive),'exact_zero_competitors':zero,'whole_source_L':str(L),'reference_ray_infinity_upper':str(R),'inner_radial_cut_exact':str(cut),'inner_radial_cut_display':float(cut),'inner_proof':proof,'selected_face_competitor':chosen,'strict_inner_proof_margin_exact':str(margin),'outer_oblique_arrangement_still_pending':True})
    assert len(results)==14
    dest.mkdir();(dest/'input.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in paths))
    result={'scope':'exact original saved canonical apex/target predicates plus retained whole-source L bound, no field/root/quadrature replay','oblique_chart_count':14,'proved_strict_inner_apex_neighborhoods':sum(r['inner_radial_cut_display']>0 for r in results),'inside_apex_count':sum(r['exact_apex_membership']=='inside' for r in results),'outside_apex_count':sum(r['exact_apex_membership']=='outside' for r in results),'boundary_apex_count':sum(r['exact_apex_membership']=='boundary' for r in results),'proof_formula':'|h(chi(A+rR))-h(chi(A))| <= ||n||1*(1+L)*R_infinity*r; retain half of exact signed margin','chart_results':results,'all14_actual_atlases_still_pending':True,'new_physical_quadrature_queries':0}
    out=dest/'contract.json';out.write_text(json.dumps(result,indent=2)+'\n');(dest/'output.sha256').write_text(f'{sha(out)}  {out}\n');print(json.dumps({k:result[k] for k in ('scope','oblique_chart_count','proved_strict_inner_apex_neighborhoods','inside_apex_count','outside_apex_count','boundary_apex_count','all14_actual_atlases_still_pending')},indent=2))
if __name__=='__main__':main()

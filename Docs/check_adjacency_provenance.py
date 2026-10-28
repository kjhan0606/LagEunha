"""Postprocess actual stored quadrature/dm and face-inventory reproduction.
No nonlinear density fit, RHS, geometry completion or evolution rerun.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_adjacency_prior_impact_20261018';sys.path[:0]=[str(base/'Codex_compiled_midpoint_modules_20261018'),str(base)]
from contact_prior_density import clip_halfplane,polygon_area,features
from density_quadrature import rule
reports={}
for label in ('stage0','midpoint'):
    folder=root/label;common=np.load(folder/'common_inputs.npz');saved=np.load(base/'Codex_compiled_midpoint_once_20261018'/label/'inputs.npz');priors=json.loads((folder/'priors.json').read_text())['old'];fits=json.loads((folder/'fit_residuals.json').read_text());offset=common['polygon_offsets'];polygons=[common['polygon_vertices'][offset[i]:offset[i+1]] for i in range(len(offset)-1)];raw=common['raw'];gravity=saved['gravity_nodes'];goff=saved['gravity_offsets'];point_exact=True;dm_exact=True;maxdm=0.;maxrel=0.;counts=[]
    for i,(poly,meta,fit) in enumerate(zip(polygons,priors,fits)):
        if meta['flat']:parts=[(poly,meta['value'])]
        else:
            normal=np.array(meta['normal']);b=meta['offset'];parts=[(clip_halfplane(poly,-normal,-b),meta['low']),(clip_halfplane(poly,normal,b),meta['high'])]
        pts=[];weights=[];prior=[]
        for part,rho in parts:
            if polygon_area(part)==0:continue
            p,w=rule(part,12);pts.append(p);weights.append(w);prior.append(np.full(len(w),rho))
        points=np.concatenate(pts);weights=np.concatenate(weights);prior=np.concatenate(prior);original=gravity[goff[i]:goff[i+1]]
        assert original.shape[1]==7 and points.shape==original[:,:2].shape
        same=np.array_equal(points,original[:,:2]);point_exact=point_exact and same
        if not same:raise RuntimeError(('stored volume point/order mismatch',label,i))
        area=float(weights.sum());origin=weights@points/area;scale=np.sqrt(area);theta=np.array(fit['old']['theta']);F=features(points,origin,scale);exponent=F@theta+np.log(prior);shift=float(exponent.max());weighted=weights*np.exp(exponent-shift);Z=float(weighted.sum());dm=raw[i,0]*(weighted/Z);source_dm=original[:,2]
        same=np.array_equal(dm,source_dm);dm_exact=dm_exact and same;maxdm=max(maxdm,float(abs(dm-source_dm).max()));maxrel=max(maxrel,float((abs(dm-source_dm)/source_dm).max()));counts.append(len(points))
    # Regenerate ONLY reference geometry to identify the old face segments.
    from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
    import importlib.util
    spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py');clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
    initial=np.load(base/'Codex_production_moment_geometry_20261018/stage0_inputs.npz');a=common['state'][:-4].reshape(-1,20);ids=common['ids'];regenerated,edges=EmbeddedLaguerreGeometry(initial['xy'],ids,clip.clip)(a[:,:2],a[:,2],float(saved['time']))
    assert np.array_equal(np.concatenate(regenerated),common['polygon_vertices'])
    nodes=saved['node_inputs'];owners=saved['owners'];pairs,inverse=np.unique(owners,axis=0,return_inverse=True);sort=np.argsort(inverse,kind='stable');number=np.bincount(inverse);start=np.r_[0,np.cumsum(number)];edge_by_pair={}
    for e,(i,j,*_) in enumerate(edges):edge_by_pair.setdefault((i,j),[]).append(e)
    assert set(map(tuple,pairs))==set(edge_by_pair),'face owner pair inventory differs'
    face_group=np.full(len(nodes),-1,dtype=np.int64);max_normal_error=0.;outside=0.;ambiguous=[]
    for k,pair in enumerate(pairs):
        loc=sort[start[k]:start[k+1]];choices=edge_by_pair[tuple(pair)]
        if len(choices)==1:face_group[loc]=choices[0]
        else:
            assert pair[1]<0,'multiple internal segments need separate review'
            # Exterior -1 is a BC role, never a unique face key. Match actual
            # geometric segments and use normal only as an independent check.
            cost=[]
            for e in choices:
                i,j,A,B,n,length=edges[e];delta=B-A;t=((nodes[loc,:2]-A)@delta)/(length*length);distance=abs((nodes[loc,:2]-A)@n)+np.maximum(np.maximum(-t,t-1),0)*length;cost.append(distance)
            cost=np.asarray(cost);chosen=cost.argmin(axis=0)
            for q in range(len(loc)):
                tied=np.flatnonzero(cost[:,q]==cost[chosen[q],q])
                if len(tied)>1:
                    # Preserve the ambiguity record; normal picks the stored
                    # segment only when its direction agrees exactly.
                    errors=[np.linalg.norm(nodes[loc[q],2:4]-edges[choices[t]][4]) for t in tied];picked=int(tied[np.argmin(errors)]);ambiguous.append(dict(node=int(loc[q]),candidate_faces=[int(choices[t]) for t in tied],normal_selected_face=int(choices[picked])));chosen[q]=picked
            face_group[loc]=np.array(choices)[chosen]
        assigned=np.array([edges[e][4] for e in face_group[loc]]);normal_error=np.linalg.norm(assigned-nodes[loc,2:4],axis=1);max_normal_error=max(max_normal_error,float(normal_error.max()));outside+=int((normal_error>0).sum())
    assert np.all(face_group>=0) and outside==0
    measured=np.bincount(face_group,weights=nodes[:,4],minlength=len(edges));expected=np.array([e[5] for e in edges]);count=np.bincount(face_group,minlength=len(edges));assert np.all(count>0)
    difference=abs(measured-expected);relative=difference/expected
    # A priori sum-roundoff plus Gauss normalization and subsegment arithmetic
    # diagnostic scale, NOT a proof of upstream clipping geometry.
    eps=np.finfo(float).eps;roundoff=(count+64)*eps*expected
    inventory=np.array([[int(ids[e[0]]),int(ids[e[1]]) if e[1]>=0 else -1,*e[2],*e[3],*e[4],e[5]] for e in edges])
    np.savez(folder/'old_face_node_index_provenance.npz',node_index=np.arange(len(nodes)),face_group_index=face_group,face_inventory=inventory,face_columns=np.array(['owner_global_id','neighbor_global_id_or_mirror','Ax','Ay','Bx','By','nx','ny','length']),source_owner_pairs=owners,source_normals=nodes[:,2:4],source_length_weights=nodes[:,4],per_face_node_count=count,per_face_weight_sum=measured)
    reports[label]=dict(volume_cells=len(polygons),volume_nodes=sum(counts),source_offsets_match=np.array_equal(np.r_[0,np.cumsum(counts)],goff),volume_points_shape_order_bitwise_equal=point_exact,weighted_mass_bitwise_equal=dm_exact,maximum_weighted_mass_absolute_difference=maxdm,maximum_weighted_mass_relative_difference=maxrel,regenerated_polygon_matches_audit_input_bitwise=True,original_vertex_bits_available=False,old_face_segments=len(edges),old_face_nodes=len(nodes),face_inventory_owner_pairs_match=True,all_face_segments_used=True,face_normal_bitwise_equal=outside==0,maximum_face_normal_difference=max_normal_error,maximum_face_weight_length_absolute_difference=float(difference.max()),maximum_face_weight_length_relative_difference=float(relative.max()),maximum_sum_diagnostic_scale_fraction=float(np.divide(difference,roundoff,out=np.zeros_like(difference),where=roundoff>0).max()),exterior_ambiguous_node_matches=ambiguous,trace_index_correction='old_saved_face_quadrature_trace.npz face_ids is the node index, not a physical face ID. old_face_node_index_provenance.npz supplies explicit node_index and actual geometric face_group_index.',scope=__doc__,limits='Stored actual quadrature/dm/face inventory reproduction; original vertex byte identity cannot be proved because no original vertex array was stored. Weight diagnostic scale does not certify upstream geometry.')
    (folder/'provenance_report.json').write_text(json.dumps(reports[label],indent=2))
report=dict(status='completed',stages=reports);(root/'provenance_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

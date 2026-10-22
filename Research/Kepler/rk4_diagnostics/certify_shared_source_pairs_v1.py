"""Candidate exact saved-geometry certificates for all source image pairs.

Uses retained L/displacement and canonical C0 topology proofs. No field,
root, quadrature, trajectory or floating-point precision search is run.
Theory/result review remains a separate adoption gate.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import csv
import json
import struct
import sys
import time
from audit_shared_regular_q_transfer_v1 import Q, rows, sha

ROOT=Path(__file__).resolve().parent


def sub(a,b):return a[0]-b[0],a[1]-b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]


def convex_union(a,b,shared):
    first,last=sorted(shared)
    tip_a=next(v for v in a if v not in shared)
    tip_b=next(v for v in b if v not in shared)
    D=sub(last,first);A=sub(tip_a,first);B=sub(tip_b,first);W=sub(tip_a,tip_b)
    if cross(D,A)*cross(D,B)>=0:return None
    denominator=cross(D,W)
    assert denominator
    parameter=cross(A,W)/denominator
    return parameter if 0<=parameter<=1 else None


def cone_separator(a,b,vertex,La,Lb):
    rays_a=[sub(v,vertex) for v in a if v!=vertex]
    rays_b=[sub(v,vertex) for v in b if v!=vertex]
    # On each L1-normal quadrant the four strict inequalities are an
    # EXACT one-variable linear feasibility problem, not sampled angles.
    for sx,sy in ((1,1),(1,-1),(-1,1),(-1,-1)):
        constraints=[]
        for sign,L,rays in ((1,La,rays_a),(-1,Lb,rays_b)):
            for x,y in rays:
                constraints.append((sign*(sx*x-sy*y),sign*sy*y+L*max(abs(x),abs(y))))
        low,high=F(0),F(1);valid=True
        for slope,offset in constraints:
            if slope==0:
                if offset>=0:valid=False;break
            elif slope>0:high=min(high,-offset/slope)
            else:low=max(low,-offset/slope)
        if not valid or low>=high:continue
        p=(low+high)/2
        margin=min(-(slope*p+offset) for slope,offset in constraints)
        assert margin>0 and 0<=p<=1
        n=(sx*p,sy*(1-p));assert abs(n[0])+abs(n[1])==1
        return n,margin
    return None


def sat(a,b,da,db):
    for triangle in (a,b):
        for j in range(3):
            D=sub(triangle[(j+1)%3],triangle[j])
            for sign in (-1,1):
                n=(sign*D[1],-sign*D[0])
                gap=min(dot(n,v) for v in b)-max(dot(n,v) for v in a)
                radius=sum(abs(n[k])*(da[k]+db[k]) for k in range(2))
                if gap>radius:return n,gap-radius
    return None


def main():
    name=sys.argv[1] if len(sys.argv)>1 else "shared_source_pairs_v1"
    dest=ROOT/name
    assert dest.parent==ROOT and not dest.exists()
    start=time.monotonic()
    packet_path=ROOT/"material_p3_finite_v2/mpi1_cpu/curved_map_control.bin"
    packet=packet_path.read_bytes();assert (len(packet)-160)%384==0
    width=(len(packet)-160)//384
    triangles=[]
    for s in range(384):
        values=struct.unpack_from("<6d",packet,160+s*width+24)
        vertices=tuple((F.from_float(values[2*i]),F.from_float(values[2*i+1])) for i in range(3))
        assert cross(sub(vertices[1],vertices[0]),sub(vertices[2],vertices[0]))>0
        triangles.append(vertices)
    image=rows(ROOT/"shared_registry_v2_mpi4/source_image_registry.csv")
    assert len(image)==384
    boxes=[];L=[];displacement=[]
    for source,row in enumerate(image):
        assert int(row["source"])==source
        boxes.append(tuple((Q(row[k])-Q(row[k+"_radius"]),Q(row[k])+Q(row[k+"_radius"])) for k in ("x","y")))
        L.append(Q(row["L"]));assert 0<=L[-1]<1
        displacement.append(tuple(Q(row["displacement_"+k]) for k in ("x","y")))
    edge=rows(ROOT/"shared_registry_v2_mpi4/canonical_edge_registry.csv")
    adjacent={frozenset((int(r["first_source"]),int(r["second_source"]))):r for r in edge if r["physical_external_boundary"]=="0"}
    assert len(adjacent)==536
    preliminary=json.loads((ROOT/"SHARED_REGISTRY_PREFLIGHT_20261021.json").read_text())
    domains=json.loads((ROOT/"SHARED_REGISTRY_DOMAINS_BRIDGED_20261021.json").read_text())
    assert preliminary["preliminary_registry_audit_pass"] and domains["effective_all9216_domain_arrangements_complete"]
    assert domains["saved_event_source_partition_and_edge_cancellation_complete"]
    proofs=[];counts=Counter();pending=[]
    for i in range(384):
        for j in range(i+1,384):
            kind=None;normal=(F(0),F(0));margin=F(0);parameter=""
            for axis in range(2):
                gap=max(boxes[j][axis][0]-boxes[i][axis][1],boxes[i][axis][0]-boxes[j][axis][1])
                if gap>0:
                    kind="strict_retained_image_box_separation";margin=gap;normal=(F(axis==0),F(axis==1));break
            if kind is None:
                a,b=triangles[i],triangles[j];shared=set(a)&set(b)
                if len(shared)==2:
                    assert frozenset((i,j)) in adjacent
                    e=adjacent[frozenset((i,j))]
                    assert int(e["first_orientation"])+int(e["second_orientation"])==0
                    p=convex_union(a,b,shared)
                    if p is not None:
                        kind="convex_union_shared_C0_max_L_less_than_one";parameter=str(p);margin=1-max(L[i],L[j])
                    else:kind="UNPROVED_nonconvex_adjacent_union"
                elif len(shared)==1:
                    solution=cone_separator(a,b,next(iter(shared)),L[i],L[j])
                    if solution:
                        normal,margin=solution;kind="same_mapped_vertex_strict_image_cone_separation"
                    else:kind="UNPROVED_vertex_image_cones"
                elif len(shared)==0:
                    solution=sat(a,b,displacement[i],displacement[j])
                    if solution:
                        normal,margin=solution;kind="strict_reference_SAT_gap_minus_retained_displacement"
                    else:kind="UNPROVED_nonadjacent_image_separation"
                else:kind="UNPROVED_duplicate_reference_triangle"
            counts[kind]+=1
            proofs.append((i,j,kind,str(normal[0]),str(normal[1]),str(margin),parameter))
            if kind.startswith("UNPROVED"):pending.append({"first_source":i,"second_source":j,"reason":kind})
    assert len(proofs)==73536
    dest.mkdir()
    (dest/"certificate_source.py").write_bytes(Path(__file__).read_bytes())
    with (dest/"source_pair_certificates.csv").open("w") as output:
        writer=csv.writer(output);writer.writerow(("first_source","second_source","kind","normal_x","normal_y","strict_margin","convex_diagonal_parameter"));writer.writerows(proofs)
    report={"scope":"Candidate exact saved-geometry global source-pair proof; retained map bounds and canonical C0 trace proofs assumed; no new physical computation",
            "pairs":73536,"certificate_counts":dict(counts),"unproved_pairs":pending,
            "all_pair_certificate_candidates_complete":not pending,
            "global_source_image_disjointness_adopted":False,"mathematical_theory_and_result_review_pending":True,
            "new_field_or_root_evaluations":0,"new_physical_quadrature_queries":0,"full_quadrature_allowed":False,
            "accepted_RK4_evolution":False,"long_Kepler_validated":False,"seconds":time.monotonic()-start}
    (dest/"contract.json").write_text(json.dumps(report,indent=2)+"\n")
    inputs=[packet_path,ROOT/"shared_trace_rhs_v1_mpi4/shared_state.bin",ROOT/"shared_registry_v2_mpi4/source_image_registry.csv",
            ROOT/"shared_registry_v2_mpi4/canonical_edge_registry.csv",ROOT/"SHARED_REGISTRY_PREFLIGHT_20261021.json",ROOT/"SHARED_REGISTRY_DOMAINS_BRIDGED_20261021.json"]
    with (dest/"input.sha256").open("w") as output:
        for path in inputs:output.write(f"{sha(path)}  {path}\n")
    with (dest/"output.sha256").open("w") as output:
        for filename in ("certificate_source.py","source_pair_certificates.csv","contract.json"):
            path=dest/filename;output.write(f"{sha(path)}  {path}\n")
    (dest/"status.txt").write_text("CANDIDATE_SOURCE_PAIR_CERTIFICATES_SAVED_REVIEW_PENDING\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()

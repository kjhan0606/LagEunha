"""Independent exact checks of saved global source-pair inequalities."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_regular_q_transfer_v1 import Q, rows, sha

ROOT=Path(__file__).resolve().parent


def sub(a,b):return a[0]-b[0],a[1]-b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]


def main():
    name=sys.argv[1] if len(sys.argv)>1 else "shared_source_pairs_v1"
    run=ROOT/name;assert run.parent==ROOT
    c=json.loads((run/"contract.json").read_text())
    for filename in ("input.sha256","output.sha256"):
        for line in (run/filename).read_text().splitlines():
            digest,path=line.split("  ",1);assert sha(Path(path))==digest
    theory=ROOT/"ASTRA_SHARED_REGISTRY_COMPLETION_PLAN_20261021.md"
    assert theory.is_file()
    data=(ROOT/"material_p3_finite_v2/mpi1_cpu/curved_map_control.bin").read_bytes()
    assert (len(data)-160)%384==0;width=(len(data)-160)//384
    triangle=[]
    for source in range(384):
        coords=struct.unpack_from("<6d",data,160+source*width+24)
        points=tuple(tuple(F.from_float(x) for x in coords[2*j:2*j+2]) for j in range(3))
        assert cross(sub(points[1],points[0]),sub(points[2],points[0]))>0
        triangle.append(points)
    images=rows(ROOT/"shared_registry_v2_mpi4/source_image_registry.csv")
    edge={frozenset((int(r["first_source"]),int(r["second_source"]))):r for r in rows(ROOT/"shared_registry_v2_mpi4/canonical_edge_registry.csv") if r["physical_external_boundary"]=="0"}
    previous=json.loads((ROOT/"SHARED_REGISTRY_PREFLIGHT_20261021.json").read_text())
    assert previous["preliminary_registry_audit_pass"]
    pairs=rows(run/"source_pair_certificates.csv")
    expected=((i,j) for i in range(384) for j in range(i+1,384))
    assert len(pairs)==73536
    counts=Counter()
    for row,(i,j) in zip(pairs,expected):
        assert (int(row["first_source"]),int(row["second_source"]))==(i,j)
        a,b=triangle[i],triangle[j];shared=set(a)&set(b)
        La,Lb=Q(images[i]["L"]),Q(images[j]["L"])
        assert 0<=La<1 and 0<=Lb<1
        n=F(row["normal_x"]),F(row["normal_y"])
        margin=F(row["strict_margin"]);assert margin>0
        kind=row["kind"];counts[kind]+=1
        if kind=="strict_retained_image_box_separation":
            assert n in ((1,0),(0,1));axis=0 if n[0] else 1;k=("x","y")[axis]
            ai,ar=Q(images[i][k]),Q(images[i][k+"_radius"])
            bi,br=Q(images[j][k]),Q(images[j][k+"_radius"])
            assert margin==max(bi-br-ai-ar,ai-ar-bi-br)
        elif kind=="convex_union_shared_C0_max_L_less_than_one":
            assert len(shared)==2 and frozenset((i,j)) in edge
            proof=edge[frozenset((i,j))];assert int(proof["first_orientation"])+int(proof["second_orientation"])==0
            first,last=sorted(shared);tip_a=next(v for v in a if v not in shared);tip_b=next(v for v in b if v not in shared)
            D=sub(last,first);A=sub(tip_a,first);B=sub(tip_b,first);W=sub(tip_a,tip_b)
            assert cross(D,A)*cross(D,B)<0
            parameter=F(row["convex_diagonal_parameter"])
            assert 0<=parameter<=1 and cross(A,W)==parameter*cross(D,W)
            assert margin==1-max(La,Lb)
        elif kind=="same_mapped_vertex_strict_image_cone_separation":
            assert len(shared)==1 and abs(n[0])+abs(n[1])==1
            vertex=next(iter(shared));margins=[]
            for ray in (sub(v,vertex) for v in a if v!=vertex):margins.append(-dot(n,ray)-La*max(abs(ray[0]),abs(ray[1])))
            for ray in (sub(v,vertex) for v in b if v!=vertex):margins.append(dot(n,ray)-Lb*max(abs(ray[0]),abs(ray[1])))
            assert all(m>0 for m in margins) and margin==min(margins)
        elif kind=="strict_reference_SAT_gap_minus_retained_displacement":
            assert not shared and n!=(0,0)
            gap=min(dot(n,v) for v in b)-max(dot(n,v) for v in a)
            radius=sum(abs(n[k])*(Q(images[i]["displacement_"+("x","y")[k]])+Q(images[j]["displacement_"+("x","y")[k]])) for k in range(2))
            assert margin==gap-radius
        else:raise AssertionError("Unproved/unrecognized pair certificate: "+kind)
    assert dict(counts)==c["certificate_counts"] and not c["unproved_pairs"]
    assert c["all_pair_certificate_candidates_complete"]
    assert c["new_field_or_root_evaluations"]==c["new_physical_quadrature_queries"]==0
    print(json.dumps({"scope":"All73536 exact saved-image/reference/cone inequalities; relies on retained whole-map derivative/displacement enclosures and exact shared C0 function/basis proofs; not independent field re-evaluation",
                      "run":str(run),"all73536_pairs_independently_accounted":True,"certificate_counts":dict(counts),
                      "mathematical_theory_review":str(theory),"theory_review_SHA256":sha(theory),
                      "global_source_image_interiors_disjoint_complete":True,
                      "shared_edge_vertex_contacts_preserved":True,"all_pair_strict_margins_positive":True,
                      "new_physical_quadrature_queries":0,"actual_all_target_atlas_complete":False,
                      "other23_cold_normalizers_complete":False,"all_target_analytic_certificates_complete":False,
                      "full_quadrature_allowed":False,"accepted_RK4_evolution":False,"long_Kepler_validated":False},indent=2))


if __name__=="__main__":main()

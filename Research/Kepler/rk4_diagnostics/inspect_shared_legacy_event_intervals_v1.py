"""ONLY inspect five existing domain/event bindings; no field/root work."""
from collections import defaultdict
from pathlib import Path
import json
from audit_shared_regular_q_transfer_v1 import Q, rows

ROOT=Path(__file__).resolve().parent
failed={246,253,255,285,287}
own=68526
sites=sorted(int(r["id"]) for r in rows(ROOT/"curved_map_algebra_v2/advected_sites.csv"))
competitors=[s for s in sites if s!=own]
nodes=defaultdict(dict)
for row in rows(ROOT/"shared_boundary_witness_v5_mpi1/actual_target_nodes.csv"):
    source=int(row["source"])
    if source in failed:nodes[source][int(row["node"])]=row
root=defaultdict(list)
for rank in range(4):
    for row in rows(ROOT/f"shared_registry_roots_v1_mpi4/canonical_event_root_packets_rank{rank}.csv"):
        if row["root_index"]!="-1":root[int(row["canonical_edge"]),int(row["first_target"]),int(row["second_target"])].append(row)
out=[]
for row in rows(ROOT/"shared_boundary_witness_v5_mpi1/actual_target_boundaries.csv"):
    source=int(row["source"]);edge=int(row["canonical_edge"])
    if source not in failed or edge<0:continue
    for which,coord,radius in (("first","parameter_lo","parameter_lo_radius"),("last","parameter_hi","parameter_hi_radius")):
        node=nodes[source][int(row[which])]
        if int(node["original_vertex_id"])>=0:continue
        active=[int(i) for i in node["active_plane_indices"].split(";") if i]
        assert len(active)==1
        other=competitors[active[0]];key=edge,min(own,other),max(own,other)
        t,e=Q(row[coord]),Q(row[radius]);candidates=[]
        for r in root[key]:
            n,z=Q(r["t"]),Q(r["t_radius"])
            if max(t-e,n-z)<=min(t+e,n+z):candidates.append(r)
        out.append({"source":source,"node":int(row[which]),"canonical_edge":edge,"competitor":other,
                    "old_t":row[coord],"old_radius":row[radius],"matching_new_roots":candidates,
                    "unique_complete_root_intersection":len(candidates)==1,
                    "bitwise_numeric_interval_equal":len(candidates)==1 and Q(candidates[0]["t"])==t and Q(candidates[0]["t_radius"])==e})
print(json.dumps({"scope":"Five original valid retained domain root intervals versus saved complete canonical root sets; no recomputation",
                  "sources":sorted(failed),"bindings":out,
                  "all_bindings_intersect_one_complete_root":all(r["unique_complete_root_intersection"] for r in out),
                  "nonidentical_saved_intervals":sum(not r["bitwise_numeric_interval_equal"] for r in out)},indent=2))

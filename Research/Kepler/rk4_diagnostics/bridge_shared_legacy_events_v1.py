"""Bind ONLY five valid old domains to complete shared event ordinals.

No geometry, interval, quadrature or state is rewritten. The output is a
separate saved-proof identity relation, not a replacement root certificate.
"""
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys
from audit_shared_regular_q_transfer_v1 import Q, rows, sha

ROOT=Path(__file__).resolve().parent
STATE="74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"
SOURCES={246,253,255,285,287}
OWN=68526


def manifest(path):
    return {name:digest for digest,name in (line.split("  ",1) for line in path.read_text().splitlines())}


def main():
    name=sys.argv[1] if len(sys.argv)>1 else "shared_legacy_events_v1"
    dest=ROOT/name
    assert dest.parent==ROOT and not dest.exists()
    old=ROOT/"shared_target_v3_mpi4"
    known=ROOT/"shared_boundary_witness_v5_mpi1"
    rootdir=ROOT/"shared_registry_roots_v1_mpi4"
    registry=ROOT/"shared_registry_v2_mpi4"
    actual=ROOT/"shared_registry_domains_v1_mpi4"
    partial=json.loads((ROOT/"SHARED_REGISTRY_DOMAINS_PARTIAL_20261021.json").read_text())
    assert partial["packet_and_failure_accounting_pass"]
    assert {r["source"] for r in partial["incomplete_domains"]}==SOURCES
    assert all(r["target"]==OWN and r["reason"]=="used edge crossing binds the actual canonical root index" for r in partial["incomplete_domains"])
    assert json.loads((ROOT/"SHARED_REGISTRY_ROOTS_20261021.json").read_text())["contract"]["all_registered_roots_complete"]
    inputs=[]
    for run in (old,rootdir,actual):
        wire=manifest(run/"retained_input.sha256")
        for relative in ("shared_trace_rhs_v1_mpi4/shared_state.bin","curved_map_algebra_v2/advected_sites.csv","curved_events_v9/site_weights.csv"):
            path=ROOT/relative
            assert wire[str(path)]==sha(path)
            inputs.append(path)
        assert wire[str(ROOT/"shared_trace_rhs_v1_mpi4/shared_state.bin")]==STATE
    nodes={}
    for filename,key in (("actual_target_nodes.csv","node"),("actual_target_boundaries.csv","boundary")):
        previous=[r for r in rows(old/filename) if int(r["source"]) in SOURCES]
        current=[r for r in rows(known/filename) if int(r["source"]) in SOURCES]
        assert previous==current  # Original serialized geometry is retained exactly.
        for run in (old,known):
            output=manifest(run/"output.sha256")
            assert output[str(run/filename)]==sha(run/filename)
            inputs.append(run/filename)
        if key=="node":nodes={(int(r["source"]),int(r[key])):r for r in current}
        else:boundaries=current
    original_source={int(r["source"]):r for r in rows(old/"source_target_contributions.csv")}
    coverage={int(r["source"]):r for r in rows(known/"source_coverage.csv")}
    original_events={(int(r["canonical_edge"]),int(r["competitor"])):r for r in rows(old/"canonical_target_edge_events.csv")}
    for source in SOURCES:
        assert original_source[source]["geometry_complete"]=="1"
        assert coverage[source]["coverage_complete"]==coverage[source]["reused_actual_domain"]=="1"
        assert coverage[source]["classification"]=="2" and coverage[source]["empty"]=="0"
    sites=sorted(int(r["target"]) for r in rows(registry/"target_registry.csv"))
    competitors=[t for t in sites if t!=OWN]
    edge={int(r["edge"]):r for r in rows(registry/"canonical_edge_registry.csv")}
    root=defaultdict(list)
    for rank in range(4):
        path=rootdir/f"canonical_event_root_packets_rank{rank}.csv";inputs.append(path)
        for r in rows(path):
            assert r["complete"]=="1"
            if r["root_index"]!="-1":root[int(r["canonical_edge"]),int(r["first_target"]),int(r["second_target"])].append(r)
    used=defaultdict(list)
    for boundary in boundaries:
        for which in ("first","last"):used[int(boundary["source"]),int(boundary[which])].append((boundary,which))
    bindings=[]
    for (source,node),row in sorted(nodes.items()):
        active=[int(i) for i in row["active_plane_indices"].split(";") if i]
        witness=[]
        if int(row["original_vertex_id"])>=0:event="v;"+row["original_vertex_id"]
        elif (source,node) not in used:event="unused;"+str(node)
        elif len(active)>=2:event="p;"+";".join(str(t) for t in sorted({OWN,*(competitors[i] for i in active)}))
        else:
            assert len(active)==1
            other=competitors[active[0]];events=set()
            for boundary,which in used[source,node]:
                canonical=int(boundary["canonical_edge"])
                if canonical<0:continue
                original=original_events[canonical,other]
                assert original["complete"]=="1" and original["reason"]=="complete_roots_with_exact_endpoint_predicates"
                assert source in (int(original["left_source"]),int(original["right_source"]))
                assert {int(original["left_source"]),int(original["right_source"])}=={int(edge[canonical]["first_source"]),int(edge[canonical]["second_source"])}
                parameter="parameter_lo" if which=="first" else "parameter_hi"
                radius=parameter+"_radius"
                n,e=Q(boundary[parameter]),Q(boundary[radius]);assert 0<n-e<=n+e<1
                key=canonical,min(OWN,other),max(OWN,other)
                matches=[]
                for r in root[key]:
                    a,b=Q(r["t"]),Q(r["t_radius"])
                    if max(n-e,a-b)<=min(n+e,a+b):matches.append(r)
                assert len(matches)==1
                r=matches[0];assert r["endpoint"]=="-1"
                event="e;"+";".join(map(str,(*key,int(r["root_index"]))))
                events.add(event)
                a,b=Q(r["t"]),Q(r["t_radius"])
                witness.append({"old_geometry_run":str(old),"old_source_geometry_complete":True,
                                "old_boundary":int(boundary["boundary"]),"which_endpoint":which,
                                "canonical_edge":canonical,"canonical_first_vertex":int(edge[canonical]["first_vertex"]),
                                "canonical_last_vertex":int(edge[canonical]["last_vertex"]),
                                "parameter_orientation_transform":"identity_original_canonical_t",
                                "old_exact_plane_target":OWN,"old_competitor":other,
                                "new_first_target":key[1],"new_second_target":key[2],
                                "plane_sign_relative_to_new":1 if OWN<other else -1,
                                "old_root_existence_complete_event_record":original,
                                "old_t":boundary[parameter],"old_radius":boundary[radius],
                                "new_complete_root_packet":r,"unique_intersection_count":1,
                                "exact_intersection_lo":str(max(n-e,a-b)),"exact_intersection_hi":str(min(n+e,a+b)),
                                "old_interval_unchanged":True})
            assert len(events)==1 and witness
            event=events.pop()
        bindings.append({"epoch":1,"target":OWN,"source":source,"node":node,
                         "original_state_and_basis_SHA256":STATE,"original_vertex":int(row["original_vertex_id"]),
                         "active_competitors":[competitors[i] for i in active],"event_id":event,"root_identity_witnesses":witness})
    # Multiple references to one original event must be the same old event,
    # never a merger of two previously distinct source-boundary roots.
    old_events=defaultdict(set)
    for binding in bindings:
        for w in binding["root_identity_witnesses"]:
            old_events[binding["event_id"]].add((Q(w["old_t"]),Q(w["old_radius"])))
    assert all(len(value)==1 for value in old_events.values())
    for source in SOURCES:
        events=[b["event_id"] for b in bindings if b["source"]==source and (source,b["node"]) in used]
        assert len(events)==len(set(events))
    inputs += [old/"source_target_contributions.csv",old/"canonical_target_edge_events.csv",known/"source_coverage.csv",known/"actual_target_loops.csv",
               registry/"target_registry.csv",registry/"power_plane_registry.csv",registry/"canonical_edge_registry.csv",
               ROOT/"SHARED_REGISTRY_DOMAINS_PARTIAL_20261021.json",ROOT/"SHARED_REGISTRY_ROOTS_20261021.json",
               ROOT/"ASTRA_SHARED_REGISTRY_EVENT_REVIEW_20261021.md"]
    for rank in range(4):inputs.append(actual/f"actual_domain_registry_rank{rank}.csv")
    for run in (old,known,rootdir,actual):inputs += [run/"source.sha256",run/"retained_input.sha256",run/"output.sha256"]
    dest.mkdir()
    (dest/"bridge_source.py").write_bytes(Path(__file__).read_bytes())
    (dest/"event_identity_bindings.json").write_text(json.dumps(bindings,indent=2)+"\n")
    contract={"scope":"ONLY five originally valid retained domain event identities; no geometry/interval rewrite or new field/root work",
              "baseline_native_domain_exit":2,"sources":sorted(SOURCES),"target":OWN,"epoch":1,
              "state_and_fixed_basis_SHA256":STATE,"exact_function_inputs":"same exact binary64 state/basis/sites/weights; arithmetic enclosures of the same analytic function",
              "old_geometry_bitwise_numeric_records_unchanged":True,"bindings":len(bindings),"used_crossing_nodes_with_existence_completeness_unique_intersection":sum(bool(b["root_identity_witnesses"]) for b in bindings),
              "all_five_event_metadata_complete":True,"effective_domain_arrangements_complete":True,
              "other9211_domain_records_reused":True,"original_target_publication_unchanged":True,
              "source_partition_proved":False,"global_source_image_disjointness_proved":False,
              "new_root_evaluations":0,"new_geometry_evaluations":0,"new_physical_quadrature_queries":0,"full_quadrature_allowed":False}
    (dest/"contract.json").write_text(json.dumps(contract,indent=2)+"\n")
    with (dest/"input.sha256").open("w") as output:
        for path in sorted(set(inputs)):output.write(f"{sha(path)}  {path}\n")
    with (dest/"output.sha256").open("w") as output:
        for name in ("bridge_source.py","event_identity_bindings.json","contract.json"):
            path=dest/name;output.write(f"{sha(path)}  {path}\n")
    (dest/"status.txt").write_text("FIVE_RETAINED_EVENT_IDENTITIES_BOUND_NO_GEOMETRY_REPLAY\n")
    print(json.dumps(contract,indent=2))


if __name__=="__main__":main()

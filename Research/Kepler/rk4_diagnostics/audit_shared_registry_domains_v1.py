"""Saved actual-domain accounting and oriented boundary incidence audit."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys
from audit_shared_regular_q_transfer_v1 import Q, rows, sha

ROOT = Path(__file__).resolve().parent
STATE = "74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "shared_registry_domains_v1_mpi4"
    run = ROOT / name
    assert run.parent == ROOT
    c = json.loads((run / "contract.json").read_text())
    native_exit = int((run / "native_exit_code.txt").read_text())
    assert native_exit in (0,2)
    for filename in ("roots_source_integrity.txt","roots_output_integrity.txt","final_source_integrity.txt","final_input_integrity.txt"):
        text = (run / filename).read_text().splitlines()
        assert text and all(line.endswith(": OK") for line in text)
    for filename in ("source.sha256","retained_input.sha256","binary.sha256","output.sha256"):
        for line in (run / filename).read_text().splitlines():
            digest, path = line.split("  ",1)
            assert sha(Path(path)) == digest
    assert sha(run / "caller_unchanged_state.bin") == STATE
    registry = ROOT / "shared_registry_v2_mpi4"
    slots = rows(registry / "source_target_registry.csv")
    targets = rows(registry / "target_registry.csv")
    indices = {int(r["target"]):int(r["target_index"]) for r in targets}
    edge = {int(r["edge"]):r for r in rows(registry / "canonical_edge_registry.csv")}
    rootsets = defaultdict(list)
    for rank in range(4):
        for row in rows(ROOT / f"shared_registry_roots_v1_mpi4/canonical_event_root_packets_rank{rank}.csv"):
            if row["root_index"] != "-1":
                key = (int(row["canonical_edge"]),int(row["first_target"]),int(row["second_target"]))
                rootsets[key].append(row)
    domains, nodes, boundaries, queries = {}, defaultdict(dict), defaultdict(list), defaultdict(list)
    for rank in range(c["MPI_ranks"]):
        for row in rows(run / f"actual_domain_registry_rank{rank}.csv"):
            slot = int(row["slot_id"])
            assert slot not in domains and indices[int(row["target"])] % c["MPI_ranks"] == rank
            domains[slot] = row
        for row in rows(run / f"domain_nodes_rank{rank}.csv"):
            slot, node = int(row["slot_id"]),int(row["node"])
            assert node not in nodes[slot]
            nodes[slot][node] = row
        for row in rows(run / f"domain_boundaries_rank{rank}.csv"):
            boundaries[int(row["slot_id"])].append(row)
        for row in rows(run / f"actual_membership_queries_rank{rank}.csv"):
            queries[int(row["slot_id"])].append(row)
    assert set(domains) == set(range(9216))
    patch_count=0
    bridge_report=None
    if len(sys.argv)>2:
        bridge=ROOT/sys.argv[2]
        assert bridge.parent==ROOT
        bc=json.loads((bridge/"contract.json").read_text())
        assert bc["all_five_event_metadata_complete"] and bc["state_and_fixed_basis_SHA256"]==STATE
        assert bc["sources"]==[246,253,255,285,287] and bc["target"]==68526
        assert bc["new_root_evaluations"]==bc["new_geometry_evaluations"]==bc["new_physical_quadrature_queries"]==0
        for filename in ("input.sha256","output.sha256"):
            for line in (bridge/filename).read_text().splitlines():
                digest,path=line.split("  ",1)
                assert sha(Path(path))==digest
        bindings=json.loads((bridge/"event_identity_bindings.json").read_text())
        binding={(b["source"],b["node"]):b for b in bindings}
        assert len(binding)==len(bindings)==32
        old=ROOT/"shared_target_v3_mpi4"
        known=ROOT/"shared_boundary_witness_v5_mpi1"
        old_nodes=[r for r in rows(old/"actual_target_nodes.csv") if int(r["source"]) in bc["sources"]]
        known_nodes=[r for r in rows(known/"actual_target_nodes.csv") if int(r["source"]) in bc["sources"]]
        old_boundaries=[r for r in rows(old/"actual_target_boundaries.csv") if int(r["source"]) in bc["sources"]]
        known_boundaries=[r for r in rows(known/"actual_target_boundaries.csv") if int(r["source"]) in bc["sources"]]
        assert old_nodes==known_nodes and old_boundaries==known_boundaries
        old_source={int(r["source"]):r for r in rows(old/"source_target_contributions.csv")}
        coverage={int(r["source"]):r for r in rows(known/"source_coverage.csv")}
        old_events={(int(r["canonical_edge"]),int(r["competitor"])):r for r in rows(old/"canonical_target_edge_events.csv")}
        ids=sorted(indices)
        competitors=[t for t in ids if t!=68526]
        old_nodes_by_key={(int(r["source"]),int(r["node"])):r for r in old_nodes}
        old_boundaries_by_key={(int(r["source"]),int(r["boundary"])):r for r in old_boundaries}
        for (source,node),b in binding.items():
            assert b["epoch"]==1 and b["target"]==68526 and b["original_state_and_basis_SHA256"]==STATE
            original=old_nodes_by_key[source,node]
            assert b["original_vertex"]==int(original["original_vertex_id"])
            active=[competitors[int(j)] for j in original["active_plane_indices"].split(";") if j]
            assert b["active_competitors"]==active
            for witness in b["root_identity_witnesses"]:
                boundary=old_boundaries_by_key[source,witness["old_boundary"]]
                endpoint=witness["which_endpoint"]
                assert int(boundary[endpoint])==node
                parameter="parameter_lo" if endpoint=="first" else "parameter_hi"
                assert witness["old_t"]==boundary[parameter] and witness["old_radius"]==boundary[parameter+"_radius"]
                key=(int(boundary["canonical_edge"]),min(68526,active[0]),max(68526,active[0]))
                original_event=old_events[key[0],active[0]]
                assert witness["old_root_existence_complete_event_record"]==original_event
                assert original_event["complete"]=="1" and source in (int(original_event["left_source"]),int(original_event["right_source"]))
                assert witness["parameter_orientation_transform"]=="identity_original_canonical_t"
                assert witness["plane_sign_relative_to_new"]==(1 if 68526<active[0] else -1)
                assert (witness["new_first_target"],witness["new_second_target"])==key[1:]
                n,e=Q(witness["old_t"]),Q(witness["old_radius"])
                matches=[r for r in rootsets[key] if max(n-e,Q(r["t"])-Q(r["t_radius"]))<=min(n+e,Q(r["t"])+Q(r["t_radius"]))]
                assert len(matches)==witness["unique_intersection_count"]==1
                assert matches[0]==witness["new_complete_root_packet"]
                assert b["event_id"]=="e;"+";".join(map(str,(*key,int(matches[0]["root_index"]))))
            if b["event_id"].startswith("e;"):assert b["root_identity_witnesses"]
        for source in bc["sources"]:
            slot=source*24+indices[68526]
            row=domains[slot]
            assert row["geometry_complete"]=="0" and row["reason"]=="used edge crossing binds the actual canonical root index"
            assert row["target"]=="68526" and row["reused"]=="1" and not boundaries[slot] and not nodes[slot]
            assert old_source[source]["geometry_complete"]=="1"
            assert coverage[source]["coverage_complete"]==coverage[source]["reused_actual_domain"]=="1" and coverage[source]["empty"]=="0"
            for n in [r for r in old_nodes if int(r["source"])==source]:
                node=int(n["node"]);b=binding[source,node]
                nodes[slot][node]={"slot_id":str(slot),"node":str(node),"global_vertex":n["original_vertex_id"],
                                   **{k:n[k] for k in ("x","x_radius","y","y_radius")},"event_id":b["event_id"],
                                   "active_competitors":"".join(str(i)+";" for i in b["active_competitors"])}
            for b in [r for r in old_boundaries if int(r["source"])==source]:
                first,last=int(b["first"]),int(b["last"]);canonical=int(b["canonical_edge"])
                a=z=0
                if canonical<0:
                    common=set(binding[source,first]["active_competitors"])&set(binding[source,last]["active_competitors"])
                    assert len(common)==1
                    other=common.pop();a,z=min(68526,other),max(68526,other)
                boundaries[slot].append({"slot_id":str(slot),"source":str(source),"target":"68526","boundary":b["boundary"],
                                         "first_node":str(first),"last_node":str(last),"first_event":binding[source,first]["event_id"],
                                         "last_event":binding[source,last]["event_id"],"canonical_edge":str(canonical),
                                         "first_plane_target":str(a),"second_plane_target":str(z),"parameter_lo":b["parameter_lo"],
                                         "lo_radius":b["parameter_lo_radius"],"parameter_hi":b["parameter_hi"],"hi_radius":b["parameter_hi_radius"],
                                         "external_support_boundary":str(int(canonical>=0 and edge[canonical]["physical_external_boundary"]=="1")),
                                         "implicit_whole_source_boundary":"0"})
            row["geometry_complete"]="1";row["reason"]="same_original_geometry_with_certified_event_identity_bridge"
            patch_count+=1
        bridge_report={"path":str(bridge),"contract":bc,"independent_saved_existence_complete_set_unique_intersection_checks_pass":True}
    failed, reasons = [],Counter()
    nempty = nwhole = nclipped = nqueries = windows = 0
    plane_interfaces = defaultdict(list)
    source_edges = defaultdict(list)
    side_interfaces = defaultdict(list)
    external = []
    event_obligations = []
    for slot, row in domains.items():
        old = slots[slot]
        source, target = int(row["source"]),int(row["target"])
        assert source == int(old["source"]) and target == int(old["target"])
        assert row["component"] == old["component"]
        assert row["coarse_classification"] == old["coarse_classification"]
        assert row["retained_resolution"] == old["resolution"]
        assert len(queries[slot]) == int(row["winding_queries"])
        for index, q in enumerate(queries[slot]):
            assert int(q["query"]) == index
            assert Q(q["x_radius"]) >= 0 and Q(q["y_radius"]) >= 0
        assert sum(int(q["root_windows"]) for q in queries[slot]) == int(row["winding_windows"])
        nqueries += len(queries[slot]);windows += int(row["winding_windows"])
        complete = row["geometry_complete"] == "1"
        empty = row["empty"] == "1"
        whole = row["coarse_classification"] == "1" or row["retained_resolution"] == "3"
        if target == 68526:
            assert row["reused"] == "1" and not queries[slot]
        if not complete:
            failed.append({"slot_id":slot,"source":source,"target":target,"reason":row["reason"]})
            reasons[row["reason"]] += 1
            assert not boundaries[slot] and not nodes[slot]
            continue
        if empty:
            nempty += 1
            assert not boundaries[slot] and not nodes[slot]
            continue
        nwhole += whole;nclipped += not whole
        assert len(boundaries[slot]) == (3 if whole else int(row["boundaries"]))
        assert whole or int(row["loops"]) > 0
        outgoing,incoming = Counter(),Counter()
        for j,b in enumerate(boundaries[slot]):
            assert int(b["boundary"]) == j and int(b["source"]) == source and int(b["target"]) == target
            a,z = int(b["first_node"]),int(b["last_node"])
            first,last = b["first_event"],b["last_event"]
            assert first != last and not first.startswith("unused;") and not last.startswith("unused;")
            outgoing[a] += 1;incoming[z] += 1
            assert b["implicit_whole_source_boundary"] == str(int(whole))
            if not whole:
                assert a in nodes[slot] and z in nodes[slot]
                assert first == nodes[slot][a]["event_id"] and last == nodes[slot][z]["event_id"]
            canonical = int(b["canonical_edge"])
            if canonical >= 0:
                e = edge[canonical]
                assert source in (int(e["first_source"]),int(e["second_source"]))
                is_external = e["physical_external_boundary"] == "1"
                assert b["external_support_boundary"] == str(int(is_external))
                lo,hi,le,he = Q(b["parameter_lo"]),Q(b["parameter_hi"]),Q(b["lo_radius"]),Q(b["hi_radius"])
                assert le >= 0 and he >= 0 and 0 <= lo-le <= lo+le <= 1 and 0 <= hi-he <= hi+he <= 1
                ascending = lo < hi
                assert lo+le < hi-he if ascending else hi+he < lo-le
                orientation = int(e["first_orientation"] if source == int(e["first_source"]) else e["second_orientation"])
                assert ascending == (orientation == 1)
                for event,t,radius in ((first,lo,le),(last,hi,he)):
                    bits = event.split(";")
                    if bits[0] == "v":
                        assert len(bits) == 2 and radius == 0 and t in (0,1)
                        assert int(bits[1]) == int(e["first_vertex"] if t == 0 else e["last_vertex"])
                    elif bits[0] == "e":
                        assert len(bits) == 5 and int(bits[1]) == canonical
                        key = (canonical,int(bits[2]),int(bits[3]))
                        root = rootsets[key][int(bits[4])]
                        if bridge_report and source in bridge_report["contract"]["sources"] and target==68526:
                            matches=[r for r in rootsets[key] if max(t-radius,Q(r["t"])-Q(r["t_radius"]))<=min(t+radius,Q(r["t"])+Q(r["t_radius"]))]
                            assert len(matches)==1 and matches[0]==root
                        else:
                            assert t == Q(root["t"]) and radius == Q(root["t_radius"])
                    else:
                        # A true multi-plane vertex on a source edge needs its
                        # own event binding; never silently call this proved.
                        event_obligations.append({"slot_id":slot,"boundary":j,"event":event})
                start,end = (first,last) if ascending else (last,first)
                low,high,lr,hr = (lo,hi,le,he) if ascending else (hi,lo,he,le)
                fragment={"slot_id":slot,"source":source,"target":target,"edge":canonical,
                          "first_event":start,"last_event":end,"lo":low,"hi":high,"lo_radius":lr,"hi_radius":hr,
                          "orientation":orientation}
                source_edges[source,canonical].append(fragment)
                side_interfaces[canonical,target,start,end].append(fragment)
                if is_external:
                    external.append({k:v for k,v in fragment.items() if k not in ("lo","hi","lo_radius","hi_radius")})
            else:
                assert not whole and b["external_support_boundary"] == "0"
                x,y = int(b["first_plane_target"]),int(b["second_plane_target"])
                assert x < y and target in (x,y)
                plane_interfaces[source,x,y,min(first,last),max(first,last)].append((target,first,last,slot,j))
        assert outgoing == incoming and all(value == 1 for value in outgoing.values())
        if not whole:
            assert set(nodes[slot]) == set(range(len(nodes[slot])))
            for node in nodes[slot].values():
                assert Q(node["x_radius"]) >= 0 and Q(node["y_radius"]) >= 0

    expected=[c[k] for k in ("incomplete_domains","complete_empty_or_excluded","complete_whole_source","complete_clipped_domains","actual_membership_queries","actual_membership_root_windows")]
    expected[0]-=patch_count;expected[3]+=patch_count
    assert (len(failed),nempty,nwhole,nclipped,nqueries,windows) == tuple(expected)
    assert len(failed)+nempty+nwhole+nclipped == 9216
    assert c["all_domain_arrangements_complete"] == (c["incomplete_domains"]==0) and native_exit == (2 if c["incomplete_domains"] else 0)
    mismatched_planes=[]
    for key, values in plane_interfaces.items():
        if not(len(values)==2 and values[0][0]!=values[1][0] and values[0][1:3]==(values[1][2],values[1][1])):
            mismatched_planes.append({"key":key,"incidences":values})
    chains=[]
    expected_source_edges={(int(e[name]),key) for key,e in edge.items() for name in ("first_source","second_source") if int(e[name]) >= 0}
    assert len(expected_source_edges) == 1152
    for key in expected_source_edges:
        source,canonical=key;e=edge[canonical]
        raw_fragments=source_edges[key]
        by_start={r["first_event"]:r for r in raw_fragments}
        # Certified event identities build the path. Midpoint sorting and
        # pairwise interval-overlap equivalence are not ownership proofs.
        fragments=[];current="v;"+e["first_vertex"];seen=set()
        while current in by_start and current not in seen:
            seen.add(current);fragment=by_start[current];fragments.append(fragment);current=fragment["last_event"]
        valid=bool(fragments) and len(by_start)==len(raw_fragments)==len(fragments)
        if valid:
            valid=fragments[0]["first_event"]=="v;"+e["first_vertex"] and fragments[-1]["last_event"]=="v;"+e["last_vertex"]
            valid=valid and fragments[0]["lo"]==0 and fragments[0]["lo_radius"]==0 and fragments[-1]["hi"]==1 and fragments[-1]["hi_radius"]==0
            for a,b in zip(fragments,fragments[1:]):
                valid=valid and a["last_event"]==b["first_event"] and max(a["hi"]-a["hi_radius"],b["lo"]-b["lo_radius"])<=min(a["hi"]+a["hi_radius"],b["lo"]+b["lo_radius"])
        if not valid:
            chains.append({"source":source,"edge":canonical,"fragments":len(fragments)})
    paired_mismatches=[]
    for key, values in side_interfaces.items():
        e=edge[key[0]]
        if e["physical_external_boundary"] == "0":
            if not(len(values)==2 and values[0]["source"]!=values[1]["source"] and values[0]["orientation"]+values[1]["orientation"]==0):
                paired_mismatches.append({"key":key,"sources":[r["source"] for r in values]})
    structural=not(failed or mismatched_planes or chains or paired_mismatches or event_obligations)
    for flag in ("source_partition_boundary_cancellation_proved","global_source_image_disjointness_proved","full_quadrature_allowed","accepted_RK4_evolution","long_Kepler_validated"):
        assert c[flag] is False
    assert c["new_physical_quadrature_queries"] == 0
    print(json.dumps({"scope":"Saved actual-domain packet accounting, closed boundary cycles, canonical event reuse and oriented source/target incidence cancellation; no independent curved-map or membership evaluation",
                      "run":str(run),"native_exit_code":native_exit,"contract":c,"original_state_SHA256":STATE,
                      "separate_legacy_event_identity_bridge":bridge_report,
                      "effective_all9216_domain_arrangements_complete":not failed,
                      "all9216_actual_domain_records_accounted":True,"known_target68526_reuse_records_accounted":True,
                      "known_target68526_reuse_metadata_complete":all(r["geometry_complete"]=="1" for r in domains.values() if r["target"]=="68526"),
                      "packet_and_failure_accounting_pass":True,"incomplete_reasons":dict(reasons),"incomplete_domains":failed,
                      "internal_power_interface_groups":len(plane_interfaces),"mismatched_power_interfaces":mismatched_planes,
                      "incomplete_source_edge_chains":chains,"mismatched_paired_edge_fragments":paired_mismatches,
                      "unbound_multiplane_edge_events":event_obligations,
                      "external_actual_boundary_fragments":external,
                      "saved_event_source_partition_and_edge_cancellation_complete":structural,
                      "global_source_image_disjointness_complete":False,"all_source_numerical_conservation_complete":False,
                      "actual_all_target_geometric_registry_complete":False,"full_quadrature_allowed":False},indent=2))


if __name__ == "__main__":
    main()

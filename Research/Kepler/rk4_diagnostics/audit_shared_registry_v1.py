"""Exact saved-input broad-phase and incidence audit; no field evaluation."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cmp_to_key
import json
from pathlib import Path
import struct
import sys

from audit_shared_regular_q_transfer_v1 import Q, rows, sha

ROOT = Path(__file__).resolve().parent
STATE = "74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "shared_registry_v1_mpi4"
    run = ROOT / name
    assert run.parent == ROOT
    c = json.loads((run / "contract.json").read_text())
    assert (run / "native_exit_code.txt").read_text().strip() == "0"
    for filename in ("final_source_integrity.txt", "final_input_integrity.txt",
                     "originating_source_integrity.txt", "originating_output_integrity.txt"):
        lines = (run / filename).read_text().splitlines()
        assert lines and all(line.endswith(": OK") for line in lines)
    for manifest in ("source.sha256", "retained_input.sha256", "binary.sha256", "output.sha256"):
        for line in (run / manifest).read_text().splitlines():
            digest, path = line.split("  ", 1)
            assert sha(Path(path)) == digest
    assert sha(run / "caller_unchanged_state.bin") == STATE
    earlier = json.loads((ROOT / "SHARED_MISSING_SOURCE_LEDGER_20261021.json").read_text())
    assert earlier["target_published"] and not earlier["ORIGINAL_tau_failures"]

    original_sites = {int(r["id"]): r for r in rows(ROOT / "curved_map_algebra_v2/advected_sites.csv")}
    original_weights = {int(r["id"]): Fraction(r["weight"]) for r in rows(ROOT / "curved_events_v9/site_weights.csv")}
    targets = rows(run / "target_registry.csv")
    assert len(targets) == 24
    ids = sorted(original_sites)
    site = {}
    for index, row in enumerate(targets):
        target = int(row["target"])
        assert target == ids[index] and int(row["target_index"]) == index
        old = original_sites[target]
        xy = (Q(row["x"]), Q(row["y"]))
        assert xy == (Fraction.from_float(float(old["mapped_x"])), Fraction.from_float(float(old["mapped_y"])))
        assert Q(row["weight"]) == original_weights[target] == 0
        assert row["strict_site_interior"] == "1"
        assert row["partition_tie_rule"] == "smallest_target_id_owns_exact_equality"
        assert row["physical_domain"] == "source_patch_intersection_only_no_fill"
        assert row["all12_tau_available"] == str(int(target == 68526))
        site[target] = xy + (Q(row["weight"]),)

    def plane(a, b):
        ax, ay, aw = site[a]
        bx, by, bw = site[b]
        return bx-ax, by-ay, (bx*bx+by*by-ax*ax-ay*ay+aw-bw)/2

    def half(n):
        return 0 if n[1] > 0 or n[1] == 0 and n[0] > 0 else 1

    def cross(a, b):
        return a[0]*b[1]-a[1]*b[0]

    def compare(a, b):
        if half(a) != half(b):
            return half(a)-half(b)
        return -1 if cross(a, b) > 0 else 1 if cross(a, b) < 0 else 0

    for row in targets:
        target = int(row["target"])
        ns = sorted([plane(target, other)[:2] for other in ids if other != target], key=cmp_to_key(compare))
        unique = []
        for n in ns:
            assert n != (0, 0)
            if not unique or half(unique[-1]) != half(n) or cross(unique[-1], n) != 0:
                unique.append(n)
        bounded = len(unique) >= 3 and all(cross(unique[j], unique[(j+1) % len(unique)]) > 0 for j in range(len(unique)))
        assert int(row["bounded"]) == bounded
        x, y, _ = site[target]
        assert all(nx*x+ny*y-b < 0 for nx, ny, b in [plane(target, o) for o in ids if o != target])
    pairplanes = rows(run / "power_plane_registry.csv")
    assert len(pairplanes) == 276
    assert {(int(r["first_target"]), int(r["second_target"])) for r in pairplanes} == {(a,b) for a in ids for b in ids if a < b}
    for row in pairplanes:
        a, b = int(row["first_target"]), int(row["second_target"])
        true = plane(a, b)
        for k, name in enumerate(("nx", "ny", "b")):
            assert abs(Q(row[name])-true[k]) <= Q(row[name+"_radius"])
        assert row["first_owns_equality"] == "1"

    packet = (ROOT / "material_p3_finite_v2/mpi1_cpu/curved_map_control.bin").read_bytes()
    h = struct.unpack_from("<8Q", packet)
    assert h[2] == 384 and (len(packet)-160) % 384 == 0
    width = (len(packet)-160) // 384
    geometry = []
    for source in range(384):
        start = 160 + source*width
        cell, component, modes = struct.unpack_from("<3Q", packet, start)
        assert modes == 12
        vertices = struct.unpack_from("<6d", packet, start+24)
        geometry.append((cell, component, tuple(Fraction.from_float(x) for x in vertices)))
    images = rows(run / "source_image_registry.csv")
    bounds = {int(r["source"]): r for r in rows(ROOT / "shared_real_bound_v1_mpi4/source_real_bounds.csv")}
    assert len(images) == 384
    box = {}
    for source, row in enumerate(images):
        assert int(row["source"]) == source
        assert int(row["cell"]) == geometry[source][0]
        assert Fraction(bounds[source]["tight_L"]) < Q(row["L"]) < 1
        intervals = []
        for k, name in enumerate(("x", "y")):
            center, radius = Q(row[name]), Q(row[name+"_radius"])
            displacement = Q(row["displacement_"+name])
            assert displacement > Fraction(bounds[source]["displacement_"+name+"_bound"])
            coords = geometry[source][2][k::2]
            assert center-radius <= min(coords)-displacement
            assert center+radius >= max(coords)+displacement
            intervals.append((center-radius, center+radius))
        box[source] = intervals

    def interval(source, target, competitor):
        nx, ny, b = plane(target, competitor)
        x, y = box[source]
        return min(nx*x[0],nx*x[1])+min(ny*y[0],ny*y[1])-b, max(nx*x[0],nx*x[1])+max(ny*y[0],ny*y[1])-b

    pairs = rows(run / "source_target_registry.csv")
    raw = (run / "pair_registry.bin").read_bytes()
    assert struct.unpack_from("<8Q", raw) == (0x4753465245473031,1,1,72,9216,24,384,c["MPI_ranks"])
    assert len(raw) == 64+72*9216 and len(pairs) == 9216
    prior = {int(r["source"]): r for r in rows(ROOT / "shared_missing_source_v1_mpi1/all384_source_completion.csv")}
    wholes = {int(r["source"]): r for r in rows(ROOT / "shared_atlas_v4_mpi4/whole_inside_source_integrals.csv")}
    active = {}
    per_target = defaultdict(Counter)
    for index, row in enumerate(pairs):
        source, ti = divmod(index, 24)
        target = ids[ti]
        coarse, resolution = int(row["coarse_classification"]), int(row["resolution"])
        expected = (index,source,geometry[source][0],target,ti,index % c["MPI_ranks"],coarse,resolution,geometry[source][1],0,int(row["retained_points"]),int(row["retained_panels"]))
        assert struct.unpack_from("<4Q2I4i2Q", raw, 64+72*index) == expected
        assert tuple(int(row[k]) for k in ("slot_id","source","cell","target","target_index","owner_rank","component")) == (index,source,geometry[source][0],target,ti,index % c["MPI_ranks"],geometry[source][1])
        values = [interval(source, target, other) for other in ids if other != target]
        if coarse == 0:
            assert any(lo > 0 for lo, hi in values)
        elif coarse == 1:
            assert all(hi < 0 for lo, hi in values)
        else:
            assert coarse == 2  # Conservatively unresolved, never inferred empty.
        if target == 68526:
            old = prior[source]
            expected_resolution = 0 if old["classification"] == "0" else 3 if old["classification"] == "1" else 4 if int(old["new_positive_points"]) else 2
            assert resolution == expected_resolution
            expected_points = int(wholes[source]["positive_points"]) if resolution == 3 else int(old["new_positive_points"])
            expected_panels = int(wholes[source]["panels"]) if resolution == 3 else int(old["planned_panels"])
            assert (int(row["retained_points"]),int(row["retained_panels"])) == (expected_points, expected_panels)
        else:
            assert resolution == (1 if coarse else 0)
            assert row["retained_points"] == row["retained_panels"] == "0"
        assert row["actual_overlap_integral_available"] == str(int(resolution == 4))
        if resolution in (1,3,4):
            active[index] = row
        per_target[target]["pending_actual_atlas"] += resolution == 1
        per_target[target]["retained_individual_clipped"] += resolution == 4
        per_target[target]["retained_whole_group_members"] += resolution == 3
        per_target[target]["coarse_whole"] += coarse == 1
        per_target[target]["coarse_candidate"] += coarse == 2
        per_target[target]["retained_positive_queries"] += int(row["retained_points"])

    outgoing = rows(run / "source_outgoing_incidence.csv")
    incoming = rows(run / "target_incoming_incidence.csv")
    assert len(outgoing) == len(incoming) == len(active) == c["unique_active_incidence_slots"]
    for direction, sign in ((outgoing,-1),(incoming,1)):
        assert len({int(r["slot_id"]) for r in direction}) == len(direction)
        for row in direction:
            index = int(row["slot_id"])
            slot = active[index]
            assert row["source"] == slot["source"] and row["target"] == slot["target"]
            assert int(row["coefficient"]) == sign and int(row["slot_byte_offset"]) == 64+72*index
            assert row["numerical_value_available"] == slot["actual_overlap_integral_available"]
    group = rows(run / "retained_reduction_groups.csv")
    assert len(group) == 1 and group[0]["members"] == "240;244;248;252"
    assert group[0]["arithmetic_value_instances"] == "1" and group[0]["individual_source_closure_proved"] == "0"

    edges = rows(run / "canonical_edge_registry.csv")
    assert len(edges) == c["canonical_edges"] == 616
    edge_lookup = {int(r["edge"]): r for r in edges}
    assert len(edge_lookup) == 616
    paired = external = 0
    for row in edges:
        assert row["interval_lo"] == "0" and row["interval_hi"] == "1"
        if row["physical_external_boundary"] == "1":
            external += 1
            assert row["second_source"] == "-1" and row["fragment_state"] == "unsubdivided_candidate_only"
        else:
            paired += 1
            assert int(row["first_orientation"])+int(row["second_orientation"]) == 0
        vertices = geometry[int(row["first_source"])][2]
        a = int(row["first_local"])
        first = vertices[2*a:2*a+2]
        last = vertices[2*((a+1) % 3):2*((a+1) % 3)+2]
        # Factory canonical orientation sorts reference vertex coordinates.
        assert int(row["first_orientation"]) == (1 if first < last else -1)
    assert (paired,external) == (c["paired_edges"],c["external_edges"])
    assert paired*2+external == 3*384
    expected_external = {(int(e["edge"]),int(e["first_source"]),ids[t]) for e in edges if e["physical_external_boundary"] == "1" for t in range(24) if int(e["first_source"])*24+t in active}
    actual_external = rows(run / "external_edge_target_candidates.csv")
    assert {(int(r["edge"]),int(r["source"]),int(r["target"])) for r in actual_external} == expected_external
    assert len(actual_external) == len(expected_external) == c["external_edge_candidate_incidences"]
    assert all(r["fragment_state"] == "pending_actual_edge_partition" for r in actual_external)

    old_roots = {(int(r["canonical_edge"]),min(68526,int(r["competitor"])),max(68526,int(r["competitor"]))) for r in rows(ROOT / "shared_coverage_v2_mpi4/canonical_event_root_packets.csv") if r["complete"] == "1"}
    jobs = rows(run / "canonical_root_work_registry.csv")
    keys = {(int(r["edge"]),int(r["first_target"]),int(r["second_target"])) for r in jobs}
    assert len(keys) == len(jobs)
    required = set()
    for edge in edges:
        key = int(edge["edge"])
        sources = [int(edge["first_source"])]
        if edge["physical_external_boundary"] == "0":
            sources.append(int(edge["second_source"]))
        for t, target in enumerate(ids):
            if not any(s*24+t in active for s in sources):
                continue
            for other in ids:
                if other == target:
                    continue
                if any(lo <= 0 <= hi for lo,hi in [interval(s,target,other) for s in sources]):
                    required.add((key,min(target,other),max(target,other)))
    assert required <= keys  # Interval arithmetic may conservatively keep extra jobs.
    for index, row in enumerate(jobs):
        assert int(row["work_id"]) == index
        key = int(row["edge"]),int(row["first_target"]),int(row["second_target"])
        reused = key in old_roots
        assert row["retained_complete_roots"] == str(int(reused))
        assert int(row["new_root_evaluation_cap"]) == (0 if reused else 256)
        assert key[0] in edge_lookup and key[1] < key[2] and key[1] in site and key[2] in site
    assert sum(r["retained_complete_roots"] == "1" for r in jobs) == c["retained_root_jobs"]
    assert len(jobs)-c["retained_root_jobs"] == c["new_deduplicated_root_jobs"]
    cost = {r["item"]: int(r["count"]) for r in rows(run / "cost_preflight.csv")}
    for item in ("new_panel_count","new_positive_query_count","new_positive_query_bytes"):
        assert cost[item] == -1
    assert cost["registry_binary_bytes"] == len(raw)
    assert cost["new_root_window_evaluation_cap"] == c["new_deduplicated_root_jobs"]*256
    assert cost["full_quadrature_launch_allowed"] == 0
    for flag in ("actual_all_target_overlap_geometry_complete","source_partition_boundary_cancellation_complete",
                 "global_source_image_disjointness_complete","all_source_numerical_conservation_complete",
                 "external_actual_fragments_complete","other23_physical_normalizers_complete",
                 "full_quadrature_launch_allowed","accepted_RK4_evolution","long_Kepler_validated"):
        assert c[flag] is False
    assert c["new_field_queries"] == 0
    print(json.dumps({"scope":"Exact saved-input power-plane/image broad phase, binary registry, incidence transpose and root-work accounting; no independent curved map/field evaluation",
                      "run":str(run),"contract":c,"original_state_SHA256":STATE,
                      "all9216_pairs_audited":True,"all276_exact_pair_planes_enclosed":True,
                      "same_immutable_slot_source_outgoing_target_incoming":True,
                      "known_target68526_reused_without_reintegration":True,
                      "retained_whole_group_once_individual_numerical_closure_pending":True,
                      "all_required_coarse_root_jobs_registered":True,
                      "cost_preflight":cost,
                      "per_target_preflight":[{"target":target,**dict(per_target[target]),"new_panels":0 if target==68526 else None,"new_positive_queries":0 if target==68526 else None,"physical_normalizer_available":target==68526} for target in ids],
                      "preliminary_registry_audit_pass":True,"complete_geometric_registry":False,
                      "full_quadrature_allowed":False},indent=2))


if __name__ == "__main__":
    main()

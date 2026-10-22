"""Read-only registry/root packet accounting; no numerical root replay."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys
from audit_shared_regular_q_transfer_v1 import Q, rows, sha

ROOT = Path(__file__).resolve().parent
STATE = "74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "shared_registry_roots_v1_mpi4"
    run = ROOT / name
    assert run.parent == ROOT
    contract = json.loads((run / "contract.json").read_text())
    native_exit = int((run / "native_exit_code.txt").read_text())
    assert native_exit in (0, 2)
    for filename in ("registry_source_integrity.txt", "registry_output_integrity.txt",
                     "final_source_integrity.txt", "final_input_integrity.txt"):
        lines = (run / filename).read_text().splitlines()
        assert lines and all(line.endswith(": OK") for line in lines)
    for filename in ("source.sha256", "retained_input.sha256", "binary.sha256", "output.sha256"):
        for line in (run / filename).read_text().splitlines():
            digest, path = line.split("  ", 1)
            assert sha(Path(path)) == digest
    assert sha(run / "caller_unchanged_state.bin") == STATE
    registry = ROOT / "shared_registry_v2_mpi4"
    plan = rows(registry / "canonical_root_work_registry.csv")
    old = defaultdict(list)
    for row in rows(ROOT / "shared_coverage_v2_mpi4/canonical_event_root_packets.csv"):
        assert row["complete"] == "1"
        other = int(row["competitor"])
        key = int(row["canonical_edge"]), min(68526, other), max(68526, other)
        old[key].append(row)
    results, packets, witnesses = {}, defaultdict(list), {}
    for rank in range(contract["MPI_ranks"]):
        for row in rows(run / f"canonical_root_results_rank{rank}.csv"):
            work = int(row["work_id"])
            assert work % contract["MPI_ranks"] == rank == int(row["owner_rank"])
            assert work not in results
            results[work] = row
        for row in rows(run / f"canonical_event_root_packets_rank{rank}.csv"):
            work = int(row["work_id"])
            assert work % contract["MPI_ranks"] == rank
            packets[work].append(row)
        for row in rows(run / f"failed_root_windows_rank{rank}.csv"):
            work = int(row["work_id"])
            assert work % contract["MPI_ranks"] == rank and work not in witnesses
            witnesses[work] = row
    assert set(results) == set(packets) == set(range(2452))
    new = reused = incomplete = windows = roots = endpoints = 0
    reason = Counter()
    failed = []
    for work, job in enumerate(plan):
        assert int(job["work_id"]) == work
        row = results[work]
        for a, b in (("edge", "canonical_edge"), ("first_target", "first_target"), ("second_target", "second_target")):
            assert job[a] == row[b]
        retained = job["retained_complete_roots"] == "1"
        assert row["retained"] == str(int(retained))
        complete = row["complete"] == "1"
        rpackets = packets[work]
        count = int(row["roots"])
        for packet in rpackets:
            for field in ("canonical_edge", "first_target", "second_target", "retained", "complete"):
                assert packet[field] == row[field]
            assert int(packet["root_count"]) == count
        if count == 0:
            assert len(rpackets) == 1 and rpackets[0]["root_index"] == "-1"
        else:
            assert len(rpackets) == count
            assert [int(p["root_index"]) for p in rpackets] == list(range(count))
            for index, packet in enumerate(rpackets):
                t, radius = Q(packet["t"]), Q(packet["t_radius"])
                assert radius >= 0 and 0 <= t-radius <= t+radius <= 1
                endpoint = int(packet["endpoint"])
                assert endpoint in (-1,0,1)
                if endpoint >= 0:
                    assert t == endpoint and radius == 0
                if complete and index:
                    previous = rpackets[index-1]
                    assert Q(previous["t"])+Q(previous["t_radius"]) < t-radius
            assert sum(int(p["endpoint"]) >= 0 for p in rpackets) == int(row["endpoint_roots"])
        if retained:
            reused += 1
            key = int(row["canonical_edge"]), int(row["first_target"]), int(row["second_target"])
            previous = old[key]
            assert complete and len(previous) == len(rpackets)
            for a, b in zip(previous, rpackets):
                assert a["root_index"] == b["root_index"] and a["root_count"] == b["root_count"]
                assert Q(a["t"]) == Q(b["t"]) and Q(a["t_radius"]) == Q(b["t_radius"])
                assert a["endpoint"] == b["endpoint"]
        else:
            new += 1
            counter = int(row["root_windows"])
            assert 0 <= counter <= 257  # Original postfix cap guard counts its rejected attempt.
            if complete:
                assert counter <= 256
            windows += counter
        roots += count
        endpoints += int(row["endpoint_roots"])
        if not complete:
            incomplete += 1
            assert work in witnesses
            witness = witnesses[work]
            assert witness["reason"] == row["reason"]
            for field in ("canonical_edge", "first_target", "second_target"):
                assert witness[field] == row[field]
            reason[row["reason"]] += 1
            failed.append({"work_id":work, "canonical_edge":int(row["canonical_edge"]),
                           "first_target":int(row["first_target"]), "second_target":int(row["second_target"]),
                           "windows":int(row["root_windows"]),"reason":row["reason"]})
        else:
            assert work not in witnesses
    assert (new,reused,incomplete,windows,roots,endpoints) == tuple(contract[k] for k in ("new_jobs","retained_jobs","incomplete_jobs","new_root_windows","retained_and_new_roots","endpoint_roots"))
    assert new == 2289 and reused == 163 and new+reused == contract["planned_jobs"] == 2452
    assert contract["all_registered_roots_complete"] == (not incomplete)
    assert native_exit == (2 if incomplete else 0)
    for field in ("actual_all_target_domains_complete", "source_partition_complete", "full_quadrature_allowed", "accepted_RK4_evolution", "long_Kepler_validated"):
        assert contract[field] is False
    assert contract["new_physical_quadrature_queries"] == 0
    print(json.dumps({"scope":"Saved canonical root packets, work ownership, retained-bit preservation, finite root intervals and fixed-cap result accounting; no independent root existence/completeness proof replay",
                      "run":str(run),"native_exit_code":native_exit,"contract":contract,
                      "original_state_SHA256":STATE,"all2452_jobs_bound_to_audited_registry":True,
                      "all163_retained_root_sets_unchanged":True,
                      "packet_and_failure_accounting_pass":True,
                      "incomplete_reasons":dict(reason),"incomplete_jobs":failed,
                      "actual_all_target_geometric_registry_complete":False,
                      "full_quadrature_allowed":False}, indent=2))


if __name__ == "__main__":
    main()

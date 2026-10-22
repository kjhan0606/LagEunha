"""One read-only progress observation; call at most once per 60 seconds.

Reads source IDs from new, fully flushed records only. Never evaluates a
field, restarts a job, or modifies frozen run files. The small index lives
outside the frozen run. A flushed point count is not an integral verdict.
"""
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "shared_clipped_target_v8_mpi4"
INDEX = ROOT / "shared_clipped_target_v8_progress_index.json"
PLAN = ROOT / "shared_clipped_target_v7_mpi4/all384_source_completion.csv"
WIDTH = 1008


def main():
    plan_hash = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    with PLAN.open() as handle:
        planned = {int(row["source"]): int(row["planned_panels"]) * 256
                   for row in csv.DictReader(handle) if int(row["planned_panels"])}
    assert len(planned) == 26 and sum(planned.values()) == 222976
    if INDEX.exists():
        index = json.loads(INDEX.read_text())
        assert index["run"] == str(RUN) and index["plan_SHA256"] == plan_hash
        assert index["record_bytes"] == WIDTH
    else:
        index = {"run": str(RUN), "plan_SHA256": plan_hash,
                 "record_bytes": WIDTH, "ranks": {}}
    counts = Counter()
    by_rank = []
    for rank in range(4):
        packet = RUN / f"NEW_clipped_positive_queries_rank{rank}.bin"
        stat = packet.stat()
        count = stat.st_size // WIDTH
        identity = [stat.st_dev, stat.st_ino]
        previous = index["ranks"].get(str(rank), {
            "identity": identity, "records": 0, "source_counts": {}})
        assert previous["identity"] == identity and count >= previous["records"]
        observed = Counter({int(s): n for s, n in previous["source_counts"].items()})
        with packet.open("rb") as handle:
            handle.seek(previous["records"] * WIDTH)
            for first in range(previous["records"], count, 4096):
                number = min(4096, count - first)
                raw = handle.read(number * WIDTH)
                assert len(raw) == number * WIDTH
                for point in range(number):
                    source = struct.unpack_from("<Q", raw, point * WIDTH)[0]
                    assert source in planned and source % 4 == rank
                    observed[source] += 1
        assert sum(observed.values()) == count
        assert all(n <= planned[s] for s, n in observed.items())
        index["ranks"][str(rank)] = {
            "identity": identity, "records": count, "source_counts": dict(observed)}
        counts.update(observed)
        by_rank.append(count)
    index["UTC"] = datetime.now(timezone.utc).isoformat()
    temporary = INDEX.with_name(INDEX.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(index, indent=2) + "\n")
    temporary.replace(INDEX)
    terminal = RUN / "native_exit_code.txt"
    status = RUN / "status.txt"
    print(json.dumps({
        "UTC": index["UTC"],
        "scope": "flushed record source index only; no integral or physical validation verdict",
        "record_estimates_by_rank": by_rank,
        "record_estimate_total": sum(by_rank),
        "planned_new_records": sum(planned.values()),
        "fully_recorded_sources": sum(counts[s] == n for s, n in planned.items()),
        "planned_sources": len(planned),
        "remaining_sources": [
            {"source": s, "rank": s % 4, "recorded": counts[s], "planned": n,
             "remaining": n - counts[s]}
            for s, n in sorted(planned.items()) if counts[s] < n],
        "native_exit_code": int(terminal.read_text()) if terminal.exists() else None,
        "terminal_status": status.read_text().strip() if status.exists() else None,
    }, indent=2))


if __name__ == "__main__":
    main()

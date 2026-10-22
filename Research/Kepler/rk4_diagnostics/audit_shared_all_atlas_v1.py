"""Read-only all-target atlas packet/count/cost accounting."""
from collections import Counter, defaultdict
from pathlib import Path
import json
import struct
import sys
from audit_shared_regular_q_transfer_v1 import decode_quad, rows, sha

ROOT=Path(__file__).resolve().parent


def main():
    name=sys.argv[1] if len(sys.argv)>1 else "shared_all_atlas_v2_mpi4"
    run=ROOT/name;assert run.parent==ROOT
    c=json.loads((run/"contract.json").read_text());exit=int((run/"native_exit_code.txt").read_text())
    for filename in ("domain_source_integrity.txt","domain_output_integrity.txt","final_source_integrity.txt","final_input_integrity.txt"):
        lines=(run/filename).read_text().splitlines();assert lines and all(line.endswith(": OK") for line in lines)
    for filename in ("source.sha256","retained_input.sha256","binary.sha256","output.sha256"):
        for line in (run/filename).read_text().splitlines():
            digest,path=line.split("  ",1);assert sha(Path(path))==digest
    assert sha(run/"caller_unchanged_state.bin")=="74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"
    expected=rows(ROOT/"shared_registry_v2_mpi4/source_target_registry.csv")
    ledger={};counts=defaultdict(Counter);partial=[];new_total=0;rank_bytes=[]
    for rank in range(c["MPI_ranks"]):
        for row in rows(run/f"actual_atlas_cost_rank{rank}.csv"):
            slot=int(row["slot_id"]);assert slot not in ledger and slot%c["MPI_ranks"]==rank
            old=expected[slot]
            assert row["source"]==old["source"] and row["target"]==old["target"] and row["component"]==old["component"]
            complete=row["atlas_complete"]=="1";reused=row["reused"]=="1";panels=int(row["stored_panels"])
            assert reused==(row["target"]=="68526")
            assert row["analytic_ready"]==str(int(reused))
            assert row["planned_positive_queries"]==(str(panels*256) if complete else "UNKNOWN")
            if row["empty"]=="1":assert complete and panels==0
            if reused:assert complete and panels==int(old["retained_panels"])
            if not complete:partial.append(row)
            target=int(row["target"]);counts[target]["incomplete"]+=not complete
            counts[target]["retained_panels" if reused else "new_panels"]+=panels
            counts[target]["chart_windows"]+=int(row["chart_root_windows"])
            ledger[slot]=row
        raw=(run/f"typed_new_atlas_rank{rank}.bin").read_bytes();rank_bytes.append(len(raw))
        h=struct.unpack_from("<8Q",raw);assert h[:4]==(0x475341544c415331,1,1,704) and h[5:]==(rank,24,384)
        assert len(raw)==64+h[4]*704;seen=Counter()
        for index in range(h[4]):
            offset=64+index*704
            slot,source,panel,target,policy,dependent,lower,upper=struct.unpack_from("<4Q4i",raw,offset)
            assert slot%c["MPI_ranks"]==rank and slot//24==source
            assert ledger[slot]["target"]==str(target) and target!=68526
            assert panel==seen[slot];seen[slot]+=1
            chart,axis,lf,uf=struct.unpack_from("<4i",raw,offset+48)
            assert 0<=chart<6 and (lf,uf)==(lower,upper)
            whole=raw[offset+192]!=0
            assert policy==(0 if whole else 1 if dependent==1 else 2)
            if not whole:assert dependent==axis and axis in (0,1)
            apex,part=struct.unpack_from("<2i",raw,offset+368)
            assert 2*apex+part==chart
            def q(at):return decode_quad(int.from_bytes(raw[offset+at:offset+at+16],"little"))
            assert q(688)>=0 and q(672)>q(688) # positive reference determinant ball
            if whole:
                ch,ir,it,nr,nt=struct.unpack_from("<5i",raw,offset+208)
                assert ch==chart and nr>0 and nt>0 and 0<=ir<nr and 0<=it<nt
            else:
                assert q(64)<q(96) and q(128)<q(160)
                for at in (80,112,144,176):assert q(at)>=0
        for slot,count in seen.items():assert count==int(ledger[slot]["stored_panels"])
        assert sum(seen.values())==h[4]
        assert h[4]==sum(int(r["stored_panels"]) for slot,r in ledger.items() if slot%c["MPI_ranks"]==rank and r["reused"]=="0")
        new_total+=h[4]
    assert set(ledger)==set(range(9216))
    cost=rows(run/"all_target_cost.csv");assert len(cost)==24
    old_total=over=0
    for row in cost:
        target=int(row["target"]);v=counts[target]
        assert int(row["incomplete_atlas_domains"])==v["incomplete"]
        assert int(row["new_stored_panels"])==v["new_panels"] and int(row["retained_panels"])==v["retained_panels"]
        assert int(row["chart_root_windows"])==v["chart_windows"] and row["integral_ready"]=="0"
        old_total+=v["retained_panels"]
        planned=(v["new_panels"]+v["retained_panels"])*256
        valid=not v["incomplete"] and planned<=817152
        assert row["original_per_target_query_limit"]=="817152" and row["cost_within_original_limit"]==str(int(valid))
        assert row["complete_planned_queries"]==("UNKNOWN" if v["incomplete"] else str(planned))
        assert row["new_R80_query_packet_bytes"]==("UNKNOWN" if v["incomplete"] else str(v["new_panels"]*256*1008))
        over+=not v["incomplete"] and not valid
    assert old_total==c["retained_original_target_panels"]==1047
    assert new_total==c["new_stored_panels"] and len(partial)==c["incomplete_atlas_domains"]
    assert over==c["complete_targets_exceeding_original_query_budget"]
    assert exit==(2 if partial or over else 0)
    assert c["new_physical_quadrature_queries"]==0 and not c["full_quadrature_allowed"] and not c["all_target_analytic_certificates_complete"]
    print(json.dumps({"scope":"Saved typed atlas ABI/identity and exact all9216 cost accounting; no independent geometric coverage or analytic field evaluation",
                      "run":str(run),"native_exit_code":exit,"contract":c,"packet_and_cost_accounting_pass":True,
                      "known1047_panels_268032_queries_reused":True,"per_rank_stored_typed_plan_bytes":rank_bytes,
                      "actual_incomplete_atlas_domains":partial,"per_target_actual_cost":cost,
                      "actual_atlas_complete":not partial,"all_costs_within_original_target_limit":not(partial or over),
                      "analytic_certificates_complete":False,"full_quadrature_allowed":False},indent=2))


if __name__=="__main__":main()

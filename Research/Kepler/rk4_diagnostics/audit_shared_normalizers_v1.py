"""Saved positive subset bounds/tau accounting and scalar heat oracle."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_regular_q_transfer_v1 import Q, rows, sha

ROOT=Path(__file__).resolve().parent


def main():
    name=sys.argv[1] if len(sys.argv)>1 else "shared_normalizers_v1_mpi4"
    run=ROOT/name;assert run.parent==ROOT
    c=json.loads((run/"contract.json").read_text());exit=int((run/"native_exit_code.txt").read_text())
    for filename in ("published_source_integrity.txt","published_output_integrity.txt","final_source_integrity.txt","final_input_integrity.txt"):
        lines=(run/filename).read_text().splitlines();assert lines and all(line.endswith(": OK") for line in lines)
    for filename in ("source.sha256","retained_input.sha256","binary.sha256","output.sha256"):
        for line in (run/filename).read_text().splitlines():
            digest,path=line.split("  ",1);assert sha(Path(path))==digest
    state="74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"
    assert sha(run/"caller_unchanged_state.bin")==sha(run/"immutable_heat_generating_state.bin")==state
    packet=(ROOT/"material_p3_finite_v2/mpi1_cpu/curved_map_control.bin").read_bytes()
    parameters=[F.from_float(x) for x in struct.unpack_from("<12d",packet,64)]
    assert (len(packet)-160)%384==0;width=(len(packet)-160)//384
    source=rows(run/"whole_source_normalizer_bounds.csv");assert len(source)==384
    registry=rows(ROOT/"shared_registry_v2_mpi4/source_target_registry.csv")
    images=rows(ROOT/"shared_registry_v2_mpi4/source_image_registry.csv")
    scalar_checks=0
    with localcontext() as context:
        context.prec=120
        def decimal(x):return Decimal(x.numerator)/Decimal(x.denominator)
        for s,row in enumerate(source):
            assert int(row["source"])==s
            if row["complete"]!="1":continue
            coords=[F.from_float(x) for x in struct.unpack_from("<6d",packet,160+s*width+24)]
            exact_area=((coords[2]-coords[0])*(coords[5]-coords[1])-(coords[3]-coords[1])*(coords[4]-coords[0]))/2
            area,L,jlo,jhi=(Q(row[k]) for k in ("reference_area_lower","L","Jlower","Jupper"))
            assert 0<area<=exact_area and L==Q(images[s]["L"]) and 0<=L<1
            assert 0<jlo<=(1-L)**2 and jhi>=(1+L)**2
            u=Q(row["frozen_reference_U_lower"])
            # Independent 120-digit scalar corroboration. Exact accounting
            # and retained native elementary-function model remain explicit.
            oracle=decimal(parameters[2]/(parameters[1]-1))*decimal(jhi)**decimal(1-parameters[1])
            assert decimal(u)>0 and decimal(u)<oracle
            scalar_checks+=1
            assert 0<Q(row["mass_subset_lower"])<=area*parameters[6]
            assert 0<Q(row["U_subset_lower"])<=area*u
            assert 0<Q(row["volume_subset_lower"])<=area*jlo
            assert Q(row["fixed_centroid_speed_lower"])>0 and Q(row["max_reference_vertex_radius_lower"])>0
    tau={(int(r["target"]),r["moment"]):r for r in rows(run/"all12_physical_tau.csv")}
    assert len(tau)==288
    original={r["moment"]:r for r in rows(ROOT/"shared_missing_source_v1_mpi1/all12_ORIGINAL_whole_target_budget.csv")}
    manifest=rows(run/"actual_target_normalizers.csv");assert len(manifest)==24
    pending=[];subsets=0
    for t,row in enumerate(manifest):
        target=int(row["target"])
        assert int(registry[t]["target"])==target
        if target==68526:
            assert row["scope"]=="ORIGINAL_normalizers_and_tau_unchanged"
            for moment,old in original.items():assert Q(tau[target,moment]["physical_tau"])==Q(old["ORIGINAL_tau"])
            continue
        members=[s for s in range(384) if registry[s*24+t]["coarse_classification"]=="1"]
        assert row["whole_subset_sources"]=="".join(str(s)+";" for s in members)
        possible=[s for s in range(384) if registry[s*24+t]["coarse_classification"]!="0"]
        assert int(row["coarse_nonexcluded_sources"])==len(possible)
        subsets+=len(members)
        if row["complete"]!="1":pending.append(row);continue
        assert members and all(source[s]["complete"]=="1" for s in possible)
        M,U,V=(Q(row[k]) for k in ("mass_lower","frozen_reference_U_lower","current_volume_lower"))
        assert 0<M<=sum(Q(source[s]["mass_subset_lower"]) for s in members)
        assert 0<U<=sum(Q(source[s]["U_subset_lower"]) for s in members)
        assert 0<V<=sum(Q(source[s]["volume_subset_lower"]) for s in members)
        v,r=(Q(row[k]) for k in ("fixed_anchor_vscale_lower","fixed_reference_rscale_lower"))
        assert v==min(Q(source[s]["fixed_centroid_speed_lower"]) for s in possible)
        assert r==min(Q(source[s]["max_reference_vertex_radius_lower"]) for s in possible)
        mass=Q(tau[target,"mass"]["physical_tau"])
        upper={"mass":M/10**12,"px":mass*v,"py":mass*v,"angular":mass*r*v,"kinetic":U/10**10,
               "physical_U":U/10**10,"volume":V/10**12,"qx":mass*r,"qy":mass*r,"polar":mass*r*r,
               "physical_Udot":U/10**9,"potential":U/10**10}
        for moment,limit in upper.items():
            entry=tau[target,moment];assert entry["complete"]=="1" and 0<Q(entry["physical_tau"])<=limit
        assert row["full_power_cell_coverage_claimed"]=="0"
    assert len(pending)==c["normalizer_pending_targets"] and subsets==c["positive_whole_source_subsets"]
    assert exit==(2 if pending else 0) and c["other23_cold_normalizers_complete"]==(not pending)
    assert c["new_physical_quadrature_queries"]==0 and not c["full_quadrature_allowed"] and not c["general_evolved_heat_supported"]
    print(json.dumps({"scope":"Exact saved subset/J/mass/heat/tau accounting; independent 120-digit scalar frozen-heat corroboration, not independent velocity-field or libm certification",
                      "run":str(run),"native_exit_code":exit,"contract":c,"saved_positive_subset_and_tau_accounting_pass":True,
                      "independent_scalar_heat_checks":scalar_checks,"original_target68526_tau_exactly_retained":True,
                      "pending_normalizers":pending,"other23_cold_normalizers_complete":not pending,
                      "general_evolved_heat_supported":False,"full_quadrature_allowed":False},indent=2))


if __name__=="__main__":main()

#!/usr/bin/env python3
"""Small native operator-result audit; no hydrodynamics computation."""
import argparse
import csv
import json
from pathlib import Path
import struct


def binary(path):
    data=path.read_bytes()
    magic,nb,nr=struct.unpack_from("=QQQ",data)
    assert magic==0x4746503348594231
    assert len(data)==24+8*(6*nb+2*nr)
    return data


def face_difference(pa,pb):
    with pa.open() as f: a=list(csv.DictReader(f))
    with pb.open() as f: b=list(csv.DictReader(f))
    assert len(a)==len(b)
    assert all((x["face"],x["q"])==(y["face"],y["q"]) for x,y in zip(a,b))
    return {key:max(abs(float(x[key])-float(y[key])) for x,y in zip(a,b)) for key in ("tx","ty","ax","ay")}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("directories",nargs="+")
    p.add_argument("--output",required=True)
    args=p.parse_args()
    destination=Path(args.output)
    assert not destination.exists()
    report={"scope":"initial fixed-P3 engineering gate; passing does not certify curved remap, shock entropy, scalable native runtime, or t6-10/paper stability", "versions":{}, "integration_comparison":{}}
    docs=[]
    for directory in args.directories:
        root=Path(directory)
        assert (root/"status.txt").read_text().strip()=="AUDIT_COMPLETE"
        a=json.loads((root/"mpi1/contract.json").read_text())
        b=json.loads((root/"mpi4/contract.json").read_text())
        assert [l["triangles"] for l in a["levels"]]==[96,384]
        comparisons={}
        for la,lb in zip(a["levels"],b["levels"]):
            level=la["level"]
            assert la["physical_cells"]==lb["physical_cells"]==24
            assert la["Schur_min_eigenvalue"]>0
            assert la["conforming_scalar_dofs"]==la["broken_scalar_dofs"]-la["trace_rows"]
            for ca,cb in zip(la["cases"],lb["cases"]):
                name=ca["case"]
                assert name==cb["case"]
                prefix=f"level{level}_{name}"
                sa=binary(root/"mpi1"/(prefix+"_state.bin"))
                sb=binary(root/"mpi4"/(prefix+"_state.bin"))
                comparisons[prefix]={"state_bitmatch":sa==sb,"sampled_face_max_absolute_differences":face_difference(root/"mpi1"/(prefix+"_traction.csv"),root/"mpi4"/(prefix+"_traction.csv"))}
        gate={}
        for c0,c1 in zip(a["levels"][0]["cases"],a["levels"][1]["cases"]):
            ratios={key:c0[key]/c1[key] if c1[key] else None for key in ("curvature_error_max","max_cell_curvature_error_L1")}
            per_cell=[]
            for x,y in zip(c0["cells"],c1["cells"]):
                assert x["id"]==y["id"]
                per_cell.append({"id":x["id"],"maximum_reduction":x["pressure_curvature_error_max"]/y["pressure_curvature_error_max"] if y["pressure_curvature_error_max"] else None,"L1_reduction":x["pressure_curvature_error_L1_volume"]/y["pressure_curvature_error_L1_volume"] if y["pressure_curvature_error_L1_volume"] else None})
            gate[c0["case"]]={"final_chi_max":c1["chi_max"],"final_max_cell_chi_L1":c1["max_cell_chi_L1"],"chi_max_bound_met":c1["chi_max"]<=2e-3,"every_cell_chi_L1_bound_met":c1["max_cell_chi_L1"]<=1e-3,"global_reduction_factors":ratios,"global_fourfold_reduction_met":all(v is not None and v>=4 for v in ratios.values()),"cell_reductions_before_numerical_uncertainty_exclusion":per_cell}
        report["versions"][root.name]={"mpi1":a,"mpi4":b,"MPI_comparison":comparisons,"engineering_gate":gate}
        docs.append((root,a))
    for (ra,a),(rb,b) in zip(docs,docs[1:]):
        differences={}
        for la,lb in zip(a["levels"],b["levels"]):
            for ca,cb in zip(la["cases"],lb["cases"]):
                prefix=f"level{la['level']}_{ca['case']}"
                differences[prefix]={"face_max_absolute_difference":face_difference(ra/"mpi1"/(prefix+"_traction.csv"),rb/"mpi1"/(prefix+"_traction.csv")),"metric_absolute_differences":{key:abs(ca[key]-cb[key]) for key in ("chi_max","max_cell_chi_L1","curvature_error_max","max_cell_curvature_error_L1","total_mass","total_U","total_K")}}
        report["integration_comparison"][f"{ra.name}_to_{rb.name}"]=differences
    destination.write_text(json.dumps(report,indent=2)+"\n")
    for name,r in report["versions"].items():
        for case in ("Kepler","Kepler_perturb"):
            g=r["engineering_gate"][case]
            print(name,case,{k:g[k] for k in ("final_chi_max","final_max_cell_chi_L1","chi_max_bound_met","every_cell_chi_L1_bound_met","global_reduction_factors","global_fourfold_reduction_met")})
    print(destination)


if __name__=="__main__":
    main()

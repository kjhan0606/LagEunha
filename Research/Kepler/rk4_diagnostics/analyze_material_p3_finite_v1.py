#!/usr/bin/env python3
"""Audit native finite-state controls. No time integration or hydro in Python."""
import argparse
import csv
import json
from pathlib import Path
import struct


def state(path):
    raw=path.read_bytes()
    magic,n,nr=struct.unpack_from("=QQQ",raw)
    assert magic==0x4746503346495331
    assert len(raw)==24+8*(6*n+2*nr)
    values=struct.unpack_from(f"={6*n+2*nr}d",raw,24)
    return raw,{k:values[i*n:(i+1)*n] for i,k in enumerate(("x","y","vx","vy","ax","ay"))}


def force(path):
    raw=path.read_bytes();magic,n=struct.unpack_from("=QQ",raw)
    assert magic==0x47465033464f5231 and len(raw)==16+16*n
    return raw,struct.unpack_from(f"={2*n}d",raw,16)


def compare(a,b,cases):
    result={}
    for name in cases:
        sa,va=state(a/f"{name}_state.bin");sb,vb=state(b/f"{name}_state.bin")
        fa,fva=force(a/f"{name}_force.bin");fb,fvb=force(b/f"{name}_force.bin")
        with (a/f"level1_{name}_traction.csv").open() as f:ta=list(csv.DictReader(f))
        with (b/f"level1_{name}_traction.csv").open() as f:tb=list(csv.DictReader(f))
        assert len(ta)==len(tb)
        assert all((x["face"],x["q"])==(y["face"],y["q"]) for x,y in zip(ta,tb))
        result[name]={"state_bitmatch":sa==sb,"force_bitmatch":fa==fb,"force_max_relative_difference":max(abs(x-y) for x,y in zip(fva,fvb))/max(abs(x) for x in fva),"coefficient_max_absolute_differences":{k:max(abs(x-y) for x,y in zip(va[k],vb[k])) for k in va},"physical_face_max_absolute_differences":{k:max(abs(float(x[k])-float(y[k])) for x,y in zip(ta,tb)) for k in ("tx","ty","ax","ay")}}
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument("directories",nargs="+");p.add_argument("--output",required=True);args=p.parse_args();destination=Path(args.output);assert not destination.exists()
    report={"scope":"finite-state P3 material-map RHS, stress-work, actual internal angular moments, CPU/OpenMP/MPI/CUDA; does not establish a timestep, remap, shock or long Kepler result","versions":{},"integration_sensitivity":{}}
    documents=[]
    for directory in args.directories:
        root=Path(directory);assert (root/"status.txt").read_text().strip()=="AUDIT_COMPLETE"
        docs={k:json.loads((root/k/"contract.json").read_text()) for k in ("mpi1_cpu","mpi4_cpu","mpi1_cuda","mpi4_cuda")}
        cases=[c["case"] for c in docs["mpi1_cpu"]["cases"]];assert len(cases)==6
        for mode,doc in docs.items():
            assert doc["physical_cells"]==24 and doc["triangles"]==384 and doc["OMP_threads"]==2
            assert [c["case"] for c in doc["cases"]]==cases
            for c in doc["cases"]:
                assert c["invalid_points"]==0 and c["sampled_min_J"]>0 and c["sampled_min_P"]>0
                assert c["point_energy_identity_relative"]<1e-11 and c["point_angular_identity_relative"]<1e-11
                assert c["max_energy_rate_residual_over_U"]<1e-8
                assert c["max_angular_rate_residual_over_2K"]<1e-11 and c["max_spin_rate_residual_over_2K"]<1e-11
                assert c["centered_matrix_moment_relative_error"]<1e-9
                assert c["mass_moment_relative_error"]<1e-11 and c["mass_matrix_kinetic_error_over_U"]<1e-8
                assert c["volume_rate_identity_over_area"]<1e-8
                if c["rigid_control"]:
                    assert c["thermal_rate_over_U_max"]<1e-9 and c["rigid_remainder_over_internalK_max"]<1e-10
            base=state(root/mode/"curved_Kepler_state.bin")[1]
            boost=state(root/mode/"curved_boost_state.bin")[1]
            assert force(root/mode/"curved_Kepler_force.bin")[0]==force(root/mode/"curved_boost_force.bin")[0]
            assert base["ax"]==boost["ax"] and base["ay"]==boost["ay"]
            initial=state(root/mode/"initial_Kepler_state.bin")[1]
            assert max(abs(x-y) for x,y in zip(base["x"],initial["x"]))>1e-5
            assert max(abs(x-y) for x,y in zip(base["vx"],initial["vx"]))>1e-3
            packet=(root/mode/"curved_map_control.bin").read_bytes();header=struct.unpack_from("=8Q",packet)
            assert header[0]==0x4746503346494e31 and header[1:4]==(1,384,24) and header[7]==864
            assert len(packet)==64+96+384*864
        cpu=compare(root/"mpi1_cpu",root/"mpi4_cpu",cases);gpu=compare(root/"mpi1_cuda",root/"mpi4_cuda",cases)
        assert all(c["state_bitmatch"] and c["force_bitmatch"] for c in list(cpu.values())+list(gpu.values()))
        cg=compare(root/"mpi1_cpu",root/"mpi1_cuda",cases)
        assert all(c["force_max_relative_difference"]<1e-11 for c in cg.values())
        packets=[(root/mode/"curved_map_control.bin").read_bytes() for mode in docs]
        assert all(x==packets[0] for x in packets)
        report["versions"][root.name]={"native":docs,"MPI_cpu":cpu,"MPI_cuda":gpu,"CPU_CUDA":cg,"curved_control_packet_all_backends_bitmatch":True,"boost_force_and_acceleration_bitmatch":True,"finite_operator_engineering_gates":"pass"}
        documents.append((root,docs))
    for (ra,da),(rb,db) in zip(documents,documents[1:]):
        fields=compare(ra/"mpi1_cpu",rb/"mpi1_cpu",cases)
        for a,b in zip(da["mpi1_cpu"]["cases"],db["mpi1_cpu"]["cases"]):
            fields[a["case"]]["cell_integral_max_absolute_differences"]={key:max(abs(x[key]-y[key]) for x,y in zip(a["cells"],b["cells"])) for key in ("mass","area","U","K","intrinsic_angular","polar_inertia","internal_kinetic","rigid_rotation_kinetic")}
        report["integration_sensitivity"][f"{ra.name}_to_{rb.name}"]=fields
    destination.write_text(json.dumps(report,indent=2)+"\n")
    for name,v in report["versions"].items():
        print(name,v["finite_operator_engineering_gates"],"CPU/CUDA force relative max",max(x["force_max_relative_difference"] for x in v["CPU_CUDA"].values()))
    print(destination)


if __name__=="__main__":main()

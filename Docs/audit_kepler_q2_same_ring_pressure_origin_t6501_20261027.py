"""Test whether torque-driving same-ring pressure jumps existed in the q2 IC."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_faceforce_pair_t6501_20261027 import rows
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
RUN = ROOT / "gpu_energy_neutral_capacity64_v5_first_ring_faces_t6501_ada_v1"
OUTPUT = Path(__file__).resolve().parent / "KEPLER_Q2_SAME_RING_PRESSURE_ORIGIN_T6501_2026-10-27.json"


def main():
    audit = json.loads((RUN / "audit_t6501.json").read_text())
    manifest = json.loads((RUN / "stage_manifest.json").read_text())
    if not all(audit["checks"].values()):
        raise ValueError("Native face replay is not accepted")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256")):
        if digest(Path(audit[key])) != audit[sha_key]:
            raise ValueError(f"Changed {key} snapshot")
    initial, _ = ordered(Path(audit["initial"]), audit["n"])
    start, _ = ordered(Path(audit["parent"]), audit["n"])
    r0 = np.hypot(initial["x"].astype(float)-3.,
                  initial["y"].astype(float)-3.)
    inputs = {int(row["index"]): row for row in rows(RUN / "face_input_r0_call1.tsv")}
    faces = rows(RUN / "face_dump_r0_call1.tsv")
    selected = manifest["selected_ids"]
    initial_force = {i: np.zeros(2) for i in selected}
    evolved_force = {i: np.zeros(2) for i in selected}
    jump0, jump6, radius_gap = [], [], []
    for face in faces:
        i, j = int(face["id_i"]), int(face["id_j"])
        if j < 0 or abs(r0[i]-r0[j]) >= 1e-9:
            continue
        a = inputs[int(face["i"])]
        b = inputs[int(face["j"])]
        xi = np.array((float(a["x"]), float(a["y"])))
        xj = np.array((float(b["x"]), float(b["y"])))
        x1 = np.array((float(face["x1"]), float(face["y1"])))
        x2 = np.array((float(face["x2"]), float(face["y2"])))
        ds = np.array((float(face["dSx"]), float(face["dSy"])))
        distance = np.linalg.norm(xj-xi)
        coeff = .5*ds+np.linalg.norm(ds)*(.5*(x1+x2-xi-xj))/distance
        p0 = float(initial["pressure"][i])-float(initial["pressure"][j])
        p6 = float(a["P"])-float(b["P"])
        initial_force[i] += p0*coeff
        evolved_force[i] += p6*coeff
        jump0.append(p0/float(initial["pressure"][i]))
        jump6.append(p6/float(initial["pressure"][i]))
        radius_gap.append(abs(r0[i]-r0[j]))
    ids = np.asarray(selected, dtype=int)
    x = start["x"][ids].astype(float)-3.
    y = start["y"][ids].astype(float)-3.
    mass = start["mass"][ids].astype(float)
    weights = mass/mass.sum()
    torques = {}
    for key, forces in (("initial", initial_force), ("evolved", evolved_force)):
        accel = np.array([forces[int(i)]/mass[k] for k, i in enumerate(ids)])
        torque = x*accel[:, 1]-y*accel[:, 0]
        torques[key] = float(np.sqrt(np.sum(weights*torque*torque)))
    jump0 = np.asarray(jump0)
    jump6 = np.asarray(jump6)
    if len(jump0) != 128:
        raise ValueError("Expected 128 native directed same-ring faces")
    report = {"face_run": str(RUN), "parent_sha256": audit["parent_sha256"],
              "face_dump_sha256": digest(RUN / "face_dump_r0_call1.tsv"),
              "same_ring_directed_faces": len(jump0),
              "max_initial_radius_gap": max(radius_gap),
              "initial_pressure_jump_relative_mean_abs": float(np.mean(np.abs(jump0))),
              "initial_pressure_jump_relative_max_abs": float(np.max(np.abs(jump0))),
              "evolved_pressure_jump_relative_mean_abs": float(np.mean(np.abs(jump6))),
              "evolved_pressure_jump_relative_p99_abs": float(np.quantile(np.abs(jump6), .99)),
              "same_ring_orbital_torque_rms_if_initial_pressure": torques["initial"],
              "same_ring_orbital_torque_rms_evolved_pressure": torques["evolved"],
              "scope": "Identical t=6.5 native face geometry, replacing only "
                       "same-ring pressure jumps by the t=0 IC values. "
                       "Selected 64 first-ring cells; no time evolution."}
    with OUTPUT.open("x") as out:
        json.dump(report, out, indent=2, allow_nan=False)
        out.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

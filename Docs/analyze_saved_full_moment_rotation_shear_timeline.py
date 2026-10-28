"""Read-only affine-rotation and shear-energy audit of saved Kepler states.

Uses the accepted Voronoi initial condition and existing coarse/fine
checkpoints at common t=.025 and t=.05. No RHS, provider, or time step runs.
Krot is the rigid-spin projection inside affine kinetic energy; the remainder
is the non-rigid affine deformation projection, not an extra energy term.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from compare_full_moment_checkpoint_errors import load


REGIONS = ("all", "cavity", "inner_contact", "disk_core")
FIELDS = ("ell", "Krot", "Kaffine", "Kaffine_nonspin")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def decode(path, physics_folder=None):
    path = Path(path)
    with np.load(path) as z:
        ids = z["ids"].copy()
        mass = z["mass"].copy()
        state = z["state"].copy()
        time = float(z["time"])
    n = len(ids)
    if state.size != 20 * n + 4 or mass.shape != (n,):
        raise SystemExit(f"unexpected native checkpoint layout: {path}")
    cell = state[:-4].reshape(n, 20)
    raw = np.column_stack((mass, cell[:, 3:])).astype(np.longdouble)
    m = raw[:, 0]
    q = raw[:, 1:3]
    Q = raw[:, 3:7].reshape(n, 2, 2)
    p = raw[:, 7:9]
    D = raw[:, 9:13].reshape(n, 2, 2)
    outer = lambda a, b: np.einsum("ni,nj->nij", a, b)
    S = Q - outer(q, q) / m[:, None, None]
    S = .5 * (S + S.transpose(0, 2, 1))
    C = D - outer(q, p) / m[:, None, None]
    det = S[:, 0, 0] * S[:, 1, 1] - S[:, 0, 1] * S[:, 1, 0]
    if np.any(m <= 0) or np.any(det <= 0):
        raise SystemExit(f"non-positive cell mass/covariance: {path}")
    invS = np.empty_like(S)
    invS[:, 0, 0] = S[:, 1, 1] / det
    invS[:, 1, 1] = S[:, 0, 0] / det
    invS[:, 0, 1] = -S[:, 0, 1] / det
    invS[:, 1, 0] = -S[:, 1, 0] / det
    B = np.einsum("nji,njk->nik", C, invS)
    Kaffine = .5 * np.einsum("nij,njk,nik->n", B, S, B)
    ell = C[:, 0, 1] - C[:, 1, 0]
    inertia = S[:, 0, 0] + S[:, 1, 1]
    Krot = ell * ell / (2 * inertia)
    result = dict(path=str(path.resolve()), sha256=sha256(path), time=time,
                  ids=ids, mass=mass, weight=cell[:, 2].copy(), ell=ell, Krot=Krot,
                  Kaffine=Kaffine, Kaffine_nonspin=Kaffine-Krot)
    if physics_folder is not None:
        fields, report, source = load(physics_folder)
        if report["time"] != time or not np.array_equal(fields["ids"], ids):
            raise SystemExit(f"checkpoint/physical-field identity mismatch: {path}")
        if not np.array_equal(fields["actual_mass"], mass):
            raise SystemExit(f"checkpoint/physical-field mass mismatch: {path}")
        result.update(fields={key: fields[key] for key in REGIONS},
                      physical_source=source)
    return result


def weighted_rms(value, mass, mask):
    v = np.asarray(value, dtype=np.longdouble)[mask]
    w = np.asarray(mass, dtype=np.longdouble)[mask]
    return float(np.sqrt(np.sum(w * v * v) / np.sum(w)))


def summarize(state, masks):
    mass = state["mass"]
    rows = {}
    for name, mask in masks.items():
        if not np.any(mask):
            continue
        ka = np.asarray(state["Kaffine"], dtype=np.longdouble)[mask]
        kr = np.asarray(state["Krot"], dtype=np.longdouble)[mask]
        weight = np.asarray(state["weight"], dtype=np.longdouble)[mask]
        rows[name] = dict(
            cells=int(mask.sum()),
            gauge_removed_weight_min=float((weight-weight.mean()).min()),
            gauge_removed_weight_max=float((weight-weight.mean()).max()),
            gauge_removed_weight_RMS=float(np.sqrt(np.mean((weight-weight.mean())**2))),
            internal_angular_momentum_sum=float(np.sum(state["ell"][mask], dtype=np.longdouble)),
            internal_angular_momentum_mass_weighted_RMS=weighted_rms(state["ell"], mass, mask),
            affine_kinetic_energy=float(np.sum(ka, dtype=np.longdouble)),
            rigid_spin_projection_energy=float(np.sum(kr, dtype=np.longdouble)),
            affine_nonspin_projection_energy=float(np.sum(ka-kr, dtype=np.longdouble)),
            rigid_spin_fraction_of_affine_energy=float(np.sum(kr, dtype=np.longdouble)/np.sum(ka, dtype=np.longdouble)),
        )
    return rows


def compare(a, b):
    if a["time"] != b["time"] or not np.array_equal(a["ids"], b["ids"]):
        raise SystemExit("coarse/fine endpoint identity mismatch")
    if not np.array_equal(a["mass"], b["mass"]):
        raise SystemExit("coarse/fine fixed masses differ")
    masks = {}
    for name in REGIONS:
        ma, mb = a["fields"][name], b["fields"][name]
        masks[name] = ma | mb
    mass = a["mass"]
    rows = {}
    for name, mask in masks.items():
        w = np.asarray(mass, dtype=np.longdouble)[mask]
        row = {"cells_union": int(mask.sum()),
               "region_membership_changes": int(np.count_nonzero(
                   a["fields"][name] != b["fields"][name])),
               "fields": {}}
        for field in FIELDS:
            av = np.asarray(a[field], dtype=np.longdouble)[mask]
            bv = np.asarray(b[field], dtype=np.longdouble)[mask]
            diff = bv-av
            rms = float(np.sqrt(np.sum(w*diff*diff)/np.sum(w)))
            ref_rms = float(np.sqrt(np.sum(w*bv*bv)/np.sum(w)))
            row["fields"][field] = dict(
                coarse_sum=float(np.sum(av, dtype=np.longdouble)),
                fine_sum=float(np.sum(bv, dtype=np.longdouble)),
                fine_minus_coarse_sum=float(np.sum(diff, dtype=np.longdouble)),
                mass_weighted_RMS_difference=rms,
                difference_over_fine_RMS=rms/ref_rms if ref_rms > 0 else None,
            )
        rows[name] = row
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(
        "/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019"))
    parser.add_argument("--initial", type=Path, default=Path(
        "/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/initial.npz"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")

    states = {"t0_voronoi": decode(args.initial)}
    source_paths = {"initial": str(args.initial.resolve())}
    for time_tag, branches in (
        ("t0025", (("coarse", "coarse_0005"), ("fine", "fine_0010"))),
        ("t005", (("coarse", "coarse_0010"), ("fine", "fine_0020"))),
    ):
        for branch, checkpoint_tag in branches:
            folder = args.run_root / branch / "checkpoints" / checkpoint_tag
            key = f"{branch}_{time_tag}"
            states[key] = decode(folder / "accepted.npz", folder / "physics")
            source_paths[key] = str((folder / "accepted.npz").resolve())

    initial = states["t0_voronoi"]
    initial_weights_all_zero = bool(np.all(initial["weight"] == 0.))
    if not initial_weights_all_zero:
        raise SystemExit("initial weights are not exactly zero")

    summaries = {"t0_voronoi": summarize(initial, {"all": np.ones(len(initial["ids"]), bool)})}
    comparisons = {}
    for time_tag in ("t0025", "t005"):
        a, b = states[f"coarse_{time_tag}"], states[f"fine_{time_tag}"]
        if a["time"] != b["time"]:
            raise SystemExit(f"coarse/fine time mismatch at {time_tag}")
        if not np.array_equal(a["mass"], initial["mass"]) or not np.array_equal(b["mass"], initial["mass"]):
            raise SystemExit(f"cell mass changed from the Voronoi initial state at {time_tag}")
        summaries[f"coarse_{time_tag}"] = summarize(a, a["fields"])
        summaries[f"fine_{time_tag}"] = summarize(b, b["fields"])
        comparisons[time_tag] = compare(a, b)

    result = dict(
        status="completed_saved_rotation_shear_timeline",
        scope="Saved-state affine/spin diagnostics only; no provider or accepted time step.",
        weight_slot=2,
        initial_weights_all_exactly_zero=initial_weights_all_zero,
        fixed_mass_bitwise_preserved_from_initial=True,
        definitions={
            "S": "Q - q q^T/m",
            "C": "D - q p^T/m",
            "B": "C^T S^-1",
            "Kaffine": "0.5 Tr(B S B^T)",
            "ell_internal": "C_xy - C_yx",
            "Krot": "ell_internal^2/(2 Tr(S)); a projection inside Kaffine",
            "Kaffine_nonspin": "Kaffine - Krot; remaining affine deformation projection",
        },
        summaries=summaries,
        coarse_fine_endpoint_differences=comparisons,
        sources={key: dict(path=path, sha256=sha256(path)) for key, path in source_paths.items()},
        limitations=[
            "Endpoint diagnostics do not establish a relaxation law or causal reduction of Kepler error.",
            "The nonspin affine projection includes deformation/shear and is not added as a new energy.",
            "No long-time t=6-10 stability or 3D validation is performed.",
        ],
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "summaries": summaries,
                      "comparisons": comparisons}, sort_keys=True))


if __name__ == "__main__":
    main()

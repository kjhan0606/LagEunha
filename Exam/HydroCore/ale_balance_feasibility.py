"""Frozen net-mass ALE balance and first-RHS thermal representation audit.

Analytic fields and original native geometry only. No fluid evolution or new
evolved moments. A small direct weight-rate solve is a reference existence
check, not a production solver or cost result. Never remove the RHS mean.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import struct
import subprocess

import numpy as np

GAMMA = float(np.float32(5/3))


def arrays(g, scale=1.):
    d, n = g["dimension"], g["cells"]
    f = g["faces"]
    i, j = np.array([v["i"] for v in f]), np.array([v["j"] for v in f])
    r = scale*np.array([v["r"] for v in f])
    face_area = scale**(d-1)*np.array([v["area"] for v in f])
    qface = np.concatenate([np.full(len(v["quad"]), k) for k, v in enumerate(f)])
    quadrature = np.concatenate([np.array(v["quad"]) for v in f])
    x, V = scale*np.array(g["positions"]), scale**d*np.array(g["volume"])
    relative = scale*quadrature[:, 1:]
    return {"dim": d, "n": n, "i": i, "j": j, "r": r, "area": face_area,
            "distance": np.linalg.norm(r, axis=1), "qface": qface,
            "owner": i[qface], "neighbor": j[qface], "weight": scale**(d-1)*quadrature[:, 0],
            "relative": relative, "point": x[i[qface]]+relative, "x": x, "V": V,
            "centroid": scale*np.array(g["centroid"]),
            "covariance": scale**(d+2)*np.array(g["central_second_moment"]),
            "cut": np.max(abs(r-(x[j]-x[i])), axis=1) > 1e-10*scale}


def aggregate(a, value, boundary=False):
    value = np.asarray(value)
    out = np.zeros((a["n"],)+value.shape[1:])
    weight = a["weight"]
    if boundary:
        weight = weight*a["cut"][a["qface"]]
    np.add.at(out, a["owner"], weight.reshape((-1,)+(1,)*(value.ndim-1))*value)
    return out


def evaluate(a, G, beta, pbase, gradient, rate=None, boost=None):
    d, n = a["dim"], a["n"]
    boost = np.zeros(d) if boost is None else boost
    beta = beta+boost
    qv = a["x"]@G.T+beta
    mean_v = a["centroid"]@G.T+beta
    i, j, face = a["owner"], a["neighbor"], a["qface"]
    normal_face = a["r"]/a["distance"][:, None]
    normal = normal_face[face]
    w0 = np.einsum("ka,ka->k", normal, (qv[i]+qv[j])/2)
    w0 -= np.einsum("ka,ka->k", qv[j]-qv[i], a["relative"]-a["r"][face]/2)/a["distance"][face]
    fluid = a["point"]@G.T+beta
    un = np.einsum("ka,ka->k", fluid, normal)
    rho = 1.
    raw_rhs = aggregate(a, rho*(un-w0))
    L = np.zeros((n, n))
    conductance = rho*a["area"]/(2*a["distance"])
    np.add.at(L, (a["i"], a["i"]), conductance)
    np.add.at(L, (a["i"], a["j"]), -conductance)
    compatibility = float(raw_rhs.sum(dtype=np.longdouble))
    tolerance = 4096*np.finfo(float).eps*max(1., float(np.sum(abs(raw_rhs))))
    assert abs(compatibility) < tolerance, ("Incompatible raw mass RHS; do not subtract its mean", compatibility)
    if rate is None:
        rate = np.zeros(n)
        rate[1:] = np.linalg.solve(L[1:, 1:], raw_rhs[1:])
        rate -= rate.mean()  # Weight gauge only. The RHS remains untouched.
    residual = L@rate-raw_rhs
    assert np.max(abs(residual)) < 1e-11
    wn = w0+(rate[i]-rate[j])/(2*a["distance"][face])
    pressure = pbase+a["point"]@gradient
    assert pressure.min() > 0
    fm = rho*(un-wn)
    fp = fm[:, None]*fluid+pressure[:, None]*normal
    fe = fm*(pressure/((GAMMA-1)*rho)+.5*np.sum(fluid*fluid, axis=1))+pressure*un
    mdot, pdot, edot = -aggregate(a, fm), -aggregate(a, fp), -aggregate(a, fe)
    Vdot = aggregate(a, wn)
    first_dot = aggregate(a, wn[:, None]*a["point"])
    cdot = (first_dot-a["centroid"]*Vdot[:, None])/a["V"][:, None]
    m, p = rho*a["V"], rho*a["V"][:, None]*mean_v
    centroid_identity = m[:, None]*cdot-p+aggregate(a, fm[:, None]*a["point"])
    qdot_identity = mean_v+(a["x"]-a["centroid"])@G.T-qv
    # Incompressible affine Euler initial data, with exact temporal derivative.
    Gdot = -G@G
    betadot = -G@beta-gradient/rho
    expected_pdot = m[:, None]*(cdot@G.T+a["centroid"]@Gdot.T+betadot)
    proxy_Udot = edot-np.sum(mean_v*pdot, axis=1)+.5*np.sum(mean_v*mean_v, axis=1)*mdot
    fluctuation = (a["point"]-a["centroid"][i])@G.T
    variance_dot = aggregate(a, .5*rho*np.sum(fluctuation*fluctuation, axis=1)*wn)
    variance_dot += rho*np.einsum("ab,nbc,ac->n", Gdot, a["covariance"], G)
    Uphysical = (pbase+a["centroid"]@gradient)*a["V"]/(GAMMA-1)
    global_ledgers = {
        "mass": float(abs(mdot.sum()+aggregate(a, fm, True).sum())),
        "momentum": float(np.max(abs(pdot.sum(axis=0)+aggregate(a, fp, True).sum(axis=0)))),
        "energy": float(abs(edot.sum()+aggregate(a, fe, True).sum()))}
    checks = {"mass_residual_max": float(np.max(abs(mdot))),
              "volume_rate_max": float(np.max(abs(Vdot))),
              "raw_weight_solve_residual_max": float(np.max(abs(residual))),
              "centroid_transport_identity_max": float(np.max(abs(centroid_identity))),
              "site_vs_cellmean_velocity_identity_max": float(np.max(abs(qdot_identity))),
              "momentum_Euler_reference_max_error": float(np.max(abs(pdot-expected_pdot))),
              "thermal_proxy_vs_variance_rate_max_error": float(np.max(abs(proxy_Udot-variance_dot))),
              "internal_plus_periodic_cut_ledgers": global_ledgers}
    assert max(v for k, v in checks.items() if not isinstance(v, dict)) < 1e-10, checks
    assert max(global_ledgers.values()) < 1e-10, checks
    summary = {"checks": checks, "raw_rhs_compatibility": compatibility, "rhs_mean_removed": False,
               "weight_gauge": "arithmetic mean zero", "rate_max": float(np.max(abs(rate))),
               "surface_mass_flux_absolute_sum": float(np.sum(a["weight"]*abs(fm))),
               "thermal_proxy_rate_L1_per_volume": float(np.sum(abs(proxy_Udot))/a["V"].sum()),
               "thermal_proxy_rate_signed_sum": float(proxy_Udot.sum()),
               "thermal_proxy_rate_max_per_physical_U": float(np.max(abs(proxy_Udot)/Uphysical)),
               "physical_Udot": "zero for the registered incompressible affine/rest-pressure initial states",
               "periodic_cut_scope": "Physical affine continuation on the owned-cell patch with fixed periodic geometry images. Cut fluxes are explicit external boundary ledgers, not a new production shear-periodic BC."}
    fields = {"rate": rate, "qdot": qv, "mdot": mdot, "pdot": pdot, "edot": edot,
              "Vdot": Vdot, "cdot": cdot, "proxy_Udot": proxy_Udot,
              "variance_dot": variance_dot, "Gdot": Gdot}
    fluxes = {"mass": fm, "momentum": fp, "energy": fe, "wn": wn}
    return summary, fields, fluxes


def finite_difference(executable, directory, a, fields, G, side, env):
    direction = directory/"direction.bin"
    with direction.open("xb") as f:
        f.write(struct.pack("<ii", a["dim"], a["n"]))
        f.write(np.asarray(fields["qdot"], dtype="<f8").tobytes())
        f.write(np.asarray(fields["rate"], dtype="<f8").tobytes())
    rows = []
    for level, epsilon in enumerate([1e-5, 5e-6]):
        sampled = []
        for sign, name in [(1, "plus"), (-1, "minus")]:
            output = directory/f"{name}_{level}.json"
            command = [str(executable), str(a["dim"]), str(side), str(sign*epsilon), str(direction), "1"]
            (directory/f"{name}_{level}.command.json").write_text(json.dumps(command)+"\n")
            with output.open("x") as stdout,(directory/f"{name}_{level}.log").open("x") as stderr:
                code = subprocess.run(command, env=env, cwd=directory, stdout=stdout, stderr=stderr).returncode
            assert code == 0, (output, code)
            sampled.append(json.loads(output.read_text()))
        plus, minus = sampled
        vdot = (np.array(plus["volume"])-np.array(minus["volume"]))/2/epsilon
        cdot = (np.array(plus["centroid"])-np.array(minus["centroid"]))/2/epsilon
        kinetic = []
        for sign, state in [(1, plus), (-1, minus)]:
            At = G+sign*epsilon*fields["Gdot"]
            kinetic.append(.5*np.einsum("ab,nbc,ac->n", At, np.array(state["central_second_moment"]), At))
        kdot = (kinetic[0]-kinetic[1])/2/epsilon
        verror = float(np.max(abs(vdot-fields["Vdot"])))
        cerror = float(np.max(abs(cdot-fields["cdot"])))
        kerror = float(np.max(abs(kdot-fields["variance_dot"])))
        kscale = max(1e-9, float(np.max(abs(fields["variance_dot"]))))
        row = {"epsilon": epsilon, "volume_rate_max_error": verror,
               "centroid_rate_max_error": cerror, "variance_rate_max_error": kerror,
               "variance_rate_relative_error": kerror/kscale}
        rows.append(row)
        assert verror < 1e-7 and cerror < 1e-7 and kerror/kscale < 1e-4, row
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", type=Path)
    args = parser.parse_args()
    out = args.stage
    assert not (out/"result.json").exists()
    registration = {"scope": __doc__, "native_cases": [[1, 16], [2, 3], [3, 3]],
                    "fields": ["translation", "rigid_rotation_Dge2", "steady_affine_shear_Dge2", "rest_linear_pressure"],
                    "rho": 1, "constant_pressure": 1e-8, "linear_pressure_base": 2,
                    "RHS_compatibility": "No RHS mean subtraction. A failed compatibility test stops.",
                    "weight_solver": "Small pinned direct reference only. Mean-zero solution gauge.",
                    "finite_difference_epsilons": [1e-5, 5e-6],
                    "geometry_FD_gate": "Volume/centroid max error <1e-7; variance derivative relative error <1e-4 with absolute scale floor 1e-9 for exact-zero controls. No extra epsilon scan.",
                    "quadrature": "Exact degree three for uniform-density affine flux, from original native faces. 2-point lines; positive 3x3 Duffy face triangles.",
                    "stop_rule": "No scalar m,p,E model adoption if physical thermal energy and unresolved kinetic variance need an undefined correction. Report FV spatial order separately; no PDE, native128, limiter or new evolved reservoir.",
                    "new_PDE_steps": 0}
    (out/"registration.json").write_text(json.dumps(registration, indent=2)+"\n")
    env = dict(os.environ, TMPDIR=str(out/"tmp"), OMP_NUM_THREADS="1")
    results = []
    for dim, side in registration["native_cases"]:
        case = out/f"D{dim}_N{side}"
        g = json.loads((case/"geometry.json").read_text())
        a = arrays(g)
        fields_to_check = [("translation", np.zeros((dim, dim)), np.linspace(.2, .1, dim), 1e-8, np.zeros(dim)),
                           ("linear_pressure", np.zeros((dim, dim)), np.zeros(dim), 2., .2*np.eye(dim)[0])]
        if dim > 1:
            rotation = np.zeros((dim, dim));rotation[0, 1]=-1;rotation[1, 0]=1
            shear = np.zeros((dim, dim));shear[0, 1]=1
            fields_to_check.extend([("rotation", rotation, np.zeros(dim), 1e-8, np.zeros(dim)),
                                    ("steady_shear", shear, np.zeros(dim), 1e-8, np.zeros(dim))])
        for name, G, beta, pbase, gradient in fields_to_check:
            directory = case/name;directory.mkdir()
            summary, fields, flux = evaluate(a, G, beta, pbase, gradient)
            boost = np.array([.31, -.23, .17])[:dim]
            boosted, boosted_fields, boosted_flux = evaluate(a, G, beta, pbase, gradient, fields["rate"], boost)
            boost_checks = {"mass": float(np.max(abs(boosted_flux["mass"]-flux["mass"]))),
                            "momentum": float(np.max(abs(boosted_flux["momentum"]-flux["momentum"]-flux["mass"][:, None]*boost))),
                            "energy": float(np.max(abs(boosted_flux["energy"]-flux["energy"]-flux["momentum"]@boost-.5*np.dot(boost, boost)*flux["mass"]))),
                            "thermal_proxy_rate": float(np.max(abs(boosted_fields["proxy_Udot"]-fields["proxy_Udot"])))}
            assert max(boost_checks.values()) < 1e-10, boost_checks
            if np.max(abs(fields["qdot"])) == 0 and np.max(abs(fields["rate"])) == 0:
                fd = [{"scope": "Zero geometry direction; exact zero derivative, no duplicate native geometry replay."}]
            else:
                fd = finite_difference(out/"geometry", directory, a, fields, G, side, env)
            row = {"dimension": dim, "side": side, "field": name, **summary,
                   "Galilean_checks": boost_checks, "native_geometry_FD": fd, "accepted_steps": 0}
            if name == "steady_shear":
                smaller = arrays(g, .5)
                small_summary, small_fields, _ = evaluate(smaller, G, beta, pbase, gradient, .25*fields["rate"])
                ratio = small_summary["thermal_proxy_rate_L1_per_volume"]/summary["thermal_proxy_rate_L1_per_volume"]
                assert abs(ratio-.25) < 1e-8
                row["half_scale_variance_rate_per_volume_ratio"] = ratio
                row["scale_control_scope"] = "Same local irregular pattern and physical affine gradient at half cell length; not a full global PDE convergence test."
            np.savez_compressed(directory/"fields.npz", **fields)
            (directory/"result.json").write_text(json.dumps(row, indent=2)+"\n")
            results.append(row)
            print(json.dumps({"dimension": dim, "field": name, "thermal_proxy_L1_rate_per_volume": summary["thermal_proxy_rate_L1_per_volume"],
                              "physical_U_rate": 0, "surface_mass_flux_absolute_sum": summary["surface_mass_flux_absolute_sum"]}), flush=True)
    result = {"scope": __doc__, "cases": results, "accepted_steps": 0,
              "production_physics_changed": False, "scalar_momentum_energy_model_adopted": False,
              "new_evolved_moments_or_reservoirs": False, "old_Hopkins_projection_repeated": False,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "interpretation": "Net cell mass balance retains face transport and correct linear pressure, but the mean-kinetic thermal proxy has the unresolved kinetic-variance time derivative. No initial thermal projection or evolved-state correction is used to hide the first-RHS term."}
    (out/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    print("All finite balance/reference checks completed. No PDE released.", flush=True)


if __name__ == "__main__":
    main()

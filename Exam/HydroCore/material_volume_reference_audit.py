"""Independent linear-algebra diagnosis of one frozen native coefficient graph.

No fluid evolution and no production solver replacement. The pinned sparse
factorization supplies a mean-zero reference; spectral values are Ritz estimates
with measured residuals, not certified bounds on every uncomputed eigenvalue.
"""
from pathlib import Path
import argparse
import json
import struct
import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu, eigsh

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
args.output.mkdir()
with args.input.open("rb") as stream:
    magic, n, nf = struct.unpack("<QQQ", stream.read(24))
    assert magic == 0x4D56434F53543131
    cells = np.fromfile(stream, dtype="<f8", count=6 * n).reshape(n, 6)
    faces = np.fromfile(stream, dtype=np.dtype([
        ("i", "<i4"), ("j", "<i4"), ("rx", "<f8"), ("ry", "<f8"),
        ("cx", "<f8"), ("cy", "<f8"), ("area", "<f8")]), count=nf)
    assert stream.read() == b""
i, j = faces["i"], faces["j"]
r = np.column_stack([faces["rx"], faces["ry"]])
c = np.column_stack([faces["cx"], faces["cy"]])
distance = np.linalg.norm(r, axis=1)
coeff = faces["area"] / (2 * distance)
diagonal = np.bincount(i, weights=coeff, minlength=n)
L = sparse.diags(diagonal) - sparse.csr_matrix((coeff, (i, j)), shape=(n, n))
bi = faces["area"][:, None] * c / distance[:, None]
bj = faces["area"][:, None] * r / distance[:, None] - bi
K = np.zeros((n, 2, 2))
for a in range(2):
    for b in range(2):
        np.add.at(K[:, a, b], i, bj[:, a] * r[:, b])
K /= cells[:, 2, None, None]
rhs = cells[:, 2, None] * np.column_stack([
    K[:, 1, 1] - K[:, 0, 0], -K[:, 0, 1] - K[:, 1, 0]]) / np.sqrt(2)
compatibility = rhs.sum(axis=0)
assert np.max(abs(compatibility)) < 4096 * np.finfo(float).eps * max(1, np.sum(abs(rhs)))
start = time.monotonic()
factor = splu(L[1:, 1:].tocsc())
psi = np.zeros((n, 2))
psi[1:] = factor.solve(rhs[1:])
psi -= psi.mean(axis=0)
reference_seconds = time.monotonic() - start
residual = L @ psi - rhs
relative = np.linalg.norm(residual, axis=0) / np.linalg.norm(rhs, axis=0)
assert np.max(relative) < 1e-10
np.save(args.output / "psi_reference.npy", psi)
np.save(args.output / "rhs.npy", rhs)
sparse.save_npz(args.output / "native_L.npz", L)

def spectrum(matrix):
    values, vectors = eigsh(matrix, k=6, sigma=-1e-6, which="LM", tol=1e-10, maxiter=256)
    order = np.argsort(values)
    values, vectors = values[order], vectors[:, order]
    errors = np.linalg.norm(matrix @ vectors - vectors * values, axis=0)
    largest, direction = eigsh(matrix, k=1, which="LM", tol=1e-10, maxiter=256)
    largest_error = float(np.linalg.norm(matrix @ direction[:, 0] - largest[0] * direction[:, 0]))
    return {"small_ritz_values": values.tolist(), "ritz_absolute_residuals": errors.tolist(),
            "largest_ritz_value": float(largest[0]), "largest_ritz_residual": largest_error,
            "estimated_positive_condition_ratio": float(largest[0] / values[1]),
            "certified_global_spectral_bound": False}

asym = L - L.T
symmetric = (L + L.T) * 0.5
inverse_root = sparse.diags(1 / np.sqrt(diagonal))
result = {"scope": __doc__, "cells": n, "directed_faces": nf,
          "rhs_compatibility": compatibility.tolist(),
          "matrix_max_asymmetry": float(np.max(abs(asym.data), initial=0)),
          "reference_relative_residuals": relative.tolist(),
          "psi_arithmetic_means": psi.mean(axis=0).tolist(),
          "psi_L2_norms": np.linalg.norm(psi, axis=0).tolist(),
          "reference_factor_and_solve_seconds": reference_seconds,
          "reference_timing_is_production_cost": False,
          "unscaled": spectrum(symmetric),
          "Jacobi_scaled": spectrum(inverse_root @ symmetric @ inverse_root),
          "accepted_steps": 0}
(args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2), flush=True)

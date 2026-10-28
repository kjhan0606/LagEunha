"""Apply the opt-in GPU polygon-moment diagnostic to an isolated worktree."""

import hashlib
import json
import os
import shutil
import sys
from pathlib import Path


def replace_once(text: str, before: str, after: str, label: str) -> str:
    count = text.count(before)
    if count != 1:
        raise RuntimeError(f"{label}: expected one source match, found {count}")
    return text.replace(before, after, 1)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: apply_kepler_polygon_moment_candidate_20261028.py WORKTREE")
    root = Path(sys.argv[1]).resolve()
    exam = root / "Exam"
    header = exam / "exam_gpu.h"
    cuda = exam / "exam_gpu.cu"
    ref_source = Path(
        "/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/"
        "gfs_polygon_moments_reference_20261028.h")
    ref_target = exam / "gfs_polygon_moments_reference_20261028.h"
    backup = root / "../Kepler_spin_moments_20261028_prepatch"

    if not (root / ".git").exists() or not header.is_file() or not cuda.is_file():
        raise RuntimeError(f"Not the expected isolated source worktree: {root}")
    if backup.exists() or ref_target.exists():
        raise RuntimeError("Refusing to overwrite a backup or reference header")
    original_h = header.read_text()
    original_cu = cuda.read_text()
    if "moment_valid" in original_h or "gfs_audit_polygon_moments" in original_cu:
        raise RuntimeError("Moment diagnostic already appears in the source")

    modified_h = replace_once(
        original_h,
        "    int    stencil_pad;                   /* adaptive linked-cell search half-width */\n"
        "    int    halo_clipped;                  /* final influence bound crosses the local grid */\n",
        "    int    stencil_pad;                   /* adaptive linked-cell search half-width */\n"
        "    int    halo_clipped;                  /* final influence bound crosses the local grid */\n"
        "    double moment_area;                    /* raw signed polygon area */\n"
        "    double moment_cx, moment_cy;           /* centroid relative to generator */\n"
        "    double moment_qxx_mass, moment_qxy_mass, moment_qyy_mass;\n"
        "    int    moment_valid;                    /* exact raw polygon and positive mass Q */\n",
        "PolygonOut fields")

    modified_cu = replace_once(
        original_cu,
        "#include \"gfs_pair.h\"\n#include \"gfs_velocity_limiter.h\"\n#include \"gfs_thermal64.h\"\n",
        "#include \"gfs_pair.h\"\n#include \"gfs_velocity_limiter.h\"\n"
        "#include \"gfs_thermal64.h\"\n"
        "#include \"gfs_polygon_moments_reference_20261028.h\"\n",
        "reference include")
    modified_cu = replace_once(
        modified_cu,
        "static int h_geometry_trial=0,h_geometry_failed=0;\n",
        "static int h_geometry_trial=0,h_geometry_failed=0;\n"
        "static int h_gfs_polygon_moment_audit_done=0;\n",
        "one-shot audit flag")
    modified_cu = replace_once(
        modified_cu,
        "    int i = blockIdx.x * blockDim.x + threadIdx.x;\n"
        "    if (i >= nbp) return;\n\n"
        "    double xi = px[i], yi = py[i];\n",
        "    int i = blockIdx.x * blockDim.x + threadIdx.x;\n"
        "    if (i >= nbp) return;\n"
        "    PolygonOut *moment_out = &polygons[i];\n"
        "    moment_out->moment_area = 0.0;\n"
        "    moment_out->moment_cx = 0.0;\n"
        "    moment_out->moment_cy = 0.0;\n"
        "    moment_out->moment_qxx_mass = 0.0;\n"
        "    moment_out->moment_qxy_mass = 0.0;\n"
        "    moment_out->moment_qyy_mass = 0.0;\n"
        "    moment_out->moment_valid = 0;\n\n"
        "    double xi = px[i], yi = py[i];\n",
        "per-cell safe defaults")
    modified_cu = replace_once(
        modified_cu,
        "    /* ---- Compute volume via shoelace formula ---- */\n"
        "    double area = 0;\n"
        "    for (int e = 0; e < nc; e++) {\n"
        "        int en = (e + 1) % nc;\n"
        "        area += poly_x[e] * poly_y[en] - poly_y[e] * poly_x[en];\n"
        "    }\n"
        "    double volume = 0.5 * fabs(area);\n",
        "    /* ---- Compute volume and signed polygon moments in one pass ---- */\n"
        "    double area = 0;\n"
        "    double cx_sum = 0.0, cy_sum = 0.0;\n"
        "    double jxx12 = 0.0, jyy12 = 0.0, jxy24 = 0.0;\n"
        "    for (int e = 0; e < nc; e++) {\n"
        "        int en = (e + 1) % nc;\n"
        "        double x0 = poly_x[e], y0 = poly_y[e];\n"
        "        double x1 = poly_x[en], y1 = poly_y[en];\n"
        "        double cross = x0 * y1 - x1 * y0;\n"
        "        area += cross;\n"
        "        cx_sum += (x0 + x1) * cross;\n"
        "        cy_sum += (y0 + y1) * cross;\n"
        "        jxx12 += (x0*x0 + x0*x1 + x1*x1) * cross;\n"
        "        jyy12 += (y0*y0 + y0*y1 + y1*y1) * cross;\n"
        "        jxy24 += (2.0*x0*y0 + x0*y1 + x1*y0 + 2.0*x1*y1) * cross;\n"
        "    }\n"
        "    double moment_area = 0.5 * area;\n"
        "    double moment_cx = 0.0, moment_cy = 0.0;\n"
        "    double moment_qxx_mass = 0.0, moment_qxy_mass = 0.0;\n"
        "    double moment_qyy_mass = 0.0;\n"
        "    int moment_valid = 0;\n"
        "    if (raw_volume >= 2 && nc >= 3 && area > 0.0 && massi > 0.0 &&\n"
        "            isfinite(area) && isfinite(massi)) {\n"
        "        moment_cx = cx_sum / (3.0 * area);\n"
        "        moment_cy = cy_sum / (3.0 * area);\n"
        "        double jxx = jxx12 / 12.0;\n"
        "        double jyy = jyy12 / 12.0;\n"
        "        double jxy = jxy24 / 24.0;\n"
        "        double qxx_area = jxx - moment_area * moment_cx * moment_cx;\n"
        "        double qyy_area = jyy - moment_area * moment_cy * moment_cy;\n"
        "        double qxy_area = jxy - moment_area * moment_cx * moment_cy;\n"
        "        double mass_scale = massi / moment_area;\n"
        "        moment_qxx_mass = mass_scale * qxx_area;\n"
        "        moment_qxy_mass = mass_scale * qxy_area;\n"
        "        moment_qyy_mass = mass_scale * qyy_area;\n"
        "        double qdet = moment_qxx_mass * moment_qyy_mass -\n"
        "                      moment_qxy_mass * moment_qxy_mass;\n"
        "        moment_valid = isfinite(moment_cx) && isfinite(moment_cy) &&\n"
        "            isfinite(moment_qxx_mass) && isfinite(moment_qxy_mass) &&\n"
        "            isfinite(moment_qyy_mass) && moment_qxx_mass > 0.0 &&\n"
        "            moment_qyy_mass > 0.0 && qdet > 0.0;\n"
        "    }\n"
        "    double volume = 0.5 * fabs(area);\n",
        "signed polygon sums")
    modified_cu = replace_once(
        modified_cu,
        "    po->n_real_faces = n_real_faces;\n",
        "    po->n_real_faces = n_real_faces;\n"
        "    po->moment_area = moment_area;\n"
        "    po->moment_cx = moment_cx;\n"
        "    po->moment_cy = moment_cy;\n"
        "    po->moment_qxx_mass = moment_qxx_mass;\n"
        "    po->moment_qxy_mass = moment_qxy_mass;\n"
        "    po->moment_qyy_mass = moment_qyy_mass;\n"
        "    po->moment_valid = moment_valid;\n",
        "write moments to PolygonOut")

    helper = r'''static void gfs_audit_polygon_moments(GPUContext *ctx, int nbp)
{
    PolygonOut *poly = (PolygonOut *)malloc((size_t)nbp * sizeof(*poly));
    if (!poly) {
        fprintf(stderr, "[POLYGON_MOMENT_AUDIT] host allocation failed\n");
        exit(1);
    }
    CUDA_CHECK(cudaMemcpy(poly, ctx->d_polygons,
                          (size_t)nbp * sizeof(*poly), cudaMemcpyDeviceToHost));
    size_t invalid = 0;
    long double max_norm_error = 0.0L;
    for (int i = 0; i < nbp; ++i) {
        const PolygonOut *p = &poly[i];
        if (!p->moment_valid || p->n_corners < 3 ||
                p->n_corners > GPU_MAX_CORNERS || !(ctx->h_parts.mass[i] > 0.0)) {
            invalid++;
            continue;
        }
        GfsPolygonMomentsReference ref;
        if (!gfs_polygon_moments_reference(p->cx, p->cy,
                                           (size_t)p->n_corners, &ref)) {
            invalid++;
            continue;
        }
        long double r2max = 0.0L;
        for (int e = 0; e < p->n_corners; ++e) {
            const long double x = p->cx[e], y = p->cy[e];
            r2max = fmaxl(r2max, x*x + y*y);
        }
        const long double mass = ctx->h_parts.mass[i];
        const long double area_scale = ref.area;
        const long double length_scale = sqrtl(ref.area);
        const long double q_scale = mass * fmaxl(r2max, LDBL_MIN);
        const long double q_factor = mass / ref.area;
        const long double errors[] = {
            fabsl((long double)p->moment_area - ref.area) / area_scale,
            fabsl((long double)p->moment_cx - ref.centroid_x) / length_scale,
            fabsl((long double)p->moment_cy - ref.centroid_y) / length_scale,
            fabsl((long double)p->moment_qxx_mass - q_factor*ref.qxx_area) / q_scale,
            fabsl((long double)p->moment_qxy_mass - q_factor*ref.qxy_area) / q_scale,
            fabsl((long double)p->moment_qyy_mass - q_factor*ref.qyy_area) / q_scale
        };
        for (int k = 0; k < 6; ++k)
            max_norm_error = fmaxl(max_norm_error, errors[k]);
    }
    free(poly);
    const long double tolerance = 1.0e-8L;
    fprintf(stderr, "[POLYGON_MOMENT_AUDIT] cells=%d invalid=%zu max_norm_error=%.6Le tolerance=%.1Le %s\n",
            nbp, invalid, max_norm_error, tolerance,
            (invalid == 0 && max_norm_error <= tolerance) ? "PASS" : "FAIL");
    if (invalid != 0 || max_norm_error > tolerance) exit(1);
}

'''
    modified_cu = replace_once(
        modified_cu,
        "/* ================================================================\n"
        " *  Orchestrator: gpu_tessellate_and_build_faces\n",
        helper +
        "/* ================================================================\n"
        " *  Orchestrator: gpu_tessellate_and_build_faces\n",
        "moment audit helper insertion")
    modified_cu = replace_once(
        modified_cu,
        "    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[5], s));\n"
        "    CUDA_CHECK(cudaStreamSynchronize(s));\n"
        "    if (tess_profile) {\n",
        "    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[5], s));\n"
        "    CUDA_CHECK(cudaStreamSynchronize(s));\n"
        "    const char *moment_audit = getenv(\"GFS_GPU_MOMENT_AUDIT\");\n"
        "    if (!h_gfs_polygon_moment_audit_done && moment_audit && atoi(moment_audit) > 0) {\n"
        "        if (raw_geometry < 2) {\n"
        "            fprintf(stderr, \"[POLYGON_MOMENT_AUDIT] requires GFS_GPU_RAW_FACES=1\\n\");\n"
        "            exit(1);\n"
        "        }\n"
        "        gfs_audit_polygon_moments(ctx, nbp);\n"
        "        h_gfs_polygon_moment_audit_done = 1;\n"
        "    }\n"
        "    if (tess_profile) {\n",
        "one-shot GPU moment audit call")

    backup.mkdir(parents=False, exist_ok=False)
    shutil.copy2(header, backup / header.name)
    shutil.copy2(cuda, backup / cuda.name)
    temp_h = header.with_name(header.name + ".moment-tmp")
    temp_cu = cuda.with_name(cuda.name + ".moment-tmp")
    if temp_h.exists() or temp_cu.exists():
        raise RuntimeError("Refusing to overwrite an existing temporary source")
    temp_h.write_text(modified_h)
    temp_cu.write_text(modified_cu)
    os.replace(temp_h, header)
    os.replace(temp_cu, cuda)
    shutil.copy2(ref_source, ref_target)

    print(json.dumps({
        "worktree": str(root),
        "backup": str(backup),
        "source_hashes": {
            header.name: digest(header),
            cuda.name: digest(cuda),
            ref_target.name: digest(ref_target),
        },
        "moment_output": "raw polygon area, centroid offset, and mass central Q in PolygonOut",
        "audit": "one first-stage GPU polygon transfer checked against independent long-double reference",
    }, indent=2))


if __name__ == "__main__":
    main()

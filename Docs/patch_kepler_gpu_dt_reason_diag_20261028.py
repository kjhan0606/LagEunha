#!/usr/bin/env python3
"""Add opt-in local-CFL reason and face-neighbor telemetry to the q1 CUDA path.

Apply only to the authenticated 2026-10-26 native source reconstruction. The
diagnostic is disabled by default and prints only when the per-rank minimum dt
falls below 1e-5. Reason codes: 1 face CFL, 2 viscous CFL, 3 acceleration CFL,
4 entropy-change CFL, 5 W0 swept-volume geometry CFL.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


EXPECTED = {
    "Exam/exam_gpu.cu": "9d2937a7bf5e9a4365d020c86cf05c2e2494d7746a7e20f3eef321deeecd97e9",
    "Exam/exam_gpu_extract.c": "3f50f417275f2778758c920f513d4b61b87cbb651e4825b516febd85a44dc365",
    "Exam/exam_gpu.h": "30e85f684cdf940e4d2eba255c998590d22c126ef1612da4124f737181bb4ca9",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replace_once(root: Path, rel: str, old: str, new: str) -> None:
    path = root / rel
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{rel}: expected one anchor, found {count}: {old[:100]!r}")
    path.write_text(text.replace(old, new, 1))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: patch_kepler_gpu_dt_reason_diag_20261028.py SOURCE_ROOT")
    root = Path(sys.argv[1]).resolve()
    before = {}
    for rel, expected in EXPECTED.items():
        digest = sha((root / rel).read_bytes())
        if digest != expected:
            raise SystemExit(f"refusing unexpected source {rel}: {digest} != {expected}")
        before[rel] = digest

    replace_once(
        root,
        "Exam/exam_gpu.h",
        "    float  *dt_out;\n    double *vsig_max_out;",
        "    float  *dt_out;\n"
        "    int    *dt_reason_out;   /* opt-in CFL limiter code, one per owner cell */\n"
        "    int    *dt_face_j_out;   /* SoA neighbor index when reason=1, else -1 */\n"
        "    double *vsig_max_out;",
    )

    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    float *__restrict__ dt_out,\n"
        "    double *__restrict__ vsig_max_out,\n"
        "    double *__restrict__ vsmoothx_out, double *__restrict__ vsmoothy_out,\n"
        "    /* Optional per-face audit output, NULL in ordinary runs. */",
        "    float *__restrict__ dt_out,\n"
        "    int *__restrict__ dt_reason_out, int *__restrict__ dt_face_j_out,\n"
        "    double *__restrict__ vsig_max_out,\n"
        "    double *__restrict__ vsmoothx_out, double *__restrict__ vsmoothy_out,\n"
        "    /* Optional per-face audit output, NULL in ordinary runs. */",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    const double *__restrict__ pV64, double *__restrict__ out_Udot64)\n{\n"
        "    int i = blockIdx.x * blockDim.x + threadIdx.x;",
        "    double *__restrict__ pV64, double *__restrict__ out_Udot64, int dt_diag)\n{\n"
        "    int i = blockIdx.x * blockDim.x + threadIdx.x;",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    double dK_diss = 0;\n"
        "    float  my_dt = 1.0e10f;\n"
        "    double my_vsig_max = 0;",
        "    double dK_diss = 0;\n"
        "    float  my_dt = 1.0e10f;\n"
        "    int my_dt_reason = 0, my_dt_face_j = -1;\n"
        "    double my_vsig_max = 0;",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "        float dt_face = (float)(2.0 * Courant * dramp_cfl / vsig);\n"
        "        if (dt_face < my_dt) my_dt = dt_face;",
        "        float dt_face = (float)(2.0 * Courant * dramp_cfl / vsig);\n"
        "        if (dt_face < my_dt) {\n"
        "            my_dt = dt_face;\n"
        "            if (dt_diag) { my_dt_reason = 1; my_dt_face_j = j; }\n"
        "        }",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "        double rate=var_sweep_abs/ibp_volume;\n"
        "        if (var_geom_cfl>0) my_dt=fmin((double)my_dt,var_geom_cfl/rate);",
        "        double rate=var_sweep_abs/ibp_volume;\n"
        "        if (var_geom_cfl>0) {\n"
        "            double dt_geom=var_geom_cfl/rate;\n"
        "            if (dt_geom < (double)my_dt) {\n"
        "                my_dt=(float)dt_geom;\n"
        "                if (dt_diag) { my_dt_reason=5; my_dt_face_j=-1; }\n"
        "            }\n"
        "        }",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "            float dt_visc = (float)(0.5 * h_i * h_i / nu_eff);\n"
        "            if (dt_visc < my_dt) my_dt = dt_visc;",
        "            float dt_visc = (float)(0.5 * h_i * h_i / nu_eff);\n"
        "            if (dt_visc < my_dt) {\n"
        "                my_dt = dt_visc;\n"
        "                if (dt_diag) { my_dt_reason = 2; my_dt_face_j = -1; }\n"
        "            }",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "            if (dt_acc < (double)my_dt) my_dt = (float)dt_acc;",
        "            if (dt_acc < (double)my_dt) {\n"
        "                my_dt = (float)dt_acc;\n"
        "                if (dt_diag) { my_dt_reason = 3; my_dt_face_j = -1; }\n"
        "            }",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "            if (dt_dK < (double)my_dt) my_dt = (float)dt_dK;",
        "            if (dt_dK < (double)my_dt) {\n"
        "                my_dt = (float)dt_dK;\n"
        "                if (dt_diag) { my_dt_reason = 4; my_dt_face_j = -1; }\n"
        "            }",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    dK_out[i]       = phase1 ? 0.f : (float)dK_val;\n"
        "    dt_out[i]       = my_dt;\n"
        "    vsig_max_out[i] = my_vsig_max;",
        "    dK_out[i]       = phase1 ? 0.f : (float)dK_val;\n"
        "    dt_out[i]       = my_dt;\n"
        "    if (dt_diag) {\n"
        "        dt_reason_out[i] = my_dt_reason;\n"
        "        dt_face_j_out[i] = my_dt_face_j;\n"
        "    }\n"
        "    vsig_max_out[i] = my_vsig_max;",
    )

    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "#define ALLOC_PAIR_F(field) \\\n"
        "    CUDA_CHECK(cudaMallocHost(&h->field, sf)); \\\n"
        "    CUDA_CHECK(cudaMalloc(&d->field, sf))",
        "#define ALLOC_PAIR_F(field) \\\n"
        "    CUDA_CHECK(cudaMallocHost(&h->field, sf)); \\\n"
        "    CUDA_CHECK(cudaMalloc(&d->field, sf))\n"
        "#define ALLOC_PAIR_I(field) \\\n"
        "    CUDA_CHECK(cudaMallocHost(&h->field, n * sizeof(int))); \\\n"
        "    CUDA_CHECK(cudaMalloc(&d->field, n * sizeof(int)))",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    ALLOC_PAIR_F(die_out);  ALLOC_PAIR_F(dK_out);  ALLOC_PAIR_F(dt_out);",
        "    ALLOC_PAIR_F(die_out);  ALLOC_PAIR_F(dK_out);  ALLOC_PAIR_F(dt_out);\n"
        "    ALLOC_PAIR_I(dt_reason_out); ALLOC_PAIR_I(dt_face_j_out);",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "#undef ALLOC_PAIR_D\n#undef ALLOC_PAIR_F",
        "#undef ALLOC_PAIR_D\n#undef ALLOC_PAIR_F\n#undef ALLOC_PAIR_I",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    FREE_PAIR(die_out);  FREE_PAIR(dK_out);  FREE_PAIR(dt_out);",
        "    FREE_PAIR(die_out);  FREE_PAIR(dK_out);  FREE_PAIR(dt_out);\n"
        "    FREE_PAIR(dt_reason_out); FREE_PAIR(dt_face_j_out);",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dt_out, ctx->d_parts.dt_out, sf, cudaMemcpyDeviceToHost, s));\n"
        "    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.vsig_max_out, ctx->d_parts.vsig_max_out, sd, cudaMemcpyDeviceToHost, s));",
        "    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dt_out, ctx->d_parts.dt_out, sf, cudaMemcpyDeviceToHost, s));\n"
        "    {\n"
        "        const char *diag=getenv(\"EUNHA_GPU_DT_DIAG\");\n"
        "        if (diag && diag[0] && atoi(diag)>0) {\n"
        "            size_t si=(size_t)n*sizeof(int);\n"
        "            CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dt_reason_out, ctx->d_parts.dt_reason_out, si, cudaMemcpyDeviceToHost, s));\n"
        "            CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dt_face_j_out, ctx->d_parts.dt_face_j_out, si, cudaMemcpyDeviceToHost, s));\n"
        "        }\n"
        "    }\n"
        "    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.vsig_max_out, ctx->d_parts.vsig_max_out, sd, cudaMemcpyDeviceToHost, s));",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "    force_call++;\n"
        "    double var_geom_cfl=0;",
        "    force_call++;\n"
        "    const char *dt_diag_env=getenv(\"EUNHA_GPU_DT_DIAG\");\n"
        "    int dt_diag=dt_diag_env && dt_diag_env[0] && atoi(dt_diag_env)>0;\n"
        "    double var_geom_cfl=0;",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "        ctx->d_parts.die_out, ctx->d_parts.dK_out,\n"
        "        ctx->d_parts.dt_out,\n"
        "        ctx->d_parts.vsig_max_out,",
        "        ctx->d_parts.die_out, ctx->d_parts.dK_out,\n"
        "        ctx->d_parts.dt_out,\n"
        "        ctx->d_parts.dt_reason_out, ctx->d_parts.dt_face_j_out,\n"
        "        ctx->d_parts.vsig_max_out,",
    )
    replace_once(
        root,
        "Exam/exam_gpu.cu",
        "        thermal64,ctx->d_parts.E_inv_xx,ctx->d_parts.E_inv_xy,ctx->d_parts.E_inv_yx);",
        "        thermal64,ctx->d_parts.E_inv_xx,ctx->d_parts.E_inv_xy,ctx->d_parts.E_inv_yx,dt_diag);",
    )

    replace_once(
        root,
        "Exam/exam_gpu_extract.c",
        "    gpu_download_results(&g_gpu_ctx, nbp);\n\n"
        "    double t6_writeback = MPI_Wtime();",
        "    gpu_download_results(&g_gpu_ctx, nbp);\n"
        "    {\n"
        "        const char *diag=getenv(\"EUNHA_GPU_DT_DIAG\");\n"
        "        if (diag && diag[0] && atoi(diag)>0) {\n"
        "            int i, argmin=-1;\n"
        "            float dtmin=1.0e30f;\n"
        "            for (i=0;i<nbp;i++) {\n"
        "                float d=g_gpu_ctx.h_parts.dt_out[i];\n"
        "                if (isfinite(d) && d<dtmin) { dtmin=d; argmin=i; }\n"
        "            }\n"
        "            if (argmin>=0 && dtmin<1.0e-5f) {\n"
        "                int reason=g_gpu_ctx.h_parts.dt_reason_out[argmin];\n"
        "                int face_j=g_gpu_ctx.h_parts.dt_face_j_out[argmin];\n"
        "                long long id=g_gpu_ctx.h_parts.indx[argmin];\n"
        "                long long face_id=-1;\n"
        "                double fx=NAN, fy=NAN;\n"
        "                const char *label=reason==1 ? \"face\" : reason==2 ? \"visc\" :\n"
        "                    reason==3 ? \"accel\" : reason==4 ? \"dK\" :\n"
        "                    reason==5 ? \"geom_swept\" : \"unknown\";\n"
        "                if (face_j>=0 && face_j<g_gpu_ctx.max_particles) {\n"
        "                    face_id=g_gpu_ctx.h_parts.indx[face_j];\n"
        "                    fx=g_gpu_ctx.h_parts.x[face_j];\n"
        "                    fy=g_gpu_ctx.h_parts.y[face_j];\n"
        "                }\n"
        "                printf(\"[GPU_DT_LIMIT] rank=%d call=%d stage=%d dt=%.9g reason=%s(code=%d) cell_soa=%d id=%lld xy=%.9g,%.9g V=%.9g rho=%.9g P=%.9g cs=%.9g face_soa=%d face_id=%lld face_xy=%.9g,%.9g\\n\",\n"
        "                    MYID(simpar),gpu_call_count+1,gpu_call_count%4,dtmin,label,reason,argmin,id,\n"
        "                    g_gpu_ctx.h_parts.x[argmin],g_gpu_ctx.h_parts.y[argmin],\n"
        "                    g_gpu_ctx.h_parts.volume[argmin],g_gpu_ctx.h_parts.den[argmin],\n"
        "                    g_gpu_ctx.h_parts.pressure[argmin],g_gpu_ctx.h_parts.csound[argmin],\n"
        "                    face_j,face_id,fx,fy);\n"
        "                fflush(stdout);\n"
        "            }\n"
        "        }\n"
        "    }\n\n"
        "    double t6_writeback = MPI_Wtime();",
    )

    after = {rel: sha((root / rel).read_bytes()) for rel in EXPECTED}
    receipt = {
        "source_root": str(root),
        "before_sha256": before,
        "after_sha256": after,
        "diagnostic_env": "EUNHA_GPU_DT_DIAG=1",
        "reason_codes": {
            "1": "face CFL",
            "2": "viscous CFL",
            "3": "acceleration CFL",
            "4": "entropy-change dK CFL",
            "5": "W0 swept-volume geometry CFL",
        },
    }
    out = root / "dt_reason_diag_patch_receipt.json"
    out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()

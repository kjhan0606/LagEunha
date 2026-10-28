"""Add diagnostic-only CUDA tessellation work counters to the guarded v3 source."""

import argparse
import hashlib
import json
from pathlib import Path


PARENT_SHA = "188c0c064deec4c21bcdf51c1b67377e58ed38aaf462fe9f241602f90b24f926"


def once(s, old, new, label):
    n = s.count(old)
    if n != 1:
        raise ValueError(f"{label}: expected once, found {n}")
    return s.replace(old, new, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_dir", type=Path)
    ap.add_argument("--seed3", action="store_true")
    args = ap.parse_args()
    path = args.source_dir / "Exam" / "exam_gpu.cu"
    original = path.read_bytes()
    if hashlib.sha256(original).hexdigest() != PARENT_SHA:
        raise ValueError("Expected guarded v3 CUDA source")
    s = original.decode()
    s = once(s,
        "static __device__ unsigned long long d_gfs_stencil_expanded_cells;\n",
        "static __device__ unsigned long long d_gfs_stencil_expanded_cells;\n"
        "/* Diagnostic-only totals: pre visits/clips, seed visits/clips, adaptive\n"
        " * visits/clips, wall clips, input edges processed, and box retries. */\n"
        "static __device__ unsigned long long d_gfs_tess_work[9];\n",
        "device counters")
    s = once(s,
        "int    nc = 4, stencil_pad = GPU_SEED_PAD, seed_candidates = 0, box_retry = 0;\n",
        "int    nc = 4, stencil_pad = GPU_SEED_PAD, seed_candidates = 0, box_retry = 0;\n"
        "unsigned int diag_pre_visits=0,diag_pre_clips=0,diag_seed_visits=0,diag_seed_clips=0;\n"
        "unsigned int diag_adapt_visits=0,diag_adapt_clips=0,diag_wall_clips=0,diag_edges=0;\n",
        "per-site counters")
    s = once(s,
        "int j=cell_particle_index[k];if(pindx[j]==indx_i)continue;\n",
        "int j=cell_particle_index[k];if(pindx[j]==indx_i)continue;\n"
        "                    diag_pre_visits++;\n",
        "preclip visits")
    s = once(s,
        "int j=nearest[s];if(j<0 || offset2[s]>maxdist2)continue;\n",
        "int j=nearest[s];if(j<0 || offset2[s]>maxdist2)continue;\n"
        "            diag_pre_clips++;diag_edges+=(unsigned int)nc;\n",
        "preclip edges")
    wall_anchor = "nc=gpu_clip_voronoi_halfplane(nx,ny,threshold,-1000-k,\n"
    if s.count(wall_anchor) != 2:
        raise ValueError("Expected early and final fixed-wall clip sites")
    s = s.replace(wall_anchor,
        "diag_wall_clips++;diag_edges+=(unsigned int)nc;\n"
        "                " + wall_anchor)
    s = once(s,
        "seed_candidates++;\n",
        "seed_candidates++;diag_seed_visits++;\n",
        "seed visits")
    s = once(s,
        "if (drad2 > maxdist2) continue;\n\n                nc = gpu_clip_voronoi_halfplane",
        "if (drad2 > maxdist2) continue;\n"
        "                diag_seed_clips++;diag_edges+=(unsigned int)nc;\n\n"
        "                nc = gpu_clip_voronoi_halfplane",
        "seed clips")
    s = once(s,
        "if (pindx[j] == indx_i) continue;\n                        const double dx = px[j] - xi;",
        "if (pindx[j] == indx_i) continue;\n"
        "                        diag_adapt_visits++;\n"
        "                        const double dx = px[j] - xi;",
        "adaptive visits")
    s = once(s,
        "if (wfrac * wfrac * dist2 > maxdist2) continue;\n\n                        nc = gpu_clip_voronoi_halfplane",
        "if (wfrac * wfrac * dist2 > maxdist2) continue;\n"
        "                        diag_adapt_clips++;diag_edges+=(unsigned int)nc;\n\n"
        "                        nc = gpu_clip_voronoi_halfplane",
        "adaptive clips")
    s = once(s,
        "    } /* adaptive initial box */\n    if (seed_candidates > GPU_MAX_NEIGH)",
        "    } /* adaptive initial box */\n"
        "    atomicAdd(&d_gfs_tess_work[0],(unsigned long long)diag_pre_visits);\n"
        "    atomicAdd(&d_gfs_tess_work[1],(unsigned long long)diag_pre_clips);\n"
        "    atomicAdd(&d_gfs_tess_work[2],(unsigned long long)diag_seed_visits);\n"
        "    atomicAdd(&d_gfs_tess_work[3],(unsigned long long)diag_seed_clips);\n"
        "    atomicAdd(&d_gfs_tess_work[4],(unsigned long long)diag_adapt_visits);\n"
        "    atomicAdd(&d_gfs_tess_work[5],(unsigned long long)diag_adapt_clips);\n"
        "    atomicAdd(&d_gfs_tess_work[6],(unsigned long long)diag_wall_clips);\n"
        "    atomicAdd(&d_gfs_tess_work[7],(unsigned long long)diag_edges);\n"
        "    atomicAdd(&d_gfs_tess_work[8],(unsigned long long)box_retry);\n"
        "    if (seed_candidates > GPU_MAX_NEIGH)",
        "aggregate work")
    s = once(s,
        "    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[0], s));\n",
        "    unsigned long long zero_work[9]={};\n"
        "    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_tess_work,zero_work,sizeof(zero_work),0,\n"
        "                  cudaMemcpyHostToDevice,s));\n"
        "    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[0], s));\n",
        "reset work")
    s = once(s,
        "        unsigned long long h_candidate_fallbacks = 0;\n",
        "        unsigned long long h_work[9]={};\n"
        "        CUDA_CHECK(cudaMemcpyFromSymbol(h_work,d_gfs_tess_work,sizeof(h_work)));\n"
        "        fprintf(stderr,\"[GPU_TESS_WORK] pre_visit=%llu pre_clip=%llu seed_visit=%llu seed_clip=%llu adapt_visit=%llu adapt_clip=%llu wall_clip=%llu clip_edges=%llu box_retry=%llu\\n\",\n"
        "                h_work[0],h_work[1],h_work[2],h_work[3],h_work[4],\n"
        "                h_work[5],h_work[6],h_work[7],h_work[8]);\n"
        "        unsigned long long h_candidate_fallbacks = 0;\n",
        "print work")
    if args.seed3:
        s = once(s,
            "#define GPU_MAX_NEIGH 512\n#ifdef GFS_FAST_SEED3\n",
            "#define GPU_MAX_NEIGH 512\n#define GFS_FAST_SEED3 1\n#ifdef GFS_FAST_SEED3\n",
            "certified seed3 selector")
    path.write_text(s)
    print(json.dumps({"parent_cuda_sha256": PARENT_SHA,
                      "patched_cuda_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                      "seed3": args.seed3,
                      "scope": "Diagnostic counters only; optional existing 3x3 seed selector."}, indent=2))


if __name__ == "__main__":
    main()

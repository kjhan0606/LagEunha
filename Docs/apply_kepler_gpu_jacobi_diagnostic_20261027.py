"""Attach the read-only GPU Jacobi probe to an isolated native worktree."""

import argparse
import hashlib
from pathlib import Path
import shutil


DOCS = Path(__file__).resolve().parent


def replace_once(path, old, new):
    text = path.read_text()
    if text.count(old) != 1:
        raise ValueError(f"Expected one unmodified anchor in {path}: {old[:70]!r}")
    path.write_text(text.replace(old, new))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("worktree", type=Path)
    args = ap.parse_args()
    root = args.worktree.resolve()
    if not (root / ".kepler_bounded_overlay.json").exists():
        raise ValueError("Accepted bounded source overlay is missing")
    source = DOCS / "kepler_gpu_jacobi_probe_20261027.cu"
    dest = root / "Exam/gfs_volume_rate_probe.cu"
    if dest.exists() and dest.read_bytes() != source.read_bytes():
        raise ValueError("Existing CUDA probe differs from archived candidate")
    shutil.copy2(source, dest)

    header = root / "Exam/exam_gpu.h"
    anchor = "void gpu_download_results(GPUContext *ctx, int n_particles);\n"
    replace_once(header, anchor, anchor +
        "/* Read-only volume-rate Jacobi probe on the current device FaceCSR. */\n"
        "int gpu_probe_volume_rate_jacobi(GPUContext *ctx,\n"
        "    const unsigned char *host_active, int n_particles, int n_total,\n"
        "    int iterations, double omega, double *host_rate,\n"
        "    unsigned long long *host_unsupported);\n")

    makefile = root / "Exam/Makefile"
    replace_once(makefile,
        "EXOBJ      += exam_gpu.o exam_gpu_extract.o\n",
        "EXOBJ      += exam_gpu.o exam_gpu_extract.o gfs_volume_rate_probe.o\n")
    anchor = "exam_gpu_extract.o: exam_gpu_extract.c exam_gpu.h gfs_thermal64.h gfs_pair.h gfs_moment.h gfs_moment_density.h gfs_moment_initial.h gfs_dg_euler.h gfs_enriched_euler.h gfs_enriched_coupling.h\n\t$(CC) $(CCFLAGS) $(INCLUDES) -c $< -o $@\n"
    replace_once(makefile, anchor, anchor +
        "\ngfs_volume_rate_probe.o: gfs_volume_rate_probe.cu exam_gpu.h\n"
        "\t$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@\n")

    extract = root / "Exam/exam_gpu_extract.c"
    anchor = "static int gpu_call_count = 0;\n"
    helper = r'''/* Optional diagnostic: one native FaceCSR weight-rate solve, no state edit. */
static void gpu_volume_rate_probe_once(SimParameters *simpar, int nbp, int n_total)
{
    static int done = 0;
    const char *path = getenv("GFS_GPU_WDOT_PROBE_MASK");
    if (!path || !path[0] || done || MYID(simpar) != 0 || NID(simpar) != 1)
        return;
    done = 1;
    FILE *input = fopen(path, "rb");
    if (!input || fseek(input, 0, SEEK_END) != 0) MPI_Abort(MPI_COMM(simpar), 97);
    long count = ftell(input);
    if (count <= 0 || count > 100000000 || fseek(input, 0, SEEK_SET) != 0)
        MPI_Abort(MPI_COMM(simpar), 97);
    unsigned char *by_id = (unsigned char *)malloc((size_t)count);
    unsigned char *active = (unsigned char *)calloc((size_t)n_total, 1);
    double *rate = (double *)malloc((size_t)nbp*sizeof(double));
    if (!by_id || !active || !rate ||
        fread(by_id, 1, (size_t)count, input) != (size_t)count)
        MPI_Abort(MPI_COMM(simpar), 97);
    fclose(input);
    size_t stride = TVORORK4_DDINFO(simpar)[0].n_size;
    char *base = (char *)VORORK4_TBP(simpar);
    int nactive = 0;
    for (int i=0; i<nbp; ++i) {
        treevorork4particletype *bp =
            (treevorork4particletype *)(base+(size_t)i*stride);
        unsigned long long id = ((unsigned long long)PINDX(bp)) & ((1ULL<<56)-1);
        if (id >= (unsigned long long)count) MPI_Abort(MPI_COMM(simpar), 97);
        active[i] = by_id[id] ? 1 : 0;
        nactive += active[i];
    }
    unsigned long long unsupported = 0;
    double started = MPI_Wtime();
    int status = gpu_probe_volume_rate_jacobi(&g_gpu_ctx, active,
        nbp, n_total, 16, .8, rate, &unsupported);
    double elapsed = MPI_Wtime()-started;
    if (status || unsupported) {
        fprintf(stderr,"[GPU_WDOT_PROBE] status=%d unsupported=%llu\n",
                status, unsupported);
        MPI_Abort(MPI_COMM(simpar), 97);
    }
    FILE *output = fopen("gpu_wdot_probe_r0_call1.tsv", "wx");
    if (!output) MPI_Abort(MPI_COMM(simpar), 97);
    fprintf(output,"id\twdot\n");
    for (int i=0; i<nbp; ++i) {
        treevorork4particletype *bp =
            (treevorork4particletype *)(base+(size_t)i*stride);
        unsigned long long id = ((unsigned long long)PINDX(bp)) & ((1ULL<<56)-1);
        fprintf(output,"%llu\t%.17g\n",id,rate[i]);
    }
    fclose(output);
    fprintf(stderr,"[GPU_WDOT_PROBE] active=%d faces=%d jacobi=16 omega=0.8 elapsed=%.6f s\n",
            nactive,g_gpu_ctx.d_faces.n_faces_total,elapsed);
    free(rate); free(active); free(by_id);
}

'''
    replace_once(extract, anchor, helper + anchor)
    anchor = "    double t4_kernel = MPI_Wtime();\n    /* --- Launch force kernel --- */\n"
    replace_once(extract, anchor,
        "    gpu_volume_rate_probe_once(simpar, nbp, n_total);\n"+anchor)
    print("Added opt-in read-only GPU Jacobi probe to", root)


if __name__ == "__main__":
    main()

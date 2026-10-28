"""Add a bounded runtime Jacobi-iteration choice to the isolated v2 source.

Run only after patch_kepler_energy_neutral_varpress_guard_v2_20261027.py.
The default remains 16 so W0 and existing weighted runs retain their mode.
"""

import argparse
import hashlib
import json
from pathlib import Path


EXPECTED_V2_CUDA_SHA = "188c0c064deec4c21bcdf51c1b67377e58ed38aaf462fe9f241602f90b24f926"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def replace_once(source, old, new, name):
    count = source.count(old)
    if count != 1:
        raise ValueError(f"{name}: expected one source anchor, found {count}")
    return source.replace(old, new)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path)
    args = parser.parse_args()
    root = args.source_root.resolve()
    cuda = root / "Exam/exam_gpu.cu"
    if digest(cuda.read_bytes()) != EXPECTED_V2_CUDA_SHA:
        raise ValueError("Expected verified v2 finite-weight pressure guard")

    extract = root / "Exam/exam_gpu_extract.c"
    original_extract = extract.read_bytes()
    source = original_extract.decode()
    anchor = "    int status=gpu_probe_energy_neutral_rate(&g_gpu_ctx,active,nbp,n_total,\n        16,.8,rate,gradient,&alpha,&unsupported);"
    replacement = """    static int iterations=0;
    if(!iterations){
        const char *setting=getenv("GFS_GPU_VOLW2_ITERS");
        if(!setting || !setting[0] || strcmp(setting,"16")==0) iterations=16;
        else if(strcmp(setting,"32")==0) iterations=32;
        else if(strcmp(setting,"64")==0) iterations=64;
        else {
            fprintf(stderr,"[GPU_VOLW2] GFS_GPU_VOLW2_ITERS must be 16, 32, or 64\\n");
            MPI_Abort(MPI_COMM(simpar),97);
        }
        fprintf(stderr,"[GPU_VOLW2] Jacobi iterations=%d\\n",iterations);
    }
    int status=gpu_probe_energy_neutral_rate(&g_gpu_ctx,active,nbp,n_total,
        iterations,.8,rate,gradient,&alpha,&unsupported);"""
    source = replace_once(source, anchor, replacement, "rate call")
    extract.write_bytes(source.encode())

    exam = root / "Exam/exam.c"
    original_exam = exam.read_bytes()
    source = original_exam.decode()
    source = replace_once(source,
        "[GPU_VOLW2] experimental pressure-work-neutral Jacobi16 RK4 weights enabled",
        "[GPU_VOLW2] experimental pressure-work-neutral configurable-Jacobi RK4 weights enabled",
        "mode banner")
    exam.write_bytes(source.encode())
    print(json.dumps({
        "v2_cuda_sha256": EXPECTED_V2_CUDA_SHA,
        "extract_parent_sha256": digest(original_extract),
        "extract_patched_sha256": digest(extract.read_bytes()),
        "exam_parent_sha256": digest(original_exam),
        "exam_patched_sha256": digest(exam.read_bytes()),
        "allowed_iterations": [16, 32, 64],
        "default_iterations": 16,
    }, indent=2))


if __name__ == "__main__":
    main()

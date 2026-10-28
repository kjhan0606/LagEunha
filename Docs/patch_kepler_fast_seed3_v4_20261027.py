"""Enable certified 3x3 seed search after the guarded v3 weighted patch."""

import argparse
import hashlib
import json
from pathlib import Path


V2_GUARDED_CUDA_SHA = "188c0c064deec4c21bcdf51c1b67377e58ed38aaf462fe9f241602f90b24f926"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_root", type=Path)
    args = ap.parse_args()
    file = args.source_root.resolve() / "Exam/exam_gpu.cu"
    original = file.read_bytes()
    if digest(original) != V2_GUARDED_CUDA_SHA:
        raise ValueError("Expected verified v2/v3 CUDA source before seed3 patch")
    source = original.decode()
    old = "#define GPU_MAX_NEIGH 512\n#ifdef GFS_FAST_SEED3\n"
    new = "#define GPU_MAX_NEIGH 512\n#define GFS_FAST_SEED3 1\n#ifdef GFS_FAST_SEED3\n"
    if source.count(old) != 1:
        raise ValueError("Expected one certified seed3 source anchor")
    patched = source.replace(old, new).encode()
    file.write_bytes(patched)
    print(json.dumps({"parent_cuda_sha256": digest(original),
                      "patched_cuda_sha256": digest(patched),
                      "change": "Enable existing 3x3 seed with conservative adaptive completion; no clipping or force law changed."},
                     indent=2))


if __name__ == "__main__":
    main()

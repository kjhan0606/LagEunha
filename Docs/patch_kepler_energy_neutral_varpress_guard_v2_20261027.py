"""Allow the derived finite-weight variational branch only in GPU_VOLW2 mode."""

import argparse
import hashlib
import json
from pathlib import Path


PARENT_SHA = "5ebc8d575bbb8bb361069923384b450682f950b36625c82c7396e50cb1c33f67"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"Expected one source anchor, got {source.count(old)}")
    return source.replace(old, new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_root", type=Path)
    args = ap.parse_args()
    root = args.source_root.resolve()
    file = root / "Exam/exam_gpu.cu"
    original = file.read_bytes()
    if sha(original) != PARENT_SHA:
        raise ValueError("Parent v1 CUDA source hash changed")
    source = original.decode()
    source = replace_once(
        source,
        'extern "C" int gpu_varpress_w0_enabled(void) {\n',
        'extern "C" int gfs_volw2_is_active(void);\n'
        'extern "C" int gpu_varpress_w0_enabled(void) {\n',
    )
    old = '''    if (!params->phase1 || params->av_mode != 5 || !params->use_muscl ||
        params->alphavis != 0.0 || params->betavis != 0.0 ||
        params->nu_phys != 0.0 || params->cd_amax != 0.0 ||
        params->hyperv_alpha != 0.0 || params->face_charge ||
        params->w2lin || params->laguerre_rot) {
        fprintf(stderr, "[VARPRESS_W0] requires phase1/MUSCL/av5, W0, Re=inf, no AV/weight rates\\n");
        exit(1);
    }
'''
    new = '''    const bool weighted_varpress = gfs_volw2_is_active() &&
        params->w2lin && params->laguerre_rot && params->geom_fv;
    if (!params->phase1 || params->av_mode != 5 || !params->use_muscl ||
        params->alphavis != 0.0 || params->betavis != 0.0 ||
        params->nu_phys != 0.0 || params->cd_amax != 0.0 ||
        params->hyperv_alpha != 0.0 || params->face_charge ||
        (!weighted_varpress && (params->w2lin || params->laguerre_rot))) {
        fprintf(stderr, "[VARPRESS_W0] requires phase1/MUSCL/av5, Re=inf, no AV; finite weights require projected GPU_VOLW2 geometry\\n");
        exit(1);
    }
'''
    source = replace_once(source, old, new)
    patched = source.encode()
    file.write_bytes(patched)
    receipt = {"parent_source_sha256": PARENT_SHA,
               "patched_source_sha256": sha(patched),
               "change": "Finite-weight variational force allowed only with active GPU_VOLW2 and weight-rate face geometry; all baseline W0 guards remain."}
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()

"""Repair the W0-only launch guard in the already packaged native pilot.

The mode is accepted only if the GPU stage was explicitly activated by
GFS_GPU_VOLW2 and supplies the matching weight-rate face kinematics. This
does not change the W0 reference path when the mode is disabled.
"""

import argparse
import hashlib
import json
from pathlib import Path


OLD_SHA = "5ebc8d575bbb8bb361069923384b450682f950b36625c82c7396e50cb1c33f67"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"Expected exactly one guard anchor: {old[:80]!r}")
    return source.replace(old, new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_tree", type=Path)
    ap.add_argument("--receipt", type=Path, required=True)
    args = ap.parse_args()
    target = args.source_tree / "Exam/exam_gpu.cu"
    original = target.read_bytes()
    if sha(original) != OLD_SHA:
        raise ValueError("Wrong packaged energy-neutral source")
    source = original.decode()
    source = once(source,
        'extern "C" int gpu_varpress_w0_enabled(void) {',
        'extern "C" int gfs_volw2_is_active(void);\n'
        'extern "C" int gpu_varpress_w0_enabled(void) {')
    source = once(source,
        '        params->w2lin || params->laguerre_rot) {',
        '        ((params->w2lin || params->laguerre_rot) &&\n'
        '         !(gfs_volw2_is_active() && params->w2lin &&\n'
        '           params->laguerre_rot && params->geom_fv))) {')
    source = once(source,
        '        fprintf(stderr, "[VARPRESS_W0] experimental smooth fixed-mass variational pressure active\\n");',
        '        fprintf(stderr, gfs_volw2_is_active() ?\n'
        '            "[VARPRESS_W0] gated weighted fixed-mass variational pressure active\\n" :\n'
        '            "[VARPRESS_W0] experimental smooth fixed-mass variational pressure active\\n");')
    target.write_text(source)
    receipt = {"target": str(target), "original_sha256": OLD_SHA,
               "patched_sha256": sha(target.read_bytes()),
               "scope": "W0-only launch guard now permits the explicitly activated weighted pilot."}
    with args.receipt.open("x") as output:
        json.dump(receipt, output, indent=2)
        output.write("\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()

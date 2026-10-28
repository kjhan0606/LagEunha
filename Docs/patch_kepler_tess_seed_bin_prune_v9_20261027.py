"""Skip certified non-cutting W=0 seed bins in the guarded native CUDA source."""

import argparse
import hashlib
import json
from pathlib import Path


PARENT_SHA = '188c0c064deec4c21bcdf51c1b67377e58ed38aaf462fe9f241602f90b24f926'
ANCHOR = '''            int cx = ix + dx_c;
            if (cx < 0 || cx >= mx) continue;
            int cell_id = cx + mx * cy;
'''
REPLACEMENT = '''            int cx = ix + dx_c;
            if (cx < 0 || cx >= mx) continue;
            /* For W=0 the owner/competitor bisector lies d/2 away. A bin
             * whose nearest possible site is beyond 2R cannot intersect
             * the current polygon (R is its largest vertex radius). Skipped
             * candidates would all fail the later individual drad2 test.
             * Never use a finite bin rectangle for clamped edge bins. */
            if (w2i == 0.0 && global_w2_max == 0.0 &&
                cx > 0 && cx < mx-1 && cy > 0 && cy < my-1) {
                double bx0 = xmin + (double)cx * cellsize;
                double by0 = ymin + (double)cy * cellsize;
                double dx_min = fmax(0.0, fmax(bx0-xi, xi-(bx0+cellsize)));
                double dy_min = fmax(0.0, fmax(by0-yi, yi-(by0+cellsize)));
                double dmin2 = dx_min*dx_min + dy_min*dy_min;
                if (dmin2 > 4.0*maxdist2*(1.0+64.0*DBL_EPSILON)) continue;
            }
            int cell_id = cx + mx * cy;
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source_dir', type=Path)
    args = ap.parse_args()
    path = args.source_dir / 'Exam' / 'exam_gpu.cu'
    original = path.read_bytes()
    if hashlib.sha256(original).hexdigest() != PARENT_SHA:
        raise ValueError('Expected guarded v3 CUDA source')
    text = original.decode()
    if text.count(ANCHOR) != 1:
        raise ValueError(f'Expected one seed-loop anchor, found {text.count(ANCHOR)}')
    text = text.replace(ANCHOR, REPLACEMENT, 1)
    path.write_text(text)
    print(json.dumps({'parent_cuda_sha256': PARENT_SHA,
                      'patched_cuda_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                      'scope': 'W=0-certified seed-bin skipping only; same 5x5 preclip and exact adaptive search.'},
                     indent=2))


if __name__ == '__main__':
    main()

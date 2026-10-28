"""Restore the validated 64-corner polygon capacity in isolated weighted source."""

import hashlib
import json
from pathlib import Path
import sys


EXPECTED='cf86b1c3e7f1704c655b7cb9842bde1e0df4c028e3fe6fb475824ca04fa624b5'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if len(sys.argv)!=2:raise SystemExit('Expected isolated source root')
    path=Path(sys.argv[1]).resolve()/'Exam/exam_gpu.h'
    if digest(path)!=EXPECTED:raise ValueError('Unexpected weighted header hash')
    text=path.read_text()
    old='''/* Maximum polygon corners per particle (box starts with 4,
   each neighbor adds at most 1 → 48 handles ~44 neighbors) */
#define GPU_MAX_CORNERS 48'''
    new='''/* Capacity for intermediate clipped polygons. This bounded weighted
 * branch uses the previously validated near32 clipping order and 64 slots. */
#define GPU_MAX_CORNERS 64'''
    if text.count(old)!=1:raise ValueError('Corner-capacity anchor changed')
    path.write_text(text.replace(old,new))
    print(json.dumps({'parent_header_sha256':EXPECTED,'header_sha256':digest(path),
                      'old_capacity':48,'new_capacity':64,
                      'scope':'Preserve existing weighted interface; rebuild all '
                              'objects that include the changed polygon layout.'},indent=2))


if __name__=='__main__':main()

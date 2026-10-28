"""Specialize the exact raw kernel for the established limiter=1 mode."""

import argparse
import hashlib
import json
from pathlib import Path


PARENT_SHA = '11d749f6451aa85b2462dc3645712c3e93cc224f5c5eff66889542032efcf8b5'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source_dir', type=Path)
    args = ap.parse_args()
    path = args.source_dir / 'Exam' / 'exam_gpu.cu'
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != PARENT_SHA:
        raise ValueError('Unexpected v5 CUDA source')
    source = data.decode()
    changes = (
        ('template<int EXACT_RAW>\n__global__ __launch_bounds__(256, 2)\nvoid voronoi_tessellate_kernel(',
         'template<int EXACT_RAW, int LIMIT_CONST>\n__global__ __launch_bounds__(256, 2)\nvoid voronoi_tessellate_kernel(', 1),
        ('int thermal64, double *__restrict__ out_V64)\n{\n    int i = blockIdx.x',
         'int thermal64, double *__restrict__ out_V64)\n{\n'
         '    const int limit_velocity = LIMIT_CONST >= 0 ? LIMIT_CONST : limit_velocity_runtime;\n'
         '    int i = blockIdx.x', 1),
        ('int gradient_method, int raw_volume, int limit_velocity, int merge_small_faces,',
         'int gradient_method, int raw_volume, int limit_velocity_runtime, int merge_small_faces,', 1),
        ('voronoi_tessellate_kernel<1><<<', 'voronoi_tessellate_kernel<1,1><<<', 2),
        ('voronoi_tessellate_kernel<0><<<', 'voronoi_tessellate_kernel<0,-1><<<', 2),
    )
    for old, new, count in changes:
        if source.count(old) != count:
            raise ValueError(f'Expected {count} occurrences of {old[:60]!r}, found {source.count(old)}')
        source = source.replace(old, new, count)
    old = 'if (raw_geometry >= 2) {\n        voronoi_tessellate_kernel<1,1><<<'
    if source.count(old) != 2:
        raise ValueError('Expected both raw specialization launches')
    source = source.replace(old,
        'if (raw_geometry >= 2 && limit_velocity == 1) {\n'
        '        voronoi_tessellate_kernel<1,1><<<', 1)
    path.write_text(source)
    print(json.dumps({
        'parent_cuda_sha256': PARENT_SHA,
        'patched_cuda_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'mode': 'exact raw convex clip and limiter=1 use a dedicated template; other modes retain legacy clipper',
    }, indent=2))


if __name__ == '__main__':
    main()

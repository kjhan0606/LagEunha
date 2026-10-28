"""Specialize the native limiter=1 Kepler path without changing clipping."""

import argparse
import hashlib
import json
from pathlib import Path


PARENT_SHA = '188c0c064deec4c21bcdf51c1b67377e58ed38aaf462fe9f241602f90b24f926'
LAUNCH = 'voronoi_tessellate_kernel<<<gridSize, blockSize, 0, s>>>(\n'


def replace_once(source, old, new, name):
    count = source.count(old)
    if count != 1:
        raise ValueError(f'{name}: expected one occurrence, found {count}')
    return source.replace(old, new, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source_dir', type=Path)
    args = ap.parse_args()
    path = args.source_dir / 'Exam' / 'exam_gpu.cu'
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != PARENT_SHA:
        raise ValueError('Unexpected v3 CUDA source')
    source = raw.decode()
    source = replace_once(source,
        '__global__ __launch_bounds__(256, 2)\nvoid voronoi_tessellate_kernel(',
        'template<int LIMIT_CONST>\n__global__ __launch_bounds__(256, 2)\nvoid voronoi_tessellate_kernel(',
        'kernel template')
    source = replace_once(source,
        'int gradient_method, int raw_volume, int limit_velocity, int merge_small_faces,',
        'int gradient_method, int raw_volume, int limit_velocity_runtime, int merge_small_faces,',
        'limiter argument')
    source = replace_once(source,
        'int thermal64, double *__restrict__ out_V64)\n{\n    int i = blockIdx.x',
        'int thermal64, double *__restrict__ out_V64)\n{\n'
        '    const int limit_velocity = LIMIT_CONST >= 0 ? LIMIT_CONST : limit_velocity_runtime;\n'
        '    int i = blockIdx.x',
        'compile-time limiter value')
    if source.count(LAUNCH) != 2:
        raise ValueError('Expected two native tessellation launches')
    # Stage/force tessellation: specialize the exact raw mode used in Kepler.
    source = source.replace(LAUNCH,
        'if (raw_geometry >= 2 && limit_velocity == 1) {\n'
        '        voronoi_tessellate_kernel<1><<<gridSize, blockSize, 0, s>>>(\n', 1)
    first = source.index('        voronoi_tessellate_kernel<1><<<')
    end = source.index('ctx->d_parts.E_inv_xy);', first) + len('ctx->d_parts.E_inv_xy);')
    block = source[first:end]
    source = source[:end] + '\n    } else {\n' + block.replace('kernel<1>', 'kernel<-1>') + '\n    }' + source[end:]
    # End-of-step raw volume uses limiter=1; retain generic path otherwise.
    source = source.replace(LAUNCH,
        'if (raw_geometry >= 2) {\n'
        '        voronoi_tessellate_kernel<1><<<gridSize, blockSize, 0, s>>>(\n', 1)
    second = source.rindex('        voronoi_tessellate_kernel<1><<<')
    end2 = source.index('ctx->d_parts.E_inv_xy);', second) + len('ctx->d_parts.E_inv_xy);')
    block2 = source[second:end2]
    source = source[:end2] + '\n    } else {\n' + block2.replace('kernel<1>', 'kernel<-1>') + '\n    }' + source[end2:]
    if source.count('voronoi_tessellate_kernel<1><<<') != 2 or source.count('voronoi_tessellate_kernel<-1><<<') != 2:
        raise ValueError('Specialized and generic launches not installed')
    path.write_text(source)
    print(json.dumps({
        'parent_cuda_sha256': PARENT_SHA,
        'patched_cuda_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'scope': 'Only compile-time limiter=1 specialization of the exact raw path; the original clipper and all geometric guards remain unchanged.',
    }, indent=2))


if __name__ == '__main__':
    main()

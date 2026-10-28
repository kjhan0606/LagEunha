"""Patch an isolated v5 CUDA tree with an opt-in validated bin certificate."""

import hashlib
import json
from pathlib import Path
import shutil
import sys

EXPECTED_CUDA = '188c0c064deec4c21bcdf51c1b67377e58ed38aaf462fe9f241602f90b24f926'
HEADER = 'gfs_seed_bin_certificate_20261027.h'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f'Expected exactly one CUDA anchor, got {source.count(old)}')
    return source.replace(old, new)


def main():
    if len(sys.argv) != 2:
        raise SystemExit('Expected isolated build tree')
    root = Path(sys.argv[1]).resolve()
    cuda = root/'Exam/exam_gpu.cu'
    header = root/'Exam'/HEADER
    if digest(cuda) != EXPECTED_CUDA or header.exists():
        raise ValueError('Unexpected parent CUDA source or preexisting header')
    source = cuda.read_text()
    source = once(source, '#include "gfs_thermal64.h"',
                  '#include "gfs_thermal64.h"\n#include "'+HEADER+'"')
    anchor = 'static __device__ unsigned long long d_gfs_pair_hits;'
    source = once(source, anchor, '''static __device__ unsigned int d_gfs_seed_bin_invalid;
static __device__ unsigned int d_gfs_seed_bin_diagnostics;
static __device__ unsigned long long d_gfs_seed_bin_skips;

/* Stream-ordered validation of every CSR entry. A single bad bin, invalid
 * index, or underestimated global weight bound disables every skip for the
 * current tessellation, preserving the original complete enumeration. */
__global__ void gfs_validate_seed_bins(
    const double *px, const double *py, const float *pw2,
    const int *offset, const int *index, int n_total, int n_entries,
    int mx, int my, double xmin, double ymin, double cellsize, double maxw)
{
    int c=blockIdx.x*blockDim.x+threadIdx.x;
    if (c>=mx*my) return;
    int begin=offset[c], end=offset[c+1];
    if (begin<0 || end<begin || end>n_entries || !(cellsize>0.) || !isfinite(maxw)) {
        atomicExch(&d_gfs_seed_bin_invalid,1u); return;
    }
    int cx=c%mx, cy=c/mx;
    double xlo=xmin+cx*cellsize, xhi=xlo+cellsize;
    double ylo=ymin+cy*cellsize, yhi=ylo+cellsize;
    double tol=64.*DBL_EPSILON*(1.+fabs(xlo)+fabs(xhi)+fabs(ylo)+fabs(yhi));
    for (int k=begin;k<end;k++) {
        int j=index[k];
        if (j<0 || j>=n_total) { atomicExch(&d_gfs_seed_bin_invalid,1u); return; }
        double x=px[j], y=py[j], w=(double)pw2[j];
        if (!isfinite(x) || !isfinite(y) || !isfinite(w) || w>maxw ||
            (cx>0 && cx<mx-1 && (x<xlo-tol || x>xhi+tol)) ||
            (cy>0 && cy<my-1 && (y<ylo-tol || y>yhi+tol))) {
            atomicExch(&d_gfs_seed_bin_invalid,1u); return;
        }
    }
}

''' + anchor)
    old = '''            int cell_id = cx + mx * cy;
            int c_begin = cell_offset[cell_id];
            int c_end = cell_offset[cell_id + 1];
            for (int k = c_begin; k < c_end && nc >= 3; k++) {'''
    new = '''            int cell_id = cx + mx * cy;
            if (d_gfs_seed_bin_invalid==0u &&
                gfs_seed_bin_cannot_cut(xi,yi,cx,cy,mx,my,xmin,ymin,cellsize,
                                        w2i,global_w2_max,maxdist2)) {
                if (d_gfs_seed_bin_diagnostics) seed_bin_skips++;
                continue;
            }
            int c_begin = cell_offset[cell_id];
            int c_end = cell_offset[cell_id + 1];
            for (int k = c_begin; k < c_end && nc >= 3; k++) {'''
    source = once(source, old, new)
    source = once(source, '''    int    nc = 4, stencil_pad = GPU_SEED_PAD, seed_candidates = 0, box_retry = 0;
    double maxdist2;''', '''    int    nc = 4, stencil_pad = GPU_SEED_PAD, seed_candidates = 0, box_retry = 0;
    unsigned int seed_bin_skips=0;
    double maxdist2;''')
    source = once(source, '''    if (seed_candidates > GPU_MAX_NEIGH)
        atomicAdd(&d_gfs_candidate_fallback_cells, 1ull);''',
        '''    if (seed_bin_skips && d_gfs_seed_bin_diagnostics)
        atomicAdd(&d_gfs_seed_bin_skips,(unsigned long long)seed_bin_skips);
    if (seed_candidates > GPU_MAX_NEIGH)
        atomicAdd(&d_gfs_candidate_fallback_cells, 1ull);''')
    old = '''    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[0], s));
    /* Kernel A reads neighbor densities'''
    new = '''    const char *seed_cert_env=getenv("GFS_GPU_SEED_BIN_CERT");
    const int seed_cert_on=seed_cert_env && atoi(seed_cert_env)>0;
    const char *seed_diag_env=getenv("GFS_GPU_SEED_BIN_CERT_DIAG");
    const int seed_cert_diag=seed_cert_on && seed_diag_env && atoi(seed_diag_env)>0;
    const char *seed_fallback_env=getenv("GFS_GPU_SEED_BIN_CERT_TEST_FALLBACK");
    const int seed_force_fallback=seed_fallback_env && atoi(seed_fallback_env)>0;
    unsigned int seed_cert_invalid=(seed_cert_on && !seed_force_fallback)?0u:1u;
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_seed_bin_invalid,&seed_cert_invalid,
        sizeof(seed_cert_invalid),0,cudaMemcpyHostToDevice,s));
    unsigned int diag_enabled=seed_cert_diag?1u:0u;
    unsigned long long zero_skips=0;
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_seed_bin_diagnostics,&diag_enabled,
        sizeof(diag_enabled),0,cudaMemcpyHostToDevice,s));
    if (seed_cert_diag)
        CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_seed_bin_skips,&zero_skips,
            sizeof(zero_skips),0,cudaMemcpyHostToDevice,s));
    if (seed_cert_on) {
        const int n_cells=mx*my;
        gfs_validate_seed_bins<<<(n_cells+255)/256,256,0,s>>>(
            ctx->d_parts.x,ctx->d_parts.y,ctx->d_parts.w2,
            ctx->d_cells.cell_offset,ctx->d_cells.particle_index,n_total,
            ctx->d_cells.n_entries,mx,my,xmin,ymin,cellsize,global_w2_max);
        CUDA_CHECK(cudaGetLastError());
    }
    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[0], s));
    /* Kernel A reads neighbor densities'''
    source = once(source, old, new)
    old = '''        CUDA_CHECK(cudaStreamSynchronize(s));
        unsigned long long h_candidate_fallbacks = 0;'''
    new = '''        CUDA_CHECK(cudaStreamSynchronize(s));
        if (seed_cert_diag) {
            unsigned int invalid=0; unsigned long long skips=0;
            CUDA_CHECK(cudaMemcpyFromSymbol(&invalid,d_gfs_seed_bin_invalid,sizeof(invalid)));
            CUDA_CHECK(cudaMemcpyFromSymbol(&skips,d_gfs_seed_bin_skips,sizeof(skips)));
            fprintf(stderr,"[GPU_SEED_BIN_CERT] invalid=%u skipped_bins=%llu forced=%d\\n",
                    invalid,skips,seed_force_fallback);
        }
        unsigned long long h_candidate_fallbacks = 0;'''
    source = once(source, old, new)
    shutil.copy2(Path(__file__).resolve().parent/HEADER, header)
    cuda.write_text(source)
    print(json.dumps({'parent_cuda_sha256': EXPECTED_CUDA,
                      'patched_cuda_sha256': digest(cuda), 'header_sha256': digest(header),
                      'mode': 'opt-in; invalid CSR or weight bound falls back to full seed scan'},
                     indent=2))


if __name__ == '__main__':
    main()

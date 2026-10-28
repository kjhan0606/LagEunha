"""Build-only experiment: exact raw clipping with one per-thread polygon."""

import argparse
import hashlib
import json
from pathlib import Path


PARENT_SHA = '188c0c064deec4c21bcdf51c1b67377e58ed38aaf462fe9f241602f90b24f926'
INSERT_BEFORE = '__global__ __launch_bounds__(256, 2)\nvoid voronoi_tessellate_kernel('
INPLACE = r'''
/* Convex half-plane clipping for the raw exact path. The input and output
 * share one polygon buffer. For a convex input, all retained vertices form
 * one cyclic interval; rotate it to the front, then insert its intersections.
 * This isolated experiment keeps the old clipper for non-raw modes. */
__device__ __forceinline__
void gpu_poly_reverse(double *x, double *y, int *support, int lo, int hi)
{
    while (lo < hi) {
        double tx=x[lo],ty=y[lo]; int ts=support[lo];
        x[lo]=x[hi];y[lo]=y[hi];support[lo]=support[hi];
        x[hi]=tx;y[hi]=ty;support[hi]=ts;
        lo++;hi--;
    }
}

__device__ __forceinline__
int gpu_clip_voronoi_halfplane_inplace(
    double nx, double ny, double wfrac, int neighbor,
    double *x, double *y, int *support, int nc)
{
    const double threshold=wfrac*(nx*nx+ny*ny);
    int ninside=0, enter=-1, leave=-1;
    for (int e=0;e<nc;e++) {
        int en=(e+1)%nc;
        bool ain=(x[e]*nx+y[e]*ny)<=threshold;
        bool bin=(x[en]*nx+y[en]*ny)<=threshold;
        ninside+=(int)ain;
        if(!ain && bin) enter=e;
        if(ain && !bin) leave=e;
    }
    if(ninside==nc) return nc;
    if(ninside==0) return 0;
    int first=(enter+1)%nc;
    if(enter<0 || leave<0 || (first+ninside-1)%nc!=leave){
        atomicAdd(&d_gfs_tess_geometry_errors,1ull);
        return -1;
    }
    int first_out=(enter+1)%nc, last_next=(leave+1)%nc;
    double da=x[enter]*nx+y[enter]*ny;
    double db=x[first_out]*nx+y[first_out]*ny;
    double dc=x[leave]*nx+y[leave]*ny;
    double dd=x[last_next]*nx+y[last_next]*ny;
    double de=db-da, df=dd-dc;
    if(de==0.0 || df==0.0 || !isfinite(de) || !isfinite(df)){
        atomicAdd(&d_gfs_tess_geometry_errors,1ull);
        return -1;
    }
    double te=fmax(0.0,fmin(1.0,(threshold-da)/de));
    double tf=fmax(0.0,fmin(1.0,(threshold-dc)/df));
    double xe=x[enter]+te*(x[first_out]-x[enter]);
    double ye=y[enter]+te*(y[first_out]-y[enter]);
    int se=support[enter];
    double xf=x[leave]+tf*(x[last_next]-x[leave]);
    double yf=y[leave]+tf*(y[last_next]-y[leave]);
    bool same_enter=(xe==x[first_out] && ye==y[first_out]);
    bool same_exit=(xf==x[leave] && yf==y[leave]);
    int out_count=ninside+(same_enter?0:1)+(same_exit?0:1);
    if(out_count>GPU_MAX_CORNERS){
        unsigned long long errors=atomicAdd(&d_gfs_tess_geometry_errors,1ull);
        if(errors<8)
            printf("[GPU_CLIP_CAPACITY] nc=%d capacity=%d neighbor_soa=%d normal=(%.17g,%.17g) threshold=%.17g\n",
                   nc,GPU_MAX_CORNERS,neighbor,nx,ny,threshold);
        return -1;
    }
    /* Rotate the first retained vertex to slot zero with constant scratch. */
    gpu_poly_reverse(x,y,support,0,first-1);
    gpu_poly_reverse(x,y,support,first,nc-1);
    gpu_poly_reverse(x,y,support,0,nc-1);
    int offset=same_enter?0:1;
    if(offset){
        for(int e=ninside-1;e>=0;e--){
            x[e+1]=x[e];y[e+1]=y[e];support[e+1]=support[e];
        }
        x[0]=xe;y[0]=ye;support[0]=se;
    }
    int used=ninside+offset;
    if(same_exit) support[used-1]=neighbor;
    else {x[used]=xf;y[used]=yf;support[used]=neighbor;used++;}
    /* Match the old exact clipper's duplicate support replacement. */
    int unique=0;
    for(int e=0;e<used;e++){
        double xx=x[e],yy=y[e];int ss=support[e];
        if(unique>0 && x[unique-1]==xx && y[unique-1]==yy)
            support[unique-1]=ss;
        else {x[unique]=xx;y[unique]=yy;support[unique]=ss;unique++;}
    }
    if(unique>1 && x[0]==x[unique-1] && y[0]==y[unique-1]) unique--;
    return unique;
}

'''


def replace_one(source, old, new, name):
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
    actual = hashlib.sha256(raw).hexdigest()
    if actual != PARENT_SHA:
        raise ValueError(f'Unexpected parent CUDA SHA-256 {actual}')
    source = raw.decode()
    source = replace_one(source, INSERT_BEFORE, INPLACE +
                         'template<int EXACT_RAW>\n' + INSERT_BEFORE,
                         'in-place function and kernel specialization')
    old = '''nc=gpu_clip_voronoi_halfplane(dx,dy,wf,j,poly_x,poly_y,poly_neigh,nc,
                                         tmp_x,tmp_y,tmp_neigh,1);'''
    new = '''nc=EXACT_RAW
                ? gpu_clip_voronoi_halfplane_inplace(dx,dy,wf,j,
                        poly_x,poly_y,poly_neigh,nc)
                : gpu_clip_voronoi_halfplane(dx,dy,wf,j,
                        poly_x,poly_y,poly_neigh,nc,
                        tmp_x,tmp_y,tmp_neigh,1);'''
    source = replace_one(source, old, new, 'directional seed clip')
    old = '''nc = gpu_clip_voronoi_halfplane(dx, dy_val, wfrac, j,
                                                poly_x, poly_y, poly_neigh, nc,
                                                tmp_x, tmp_y, tmp_neigh, raw_volume >= 2);'''
    new = '''nc = EXACT_RAW
                    ? gpu_clip_voronoi_halfplane_inplace(dx, dy_val, wfrac, j,
                          poly_x, poly_y, poly_neigh, nc)
                    : gpu_clip_voronoi_halfplane(dx, dy_val, wfrac, j,
                          poly_x, poly_y, poly_neigh, nc,
                          tmp_x, tmp_y, tmp_neigh, raw_volume >= 2);'''
    source = replace_one(source, old, new, 'main seed clip')
    old = '''nc = gpu_clip_voronoi_halfplane(dx, dy_val, wfrac, j,
                                                        poly_x, poly_y, poly_neigh, nc,
                                                        tmp_x, tmp_y, tmp_neigh, raw_volume >= 2);'''
    new = '''nc = EXACT_RAW
                            ? gpu_clip_voronoi_halfplane_inplace(dx, dy_val, wfrac, j,
                                  poly_x, poly_y, poly_neigh, nc)
                            : gpu_clip_voronoi_halfplane(dx, dy_val, wfrac, j,
                                  poly_x, poly_y, poly_neigh, nc,
                                  tmp_x, tmp_y, tmp_neigh, raw_volume >= 2);'''
    source = replace_one(source, old, new, 'adaptive clip')
    source = replace_one(source,
                         'if (!raw_volume && merge_small_faces) {',
                         'if (!EXACT_RAW && !raw_volume && merge_small_faces) {',
                         'legacy temporary polygon merge')
    old = 'voronoi_tessellate_kernel<<<gridSize, blockSize, 0, s>>>(\n'
    if source.count(old) != 2:
        raise ValueError('Expected two kernel launches')
    source = source.replace(old,
        '''if (raw_geometry >= 2) {
        voronoi_tessellate_kernel<1><<<gridSize, blockSize, 0, s>>>(
''', 1)
    # Duplicate the first launch argument block for the legacy specialization.
    first = source.index('        voronoi_tessellate_kernel<1><<<')
    end = source.index('ctx->d_parts.E_inv_xy);', first) + len('ctx->d_parts.E_inv_xy);')
    launch = source[first:end]
    source = source[:end] + '\n    } else {\n' + launch.replace('kernel<1>', 'kernel<0>') + '\n    }' + source[end:]
    source = source.replace(old,
        '''if (raw_geometry >= 2) {
        voronoi_tessellate_kernel<1><<<gridSize, blockSize, 0, s>>>(
''', 1)
    second = source.index('        voronoi_tessellate_kernel<1><<<', end + len(launch))
    end2 = source.index('ctx->d_parts.E_inv_xy);', second) + len('ctx->d_parts.E_inv_xy);')
    launch2 = source[second:end2]
    source = source[:end2] + '\n    } else {\n' + launch2.replace('kernel<1>', 'kernel<0>') + '\n    }' + source[end2:]
    if source.count('voronoi_tessellate_kernel<1><<<') != 2 or source.count('voronoi_tessellate_kernel<0><<<') != 2:
        raise ValueError('Kernel specializations were not installed twice')
    path.write_text(source)
    result = {
        'parent_cuda_sha256': actual,
        'patched_cuda_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'change': 'Experimental exact raw convex clipping in one polygon buffer; legacy non-raw template retained.',
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

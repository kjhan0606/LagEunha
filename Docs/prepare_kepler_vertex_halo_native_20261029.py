"""Freeze an opt-in native W=0 certificate and read-only force diagnostic."""
import difflib
import hashlib
import json
from pathlib import Path

DOCS = Path(__file__).resolve().parent
REPO = Path('/gpfs/kjhan/LagForce/Git_publish_20261028')
OUT = DOCS / 'KEPLER_VERTEX_HALO_NATIVE_SOURCE_2026-10-29'


def replace(data, old, new, count=1):
    old, new = old.encode('ascii'), new.encode('ascii')
    if data.count(old) != count:
        raise ValueError(f'Expected {count} occurrences: {old[:100]!r}, found {data.count(old)}')
    return data.replace(old, new)


def main():
    OUT.mkdir()
    (OUT / 'Exam').mkdir()
    before = {name: (REPO / 'Exam' / name).read_bytes() for name in
              ['exam.c', 'exam_gpu.cu', 'exam_gpu.h', 'exam_gpu_extract.c']}
    frozen = REPO / 'Docs/KEPLER_ACTIVE_NATIVE_SOURCE_2026-10-28'
    for name in ['exam.c', 'exam_gpu.cu']:
        assert before[name] == (frozen / 'active_clip/Exam' / name).read_bytes()
    assert before['exam_gpu.h'] == (frozen / 'baseline/Exam/exam_gpu.h').read_bytes()
    assert before['exam_gpu_extract.c'] == (frozen / 'padding_reuse/exam_gpu_extract.c').read_bytes()
    after = dict(before)
    helper = '''static int gfs_kepler_vertex_halo_requested(SimParameters *simpar,
    void (*padding)(SimParameters *,postype)){
    const char *setting=getenv("GFS_GPU_VERTEX_HALO_CERT");
    const char *ic=getenv("EUNHA_IC");
    /* Only this callback supplies the complete periodic rectangle used in
     * the initial qualification. MPI and weighted paths retain the old test. */
    return setting && atoi(setting)>0 && ic && strcmp(ic,"kepler")==0 &&
        NID(simpar)==1 && GAS_Kappa(simpar)==0 &&
        padding==paddingTreeVorork4Particles;
}

'''
    after['exam.c'] = replace(after['exam.c'], 'postype gfs_kepler_halo_width(postype cellsize){',
                              helper + 'postype gfs_kepler_halo_width(postype cellsize){')
    after['exam.c'] = replace(after['exam.c'],
        'GAS_ENTROPY_MODE(simpar), GAS_K_FLOOR(simpar));',
        'GAS_ENTROPY_MODE(simpar), GAS_K_FLOOR(simpar),\n'
        '\t\t\tgfs_kepler_vertex_halo_requested(simpar,paddingAllTreeParticles) && all_weights_zero, xmax, ymax);')
    after['exam.c'] = replace(after['exam.c'],
        '                                  boxsize, VoroAccuracyOrder(simpar));',
        '                                  boxsize, VoroAccuracyOrder(simpar), xmax, ymax,\n'
        '                                  gfs_kepler_vertex_halo_requested(simpar,paddingAllTreeParticles));')
    for name in ['exam_gpu.h', 'exam_gpu.cu']:
        after[name] = replace(after[name], '    int entropy_mode, double K_floor)' + (';' if name.endswith('.h') else ''),
            '    int entropy_mode, double K_floor, int vertex_halo, double complete_xmax, double complete_ymax)' + (';' if name.endswith('.h') else ''))
        # Restrict the replacement to the complete raw-volume prototype/function.
        old = '''    double boxsize, double global_w2_max, double OoA, int av_mode,
    int gradient_method)'''
        new = '''    double boxsize, double global_w2_max, double OoA, int av_mode,
    int gradient_method, int vertex_halo, double complete_xmax, double complete_ymax)'''
        after[name] = replace(after[name], old, new)
    for name in ['exam_gpu.h', 'exam_gpu_extract.c']:
        after[name] = replace(after[name],
            '    double cellsize, double xmin, double ymin, double boxsize, double OoA)',
            '    double cellsize, double xmin, double ymin, double boxsize, double OoA,\n'
            '    double complete_xmax, double complete_ymax, int vertex_halo)')
    after['exam_gpu_extract.c'] = replace(after['exam_gpu_extract.c'],
        '    int failed = gpu_tessellate_raw_volume(&g_gpu_ctx, nbp, mx, my,',
        '''    if(vertex_halo) for(int k=0;k<n_total;k++)
        if(g_gpu_ctx.h_parts.w2[k]!=0) {vertex_halo=0;break;}
    int failed = gpu_tessellate_raw_volume(&g_gpu_ctx, nbp, mx, my,''')
    after['exam_gpu_extract.c'] = replace(after['exam_gpu_extract.c'],
        '        GAS_GRADIENT_METHOD(simpar));\n    for (int k = 0; k < nbp; k++)',
        '        GAS_GRADIENT_METHOD(simpar), vertex_halo, complete_xmax, complete_ymax);\n    for (int k = 0; k < nbp; k++)')
    dump = (DOCS / 'gfs_gpu_geometry_force_dump_20261029.inc').read_bytes()
    after['exam_gpu_extract.c'] = replace(after['exam_gpu_extract.c'],
        'static int gpu_call_count = 0;\n', 'static int gpu_call_count = 0;\n')
    after['exam_gpu_extract.c'] = after['exam_gpu_extract.c'].replace(
        b'static int gpu_call_count = 0;\n', b'static int gpu_call_count = 0;\n\n' + dump + b'\n', 1)
    after['exam_gpu_extract.c'] = replace(after['exam_gpu_extract.c'],
        '    double t7_done = MPI_Wtime();\n',
        '    double t7_done = MPI_Wtime();\n'
        '    gfs_dump_geometry_force(simpar,&g_gpu_ctx,nbp,n_total,n_faces,gpu_call_count+1);\n')
    cu = after['exam_gpu.cu']
    cu = replace(cu, '#include "gfs_pair.h"', '#include "gfs_pair.h"\n#include "gfs_vertex_halo.h"')
    cu = replace(cu, 'static __device__ unsigned long long d_gfs_stencil_halo_clipped_cells;',
        'static __device__ unsigned long long d_gfs_stencil_halo_clipped_cells;\n'
        'static __device__ unsigned long long d_gfs_vertex_halo_recertified_cells;')
    cu = replace(cu, '    int thermal64, double *__restrict__ out_V64)',
        '    int thermal64, double *__restrict__ out_V64,\n'
        '    int vertex_halo, double complete_xmax, double complete_ymax)')
    cu = replace(cu, '        thermal64, ctx->d_parts.E_inv_xy);',
        '        thermal64, ctx->d_parts.E_inv_xy, vertex_halo, complete_xmax, complete_ymax);', count=2)
    cu = replace(cu, '        if (po->halo_clipped)\n            atomicAdd(&d_gfs_stencil_halo_clipped_cells, 1ull);',
        '''        if(vertex_halo){
            int old_clipped=po->halo_clipped, complete=1;
            /* Use actual padding bounds, never rounded-up grid allocation.
             * Strict directed bounds protect the represented vertices;
             * search, clipping, gradients and force arithmetic are unchanged. */
            for(int e=0;e<nc && complete;e++)
                complete=gfs_vertex_halo_disk_inside(xi,yi,poly_x[e],poly_y[e],
                    xmin,ymin,complete_xmax,complete_ymax);
            po->halo_clipped=!complete;
            if(old_clipped && complete)
                atomicAdd(&d_gfs_vertex_halo_recertified_cells,1ull);
        }
        if (po->halo_clipped)
            atomicAdd(&d_gfs_stencil_halo_clipped_cells, 1ull);''')
    # Reset and publish a count independently of the unchanged old guard.
    cu = replace(cu, '    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[0], s));',
        '''    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_vertex_halo_recertified_cells,&zero_fallbacks,
        sizeof(zero_fallbacks),0,cudaMemcpyHostToDevice,s));
    if(vertex_halo){
        static int banner=0;
        if(!banner++) fprintf(stderr,"[VERTEX_HALO] force=1 actual_rectangle=[%.17g,%.17g]x[%.17g,%.17g] W0_only=1\\n",
            xmin,complete_xmax,ymin,complete_ymax);
    }
    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[0], s));''')
    cu = replace(cu, '        if (h_stencil_expanded > 0 || h_halo_clipped > 0)',
        '''        if(vertex_halo){
            unsigned long long recertified=0;
            CUDA_CHECK(cudaMemcpyFromSymbol(&recertified,d_gfs_vertex_halo_recertified_cells,sizeof(recertified)));
            fprintf(stderr,"[VERTEX_HALO] recertified=%llu incomplete=%llu\\n",recertified,h_halo_clipped);
        }
        if (h_stencil_expanded > 0 || h_halo_clipped > 0)''')
    cu = replace(cu, '    int raw_geometry=(raw_env && atoi(raw_env)>0) ? 2 : 1;',
        '''    int raw_geometry=(raw_env && atoi(raw_env)>0) ? 2 : 1;
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_vertex_halo_recertified_cells,&zero,sizeof(zero),0,cudaMemcpyHostToDevice,s));
    if(vertex_halo){
        static int banner=0;
        if(!banner++) fprintf(stderr,"[VERTEX_HALO] final_volume=1 actual_rectangle=[%.17g,%.17g]x[%.17g,%.17g] W0_only=1\\n",
            xmin,complete_xmax,ymin,complete_ymax);
    }''')
    after['exam_gpu.cu'] = cu
    manifest = {'base_repository': str(REPO), 'source_basis': 'Frozen accepted clip source; existing default-off padding-reuse wrapper retained',
                'enabled_by': 'GFS_GPU_VERTEX_HALO_CERT=1; only one MPI rank, Kepler IC, kappa=0, exact W=0 and the standard periodic-padding callback',
                'physics_changed': False, 'production_promoted': False,
                'native_parity_required': True, 'files': {}}
    patches = []
    for name, data in after.items():
        (OUT / 'Exam' / name).write_bytes(data)
        manifest['files'][name] = {'before_sha256': hashlib.sha256(before[name]).hexdigest(),
                                  'after_sha256': hashlib.sha256(data).hexdigest()}
        patches.append(''.join(difflib.unified_diff(before[name].decode('latin1').splitlines(True),
            data.decode('latin1').splitlines(True), fromfile='a/Exam/'+name,tofile='b/Exam/'+name)))
    header = (DOCS / 'gfs_vertex_halo.h').read_bytes()
    (OUT / 'Exam/gfs_vertex_halo.h').write_bytes(header)
    manifest['files']['gfs_vertex_halo.h'] = {'after_sha256': hashlib.sha256(header).hexdigest()}
    (OUT / 'native.patch').write_bytes(''.join(patches).encode('latin1'))
    (OUT / 'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))


if __name__ == '__main__':
    main()

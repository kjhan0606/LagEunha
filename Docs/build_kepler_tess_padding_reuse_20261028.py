"""Build an opt-in host-wrapper optimization from the proven native objects."""
import difflib
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import socket
import subprocess

BASE = Path('/tmp/Kepler_OMP_RHS_static_corrected_20261028')
CLIP = Path('/tmp/clip_swap_build_20261028')
WORK = Path('/tmp/Kepler_tess_padding_reuse_20261028')
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/tess_padding_reuse_pair_20261028')
EXPECTED = '6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args, cwd):
    result = subprocess.run(args, cwd=cwd, check=True, capture_output=True, text=True)
    return result.stdout


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'Expected exactly one patch location: {old[:80]}')
    return text.replace(old, new, 1)


def main():
    if socket.gethostname().split('.')[0] != 'syntax':
        raise SystemExit('Run on syntax, where the verified original build lives')
    assert digest(CLIP / 'eunha2_clip_swap.exe') == EXPECTED
    assert digest(BASE / 'eunha2.exe') == 'a6f16e8b1e143f70e8b8c32155008077bd84e813f7f4e4aff2f04b23dfd04572'
    original = BASE / 'Exam/exam_gpu_extract.c'
    assert digest(original) == '3f50f417275f2778758c920f513d4b61b87cbb651e4825b516febd85a44dc365'
    assert digest(BASE / 'Exam/exam_gpu.h') == '45ceb4cb8f9b153ef8209da7c65cba53468d5f0f9cb1586803de0b5c0c6b7eb0'
    if ROOT.exists() or WORK.exists():
        raise FileExistsError('Refusing to overwrite an existing candidate workspace')
    ROOT.mkdir()
    (ROOT / 'source_overlay').mkdir()
    WORK.mkdir()
    (WORK / 'Exam').mkdir()
    (WORK / 'objects').mkdir()
    # Mirror only two immediate source-directory levels through symlinks.
    # The compiler and linker write exclusively to the new objects directory.
    for path in BASE.iterdir():
        if path.name not in ('Exam', 'objects'):
            (WORK / path.name).symlink_to(path)
    for path in (BASE / 'Exam').iterdir():
        if path.name != 'exam_gpu_extract.c':
            (WORK / 'Exam' / path.name).symlink_to(path)

    before = original.read_text()
    after = replace_once(before,
        '    int fused = gpu_fused_ready();\n',
        '''    int fused = gpu_fused_ready();
    int tess_on_device = g_gpu_ctx.initialized && g_gpu_ctx.tess_faces_on_device;
    const char *reuse_env = getenv("GFS_REUSE_TESS_PADDING");
    const char *sync_env_reuse = getenv("GFS_GPU_SYNC_HALO_STATE");
    /* The cached device CSR indexes the tessellation's padding. Reuse it
     * only when the existing post-gradient halo refresh remains enabled;
     * that refresh verifies every image ID and position before upload. */
    int reuse_tess_padding = tess_on_device && reuse_env && atoi(reuse_env)>0
        && (!sync_env_reuse || atoi(sync_env_reuse)>0);
''')
    after = replace_once(after,
        '    if (fused) {\n        /* Fused path: cells already freed by updateDenW2Pressure2DBlend,\n',
        '''    if (fused || reuse_tess_padding) {
        /* GPU tessellation already freed the linked-cell array and retained
         * the exact padding indexed by its device face CSR. */
        /* Fused path: cells already freed by updateDenW2Pressure2DBlend,
''')
    after = replace_once(after,
        '    int tess_on_device = g_gpu_ctx.tess_faces_on_device;\n',
        '''    if (reuse_tess_padding && !g_gpu_ctx.tess_faces_on_device) {
        fprintf(stderr,"[TESS_PADDING_REUSE] unexpected context rebuild; refusing stale CSR\\n");
        MPI_Abort(MPI_COMM(simpar), 98);
    }
    tess_on_device = g_gpu_ctx.tess_faces_on_device;
    if (reuse_tess_padding) {
        static int reuse_banner = 0;
        if (!reuse_banner++) fprintf(stderr,
            "[TESS_PADDING_REUSE] rank=%d np=%d npad=%d reuse=1 refresh_validation=1\\n",
            MYID(simpar), nbp, npad);
    }
''')
    after = replace_once(after,
        '    if (!fused)\n        my_free(VORO_BASICCELL(simpar));\n',
        '    if (!fused && !reuse_tess_padding)\n        my_free(VORO_BASICCELL(simpar));\n')
    candidate_source = WORK / 'Exam/exam_gpu_extract.c'
    candidate_source.write_text(after)
    (ROOT / 'source_overlay/exam_gpu_extract.original.c').write_text(before)
    (ROOT / 'source_overlay/exam_gpu_extract.c').write_text(after)
    shutil.copy2(BASE / 'Exam/exam_gpu.h', ROOT / 'source_overlay/exam_gpu.h')
    patch = ''.join(difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                    fromfile='Exam/exam_gpu_extract.c', tofile='Exam/exam_gpu_extract.c'))
    (ROOT / 'source_overlay/reuse_padding.patch').write_text(patch)

    flags = ['USE_CUDA=1', 'GFS_VARPRESS_W0=1']
    cc_text = run(['make', '--no-print-directory', '-n', '-W', 'exam_gpu_extract.c',
                   'exam_gpu_extract.o', *flags], BASE / 'Exam')
    cc = shlex.split(next(line for line in cc_text.splitlines() if line.startswith('mpiicx ')))
    cc[cc.index('exam_gpu_extract.c')] = str(candidate_source)
    cc[cc.index('-o') + 1] = str(WORK / 'objects/exam_gpu_extract.o')
    ld_text = run(['make', '--no-print-directory', '-n', '-B', 'eunha2.exe', *flags], BASE)
    ld = shlex.split(next(line for line in ld_text.splitlines() if line.startswith('mpiifx ')))
    assert ld.count('-L./Exam') == 1
    ld[ld.index('-L./Exam')] = '-L' + str(WORK / 'objects')
    library = WORK / 'objects/libexam.a'
    shutil.copy2(CLIP / 'libexam.a', library)
    ld[ld.index('-o') + 1] = str(WORK / 'objects/eunha2.control_relink.exe')
    (ROOT / 'control_link.log').write_text(run(ld, BASE))
    control_hash = digest(WORK / 'objects/eunha2.control_relink.exe')
    # Preserve the old GPU and RK4 objects exactly; rebuild one C wrapper.
    cc_result = subprocess.run(cc, cwd=BASE / 'Exam', text=True, capture_output=True)
    (ROOT / 'compile.log').write_text(cc_result.stdout + cc_result.stderr)
    cc_result.check_returncode()
    run(['ar', 'rcs', str(library), str(WORK / 'objects/exam_gpu_extract.o')], WORK)
    for member in ('exam_gpu.o', 'exam.o'):
        old = subprocess.run(['ar', 'p', str(CLIP / 'libexam.a'), member], check=True, capture_output=True).stdout
        new = subprocess.run(['ar', 'p', str(library), member], check=True, capture_output=True).stdout
        assert old == new, f'Unexpected change in {member}'
    ld[ld.index('-o') + 1] = str(WORK / 'objects/eunha2.reuse.exe')
    (ROOT / 'candidate_link.log').write_text(run(ld, BASE))
    for name in ('eunha2.control_relink.exe', 'eunha2.reuse.exe'):
        shutil.copy2(WORK / 'objects' / name, ROOT / name)
    receipt = {
        'purpose': 'Opt-in reuse of padding indexed by an already valid GPU tessellation CSR',
        'baseline_project': str(BASE), 'baseline_binary_sha256': EXPECTED,
        'control_relink_sha256': control_hash, 'control_relink_bytewise_equal_to_proven_binary': control_hash == EXPECTED,
        'candidate_sha256': digest(ROOT / 'eunha2.reuse.exe'),
        'original_wrapper_sha256': digest(original), 'candidate_wrapper_sha256': digest(candidate_source),
        'compile_command': cc, 'link_command': ld,
        'GPU_and_RK4_object_bytes_unchanged': True,
        'enabled_by': 'GFS_REUSE_TESS_PADDING=1; default remains off',
        'retained_safety': 'Full gradient-state refresh; every halo ID and image position checked before force upload',
        'runtime_validation': 'Not yet performed; keep existing long-time trajectory unchanged',
        'work': str(WORK), 'root': str(ROOT),
    }
    (ROOT / 'build_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()

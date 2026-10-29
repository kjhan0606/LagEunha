"""Isolated native build: retain baseline potential diagnostic, cache thermal mode only."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(root):
    source=Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/project')
    base=source.parent/'candidate_v2/Exam/exam.c'
    header=source.parent/'candidate_const_cache_v1/Exam/gfs_thermal64.h'
    if sha(base)!='328982249341cc46b365eb89850683d575ebe54029e453159e73a35297b3ece7':
        raise ValueError('Original wall-v2 host source changed')
    if sha(header)!='d51aca8f82aa6345e181ff37d574e8224e206fe364d0c71bf6efe861537950c9':
        raise ValueError('Qualified thermal-cache header changed')
    frozen={'exam_gpu.cu':'09a2d6a16d44f03534d38ea6fdfbaa556736d62fad7504981258a11ffbf6f5a3',
            'gfs_wall_events.inc':'719b59918cdfad5e8a041eb32090aefe523c8a208039da8f6ea559ed4f3bd450',
            'gfs_wall_dense.h':'01ca236d0eb731b43c26679b84509bef3e92cfc5872adadb4f304dfba5adb323',
            'kepler_slip_wall.h':'4032460b4631bb430dea0c61d36de51118918fa120efb826242b38fbf942cc25'}
    for name,expected in frozen.items():
        if sha(source/'Exam'/name)!=expected:raise ValueError('Frozen native source changed: '+name)
    dest=root/'project'
    if dest.exists():raise FileExistsError('Refusing to overwrite any build')
    dest.mkdir()
    # Existing libraries/source are read-only symlinks. Objects and output
    # that the build can mutate are independent copies under this directory.
    for p in source.iterdir():
        if p.name in ('Exam','eunha2.exe'):continue
        target=dest/p.name
        if p.name=='eunha2.o':shutil.copy2(p,target)
        else:target.symlink_to(p,target_is_directory=p.is_dir())
    exam=dest/'Exam';exam.mkdir()
    for p in (source/'Exam').iterdir():
        target=exam/p.name
        if p.is_dir():target.symlink_to(p,target_is_directory=True)
        elif p.is_file():shutil.copy2(p,target)
    shutil.copyfile(base,exam/'exam.c');shutil.copyfile(header,exam/'gfs_thermal64.h')
    os.utime(exam/'exam.c',None);os.utime(exam/'gfs_thermal64.h',None)
    assert sha(exam/'exam.c')==sha(base) and sha(exam/'gfs_thermal64.h')==sha(header)
    for name,cmd in [('exam',['make','-C','Exam','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','libexam.a']),
                     ('link',['make','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','this'])]:
        with (root/('build_'+name+'.log')).open('x') as log:
            subprocess.run(cmd,cwd=dest,stdout=log,stderr=subprocess.STDOUT,check=True)
    exe=dest/'eunha2.exe'
    result={'passed':True,'host_source_sha256':sha(exam/'exam.c'),'thermal_header_sha256':sha(exam/'gfs_thermal64.h'),
            'gpu_source_sha256':sha(exam/'exam_gpu.cu'),'candidate_binary':str(exe),'candidate_binary_sha256':sha(exe),
            'source_project':str(source),'baseline_native_binary_sha256':'f26816122df101ea71f72b4f3fb9233e965d06784ed94e88cfb46245fa251444',
            'scope':'Separate build with original v2 host potential code and opt-in process-constant thermal-mode cache. Existing projects and running executables untouched. No GPU parity or speed qualification has run.'}
    with (root/'build_result.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':build(Path(sys.argv[1]).resolve())

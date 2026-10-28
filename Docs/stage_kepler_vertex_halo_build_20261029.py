"""Stage declared native build dependencies for a CPU-only Slurm compile."""
import datetime
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import socket
import subprocess

DOCS = Path(__file__).resolve().parent
BASE = Path('/tmp/Kepler_OMP_RHS_static_corrected_20261028')
CLIP = Path('/tmp/clip_swap_build_20261028')
MIRROR = Path('/tmp/Kepler_vertex_halo_dependencies_20261029')
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/vertex_halo_native_build_20261029')


def run(args, cwd):
    result = subprocess.run(args, cwd=cwd, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f'{args}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert socket.gethostname().split('.')[0] == 'syntax'
    assert digest(CLIP/'eunha2_clip_swap.exe') == '6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0'
    assert not (ROOT/'job.slurm').exists(), 'Prepared/submitted build exists; inspect it instead of submitting again'
    assert not (ROOT/'eunha2.vertex.exe').exists()
    ROOT.mkdir(exist_ok=True)
    project = ROOT/'project'
    project.mkdir(exist_ok=True)
    (project/'Exam').mkdir(exist_ok=True)
    (project/'objects').mkdir(exist_ok=True)
    MIRROR.mkdir(exist_ok=True)
    (MIRROR/'Exam').mkdir(exist_ok=True)
    source = DOCS/'KEPLER_VERTEX_HALO_NATIVE_SOURCE_2026-10-29/Exam'
    changed = {p.name for p in source.iterdir()}
    for p in BASE.iterdir():
        if p.name != 'Exam':
            link=MIRROR/p.name
            if link.is_symlink(): assert link.resolve()==p.resolve()
            else: link.symlink_to(p)
    for p in (BASE/'Exam').iterdir():
        if p.name not in changed:
            link=MIRROR/'Exam'/p.name
            if link.is_symlink(): assert link.resolve()==p.resolve()
            else: link.symlink_to(p)
    for p in source.iterdir(): shutil.copy2(p, MIRROR/'Exam'/p.name)
    flags = ['USE_CUDA=1','GFS_VARPRESS_W0=1']
    commands=[]
    dependencies={}
    for filename, compiler in [('exam.c','mpiicx'),('exam_gpu.cu','/opt/ohpc/pub/cuda/13.0.2/bin/nvcc'),('exam_gpu_extract.c','mpiicx')]:
        target=filename.rsplit('.',1)[0]+'.o'
        dry=run(['make','--no-print-directory','-n','-W',filename,target,*flags],BASE/'Exam')
        cc=shlex.split(next(line for line in dry.splitlines() if line.startswith(compiler+' ')))
        cc[cc.index(filename)]=str(MIRROR/'Exam'/filename)
        if '-o' in cc:
            cc[cc.index('-o')+1]=str(project/'objects'/target)
        else: cc += ['-o',str(project/'objects'/target)]
        dep=cc[:]
        out=dep.index('-o'); del dep[out:out+2]
        dep.remove('-c'); dep.append('-M')
        listing=run(dep,MIRROR/'Exam')
        # Compiler-generated make dependencies; no recursive filesystem walk.
        target_text, separator, dependency_text=listing.replace('\\\n',' ').partition(':')
        assert separator and target_text.strip(), listing
        tokens=shlex.split(dependency_text)
        for token in tokens:
            path=Path(token)
            if not path.is_absolute(): path=MIRROR/'Exam'/path
            path=path.resolve()
            if path.is_relative_to(BASE): relative=path.relative_to(BASE)
            elif path.is_relative_to(MIRROR): relative=path.relative_to(MIRROR)
            else: continue  # shared compiler/system/external include path
            dest=project/relative
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(path,dest)
            dependencies[str(relative)]=digest(dest)
        cc=[part.replace(str(MIRROR),str(project)) for part in cc]
        commands.append({'cwd':str(project/'Exam'),'argv':cc})
    for p in source.iterdir():
        shutil.copy2(p,project/'Exam'/p.name)
        dependencies['Exam/'+p.name]=digest(p)
    dry=run(['make','--no-print-directory','-n','-B','eunha2.exe',*flags],BASE)
    ld=shlex.split(next(line for line in dry.splitlines() if line.startswith('mpiifx ')))
    libraries=[]
    for arg in ld:
        if arg.startswith('-L') and len(arg)>2:
            directory=Path(arg[2:])
            if not directory.is_absolute(): directory=BASE/directory
            libraries.append(directory.resolve())
    for arg in ld:
        if arg.startswith('-l') and arg!='-lexam':
            for directory in libraries:
                p=directory/('lib'+arg[2:]+'.a')
                if p.is_file():
                    if p.is_relative_to(BASE):
                        dest=project/p.relative_to(BASE)
                        dest.parent.mkdir(parents=True,exist_ok=True)
                        shutil.copy2(p,dest)
                        dependencies[str(dest.relative_to(project))]=digest(dest)
                    break
        elif arg.endswith('.o'):
            p=BASE/arg
            if p.is_file():
                dest=project/arg;dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(p,dest);dependencies[arg]=digest(dest)
    shutil.copy2(CLIP/'libexam.a',project/'objects/libexam.a')
    shutil.copy2(CLIP/'libexam.a',project/'objects/libexam.input.a')
    dependencies['objects/libexam.input.a']=digest(project/'objects/libexam.input.a')
    ld[ld.index('-L./Exam')]='-L'+str(project/'objects')
    ld[ld.index('-o')+1]=str(ROOT/'eunha2.vertex.exe')
    commands += [{'cwd':str(project),'argv':['ar','rcs',str(project/'objects/libexam.a'),
                    *[str(project/'objects'/x) for x in ['exam.o','exam_gpu.o','exam_gpu_extract.o']]]},
                 {'cwd':str(project),'argv':ld}]
    (ROOT/'build_commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    (ROOT/'inputs.sha256').write_text(''.join(f'{sha}  project/{name}\n' for name,sha in sorted(dependencies.items()))
                                    +f'{digest(ROOT/"build_commands.json")}  build_commands.json\n')
    runner='''import hashlib,json,subprocess,time,shutil
from pathlib import Path
root=Path(__file__).resolve().parent
commands=json.loads((root/'build_commands.json').read_text())
start=time.perf_counter()
shutil.copy2(root/'project/objects/libexam.input.a',root/'project/objects/libexam.a')
for i,item in enumerate(commands):
    result=subprocess.run(item['argv'],cwd=item['cwd'],capture_output=True)
    (root/f'command_{i}.stdout').write_bytes(result.stdout)
    (root/f'command_{i}.stderr').write_bytes(result.stderr)
    print(f"command {i} exit={result.returncode}",flush=True)
    result.check_returncode()
binary=root/'eunha2.vertex.exe'
result={'binary':str(binary),'sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'build_seconds':time.perf_counter()-start,'production_promoted':False,
        'runtime_validation':'Not performed; existing long trajectory untouched'}
(root/'build_result.json').write_text(json.dumps(result,indent=2)+'\\n')
print(json.dumps(result,indent=2))
'''
    (ROOT/'run_build.py').write_text(runner)
    job=f'''#!/usr/bin/env bash
#SBATCH --job-name=KepVertexBuild
#SBATCH --partition=a100_pcie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --chdir={ROOT}
#SBATCH --output={ROOT}/slurm-%j.out
#SBATCH --error={ROOT}/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 run_build.py
'''
    (ROOT/'job.slurm').write_text(job)
    run(['bash','-n',str(ROOT/'job.slurm')],ROOT)
    submitted=run(['sbatch','--parsable',str(ROOT/'job.slurm')],ROOT).strip().split(';')[0]
    assert submitted.isdecimal(),submitted
    result={'recorded_at':datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
            'job_id':int(submitted),'root':str(ROOT),'base':str(BASE),
            'copied_declared_dependencies':len(dependencies),'dependency_bytes':sum((project/p).stat().st_size for p in dependencies),
            'GPU_request':0,'purpose':'CPU-only compilation of an isolated native W0 vertex-certificate candidate; no simulation or production replacement',
            'source_manifest':str(DOCS/'KEPLER_VERTEX_HALO_NATIVE_SOURCE_2026-10-29/MANIFEST.json'),
            'commands':commands,'production_modified':False}
    (ROOT/'staging_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    (DOCS/'KEPLER_VERTEX_HALO_BUILD_STAGE_2026-10-29.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='commands'},indent=2))


if __name__ == '__main__':
    main()

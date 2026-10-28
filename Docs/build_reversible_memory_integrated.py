from pathlib import Path
import subprocess,shutil,json,hashlib
root=Path('/gpfs/kjhan/LagForce/Codex_reversible_memory_integrated_build_20261018')
base=Path('/gpfs/kjhan/LagForce/Codex_accepted_ledger_build_20261017/src')
c=json.loads((root/'commands.json').read_text())
with (root/'build.log').open('x') as log:
    for name in ('compile','kh_compile'):
        subprocess.run(c[name],cwd=base/'Exam',stdout=log,stderr=subprocess.STDOUT,check=True)
    shutil.copy2('/gpfs/kjhan/LagForce/Codex_integrated_cpu_cuda_20261017/libexam.a',root/'libexam.a')
    shutil.copy2('/gpfs/kjhan/LagForce/Codex_OMP_particle_parity_20261017/libkh.a',root/'libkh.a')
    for library,obj in [('libexam.a','exam.o'),('libkh.a','kh.o')]:
        subprocess.run(['ar','r',str(root/library),str(root/obj)],stdout=log,stderr=subprocess.STDOUT,check=True)
    subprocess.run(c['link'],cwd=base,stdout=log,stderr=subprocess.STDOUT,check=True)
(root/'binary.sha256').write_text(hashlib.sha256((root/'eunha2.exe').read_bytes()).hexdigest()+'\n')
(root/'build.exit').write_text('0\n')

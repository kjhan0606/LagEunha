"""Build a separate opt-in OpenMP row implementation; preserve live binaries."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
b=Path('/gpfs/kjhan/LagForce')
old=b/'Codex_reversible_memory_integrated_build_20261018'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,default=b/'Codex_memory_parallel_build_20261018')
root=parser.parse_args().root
base=b/'Codex_accepted_ledger_build_20261017/src'
commands=json.loads((old/'commands.json').read_text())
commands={key:[s.replace(str(old),str(root)) for s in value] for key,value in commands.items()}
(root/'commands.json').write_text(json.dumps(commands,indent=2))
with (root/'build.log').open('x') as log:
    subprocess.run(commands['compile'],cwd=base/'Exam',stdout=log,stderr=subprocess.STDOUT,check=True)
    for name in ('libexam.a','libkh.a'):shutil.copyfile(old/name,root/name)
    subprocess.run(['ar','r',str(root/'libexam.a'),str(root/'exam.o')],check=True,stdout=log,stderr=subprocess.STDOUT)
    subprocess.run(commands['link'],cwd=base,check=True,stdout=log,stderr=subprocess.STDOUT)
(root/'binary.sha256').write_text(hashlib.sha256((root/'eunha2.exe').read_bytes()).hexdigest()+'\n')
(root/'build.exit').write_text('0\n')

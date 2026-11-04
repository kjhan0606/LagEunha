from pathlib import Path
import json,subprocess,os,hashlib,socket
p=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1');exe=p/'smooth_native_build_v2/unified_native_general.exe'
e=json.loads((p/'smooth_preflight_v1/kh/environment.json').read_text());q=subprocess.run(['ldd',str(exe)],env={**os.environ,**e},text=True,capture_output=True)
missing=[s for s in q.stdout.splitlines() if 'not found' in s]
print(json.dumps({'host':socket.gethostname(),'project':str(p),'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'ldd_exit':q.returncode,'missing_libraries':missing}))
assert q.returncode==0 and not missing

"""Stage spin/strain postprocessing of the accepted fresh weighted t~7.425 state."""
import json
from pathlib import Path
import shutil
import subprocess
from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RESULTS=DOCS/'KEPLER_FRESH_WEIGHTED_T7425_2026-10-27'
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/fresh_spin_t7425_cpu_v1')
FILES=('kepler_fresh_spin_t7425_cpu_20261027.slurm',
       'audit_kepler_affine_spin_energy_20261027.py',
       'test_kepler_affine_spin_energy_20261027.py',
       'audit_kepler_polygon_moment_native_faces_20261027.py',
       'audit_kepler_sparse_area_t0825_20261027.py',
       'laguerre_nonlinear_volume_probe_t6_20261026.py',
       'prepare_kepler_sparse_area_restart_20261027.py')

def main():
    receipt=RESULTS/'weighted_receipt.json';data=json.loads(receipt.read_text())
    if len(data['checks'])!=21 or not all(data['checks'].values()) or abs(data['time']-7.425299506861853)>1e-12:
        raise ValueError('Weighted t≈7.425 parent gate failed')
    if RUN.exists():raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(DOCS/FILES[0])],check=True)
    RUN.mkdir();hashes={}
    for name in FILES:
        source=DOCS/name;shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Changed {name}')
    shutil.copy2(receipt,RUN/'weighted_receipt.json');hashes['weighted_receipt.json']=digest(RUN/'weighted_receipt.json')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'parent_receipt':str(receipt),'parent_receipt_sha256':digest(receipt),
                   'time':data['time'],'input_sha256':hashes,
                   'scope':'Postprocess saved polygon and velocity-gradient state; no spin variable is evolved.'},out,indent=2)
        out.write('\n')
    print(RUN)

if __name__=='__main__':main()

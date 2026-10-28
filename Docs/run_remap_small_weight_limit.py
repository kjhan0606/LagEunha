"""Reuse audited remap implementation at smaller displacements in a fresh root."""
from pathlib import Path
import hashlib,json,shutil
base=Path('/gpfs/kjhan/LagForce')
src=base/'Codex_checkpoint_weight_remap_20261018'
out=base/'Codex_remap_small_weight_limit_20261018'
out.mkdir(exist_ok=False)
for name in ('cpu.so','cuda.so','weight_remap_moments.h','probe.c','probe.cu'):
    (out/name).symlink_to(src/name)
s= (src/'audit.py').read_text()
assert s.count('for fraction in (.25,1.):')==1
s=s.replace('for fraction in (.25,1.):','for fraction in (.0625,.015625,.00390625):')
s=s.replace('for order in (4,8):','for order in (8,):')
s=s.replace('    point=np.array(',"    np.savez(root/'reference_q8.npz',ids=ids,state=reference[8],volume=volume)\n    point=np.array(")
(out/'audit.py').write_text(s)
(out/'provenance.json').write_text(json.dumps({'original_script':str(src/'audit.py'),'original_sha256':hashlib.sha256((src/'audit.py').read_bytes()).hexdigest(),'changes':['fractions .0625,.015625,.00390625','quadrature8 only','save pre-remap reference'],'scope':'Finite remap diagnostic; no physical time evolution'},indent=2)+'\n')
exec(compile(s,str(out/'audit.py'),'exec'),{'__name__':'__main__','__file__':str(out/'audit.py')})

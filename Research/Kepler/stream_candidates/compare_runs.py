from pathlib import Path
import json,re,sys
import numpy as np
from audit_snapshot import read_snapshot
ROOT=Path(__file__).resolve().parent

def compare(left,right,call):
    a,b=ROOT/left,ROOT/right
    result={'left':left,'right':right,'timing_scope':'GPU shared with production; elapsed times are not isolated speed measurements'}
    for d in (a,b):
        t=json.loads((d/'timing.json').read_text())
        if t['exit_code']!=0: raise ValueError(f'{d} did not complete')
        result[d.name+'_timing']=t
    def snapshot(d):
        files=[p for p in d.glob('khout.*.dat') if not p.is_symlink()]
        if len(files)!=1: raise ValueError(f'{d}: expected one output snapshot, got {files}')
        data,meta=read_snapshot(files[0]); ids=data['u4if']&((1<<56)-1)
        order=np.argsort(ids)
        if len(np.unique(ids))!=len(ids): raise ValueError('Duplicate particle IDs')
        return data[order],meta,ids[order]
    pa,ma,ia=snapshot(a); pb,mb,ib=snapshot(b)
    result['times']=[ma['t'],mb['t']]
    result['same_particle_ids']=bool(np.array_equal(ia,ib))
    result['state']={}
    for field in ('x','y','vx','vy','mass','volume','den','pressure','ie','w2','stress.gUxx','stress.gUxy','stress.gUyx','stress.gUyy','stress.dRhodx'):
        x,y=pa[field].astype(float),pb[field].astype(float)
        result['state'][field]={'max_abs':float(np.max(abs(x-y))), 'relative_l2':float(np.linalg.norm(x-y)/max(np.linalg.norm(x),1e-300)), 'identical_count':int((x==y).sum()), 'finite':bool(np.isfinite(x).all() and np.isfinite(y).all())}
    fa=np.vstack([np.loadtxt(p,skiprows=1,ndmin=2) for p in sorted(a.glob(f'face_dump_r*_call{call}.tsv'))])
    fb=np.vstack([np.loadtxt(p,skiprows=1,ndmin=2) for p in sorted(b.glob(f'face_dump_r*_call{call}.tsv'))])
    fa=fa[np.lexsort((fa[:,5],fa[:,4]))];fb=fb[np.lexsort((fb[:,5],fb[:,4]))]
    keys_same=fa.shape==fb.shape and np.array_equal(fa[:,4:6],fb[:,4:6])
    result['same_face_keys']=bool(keys_same)
    result['face_counts']=[len(fa),len(fb)]
    if keys_same:
        result['faces']={}
        for k,name in [(7,'pstar'),(8,'udotds'),(9,'work'),(10,'dSx'),(11,'dSy'),(12,'x1'),(13,'y1'),(14,'x2'),(15,'y2')]:
            x,y=fa[:,k],fb[:,k]
            result['faces'][name]={'max_abs':float(np.max(abs(x-y))), 'relative_l2':float(np.linalg.norm(x-y)/max(np.linalg.norm(x),1e-300))}
    for d in (a,b):
        log=(d/'log').read_text()
        rows=re.findall(r'\[GPU tess\].*?time=([0-9.]+) ms',log)
        vals=[float(v) for v in rows[8:]]
        result[d.name+'_tess_ms_median']=float(np.median(vals)) if vals else None
    return result

if __name__=='__main__':
    out=compare(sys.argv[1],sys.argv[2],int(sys.argv[3]))
    text=json.dumps(out,indent=2)+'\n'
    (ROOT/sys.argv[4]).write_text(text)
    print(text)

"""Compare explicit khout snapshots by persistent ID using verified C offsets."""
import argparse
import hashlib
import json
import struct
from pathlib import Path
import numpy as np


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--layout',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('reference')
    ap.add_argument('candidate')
    args=ap.parse_args()
    layout=json.loads(Path(args.layout).read_text())
    fields=[k for k in ('x','y','mass','vx','vy','den','pressure','ie','volume','w2','w2ceil','csound') if k in layout]
    def read(path):
        raw=Path(path).read_bytes()
        magic,n,size,stress=struct.unpack_from('<iiii',raw)
        assert magic==-1 and size==layout['size'] and stress==1
        assert len(raw)==16+n*size+52
        dtype=np.dtype({'names':['id']+fields,'formats':['<u8']+['<f8']*len(fields),
                        'offsets':[layout['u4if']]+[layout[k] for k in fields],'itemsize':size})
        a=np.frombuffer(raw,dtype=dtype,count=n,offset=16)
        ids=a['id'] & layout['id_mask']
        assert len(np.unique(ids))==n
        order=np.argsort(ids)
        return a[order],ids[order],struct.unpack_from('<dd',raw,16+n*size),hashlib.sha256(raw).hexdigest()
    a,ids,ta,ha=read(args.reference)
    b,ids_b,tb,hb=read(args.candidate)
    assert np.array_equal(ids,ids_b),'particle ID mismatch'
    assert ta==tb,'time or dt mismatch: not the same state time'
    reports={}
    for k in fields:
        x,y=a[k],b[k]
        assert np.isfinite(x).all() and np.isfinite(y).all(),('nonfinite',k)
        diff=np.abs(x-y)
        bitdiff=x.copy().view(np.uint64)!=y.copy().view(np.uint64)
        reports[k]={'bitwise_different_particles':int(np.count_nonzero(bitdiff)),
                    'max_abs_difference':float(diff.max()),
                    'worst_particle_id':int(ids[np.argmax(diff)])}
        if k == 'w2ceil':
            assert np.all(x>=0) and np.all(y>=0), 'negative squared distance'
            # Nonnegative finite IEEE doubles preserve integer bit ordering.
            # Canonicalize signed zero before interpreting the integer bits.
            xx=np.where(x==0,0.,x).copy().view(np.int64)
            yy=np.where(y==0,0.,y).copy().view(np.int64)
            reports[k]['max_ulp_difference']=int(np.max(np.abs(xx-yy)))
            nz=x!=0
            reports[k]['max_relative_difference_nonzero_reference']=float(np.max(diff[nz]/x[nz])) if np.any(nz) else 0.
            reports[k]['zero_reference_nonzero_candidate']=int(np.count_nonzero((x==0)&(y!=0)))
    result={'reference':args.reference,'candidate':args.candidate,'reference_sha256':ha,'candidate_sha256':hb,
            'n':len(ids),'t':ta[0],'dt':ta[1],'fields':reports,
            'all_physical_fields_bitwise_equal':all(v['bitwise_different_particles']==0 for v in reports.values()),
            'scope':'Listed physical fields at this snapshot only; pointer bytes and unlisted diagnostics excluded.'}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'n':len(ids),'t':ta[0],'all_physical_fields_bitwise_equal':result['all_physical_fields_bitwise_equal']}))

if __name__=='__main__':main()

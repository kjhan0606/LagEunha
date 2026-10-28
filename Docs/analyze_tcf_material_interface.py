"""Material-ID interface heat and reverse-face Q symmetry audit."""
import argparse
import json
import struct
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser()
ap.add_argument('root');ap.add_argument('--layout',required=True)
ap.add_argument('--initial',required=True);ap.add_argument('--output',required=True)
args=ap.parse_args();layout=json.loads(Path(args.layout).read_text())
raw=Path(args.initial).read_bytes();magic,n,size,stress=struct.unpack_from('<iiii',raw)
assert magic==-1 and stress==1 and size==layout['size'] and len(raw)==16+n*size+52
dtype=np.dtype({'names':['id','x','y','den'],'formats':['<u8','<f8','<f8','<f8'],
                'offsets':[layout[k] for k in ['u4if','x','y','den']], 'itemsize':size})
init=np.frombuffer(raw,dtype=dtype,count=n,offset=16)
ids=init['id']&layout['id_mask'];init=init[np.argsort(ids)]
assert np.array_equal(np.sort(ids),np.arange(n))
r0=np.hypot(init['x']-3,init['y']-3)
disk=init['den']>.5
group=np.full(n,3)
group[(~disk)&(r0<=.5)]=0
group[disk&(r0<.6)]=1
group[disk&(r0>=.6)]=2
labels=['inner_background','disk_inner_edge','remaining_disk','other']
out={}
for case in ['t4','t6']:
    f=np.loadtxt(Path(args.root)/case/'rhs.dat.faces')
    i=f[:,0].astype(int);j=f[:,1].astype(int)
    assert np.all((i>=0)&(i<n)) and np.all((j>=0)&(j<n))
    heat=.5*f[:,4]*f[:,14]*f[:,9]
    table={}
    for g in [0,1]:
        select=group[i]==g;total=heat[select].sum()
        table[labels[g]]={'total_heat':float(total),'neighbor_groups':{
            labels[h]:{'faces':int(np.count_nonzero(select&(group[j]==h))),
                       'heat':float(heat[select&(group[j]==h)].sum()),
                       'fraction':float(heat[select&(group[j]==h)].sum()/total)} for h in range(4)}}
    # Only the internal material interface; no periodic image ambiguity.
    cross=((group[i]==0)&(group[j]==1))|((group[i]==1)&(group[j]==0))
    rows=np.flatnonzero(cross);lookup={}
    for k in rows:
        key=(int(i[k]),int(j[k]));assert key not in lookup
        lookup[key]=k
    force=-f[:,4,None]*f[:,14,None]*f[:,6:8]
    ferr=[];herr=[];scale=[];pairs=0
    for k in rows:
        if i[k]>=j[k]:continue
        assert (int(j[k]),int(i[k])) in lookup
        q=lookup[(int(j[k]),int(i[k]))];pairs+=1
        ferr.append(np.linalg.norm(force[k]+force[q]));scale.append(np.linalg.norm(force[k])+np.linalg.norm(force[q]))
        herr.append(abs(heat[k]-heat[q]))
    assert 2*pairs==len(rows)
    out[case]={'heat_by_material':table,'interface_pairs':pairs,
               'Q_force_pair_residual_L1':float(sum(ferr)),
               'Q_force_pair_relative_L1':float(sum(ferr)/max(sum(scale),1e-300)),
               'half_heat_pair_residual_L1':float(sum(herr)),
               'half_heat_pair_relative_L1':float(sum(herr)/max(heat[cross].sum(),1e-300))}
result={'cases':out,'scope':'Material membership from actual initial density and radius, persistent IDs. Q force only, printed precision, two snapshots. Does not establish pressure-force symmetry, physical interface stability, or later topology correctness.'}
Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

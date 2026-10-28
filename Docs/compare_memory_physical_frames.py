"""Compare cached physical diagnostics only at matching saved times."""
import json
from pathlib import Path

base=Path('/gpfs/kjhan/LagForce')
roots=[base/'Codex_reversible_memory_physical_to_t10_20261018',
       base/'Codex_corotating_memory_physical_to_t10_20261018']
frames=[]
for root in roots:
    folder=root/'physical_frames/memory_from_t0'
    frames.append(sorted([json.loads(p.read_text()) for p in folder.glob('khout.*.dat.json')],key=lambda x:x['time']))
rows=[]
pairs=[];i=j=0
for series in frames:
    assert all(series[k+1]['time']-series[k]['time']>1e-12 for k in range(len(series)-1)), 'Ambiguous duplicate physical times'
while i<len(frames[0]) and j<len(frames[1]):
    a,b=frames[0][i],frames[1][j]
    if abs(a['time']-b['time'])<1e-12:
        pairs.append((a,b));i+=1;j+=1
    elif a['time']<b['time']:i+=1
    else:j+=1
for a,b in pairs:
    assert all(x['mass_unchanged'] and x['weights_zero'] for x in (a,b))
    regions={}
    for region,am in a['regions'].items():
        bm=b['regions'][region];assert am['cells']==bm['cells']
        regions[region]={key:dict(original=av,corotating=bm[key],
            relative_change_percent=100*(bm[key]/av-1) if av!=0 else None)
            for key,av in am.items() if key!='cells'}
    rows.append(dict(time=a['time'],sources=[dict(snapshot=x['snapshot'],sha256=x['sha256']) for x in (a,b)],regions=regions))
assert rows
print(json.dumps(dict(scope=__doc__,matching_frames=len(rows),latest=rows[-1],history=rows,
    matching_policy='Physical time within 1e-12, independent of file names; no interpolation.',
    available_frames=[len(f) for f in frames],unmatched_frames=[len(f)-len(rows) for f in frames],
    limits='Auxiliary-memory corotation comparison; not independent cell spin. Mixed diagnostic changes do not prove long-term stability. Binary and parallel-row settings differ; prior component controls are separate evidence.'),indent=2))

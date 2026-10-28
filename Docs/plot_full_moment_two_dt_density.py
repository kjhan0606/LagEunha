"""Density fields of two accepted checkpoints, not particle scatter.

Arc sampling is for display only. Show cell averages and errors relative to
the steady density integrated on each actual moving cell; no data correction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm,Normalize
from compare_full_moment_checkpoint_errors import load


def polygons(geometry,ids):
    sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_full_disk_arc_roots_modules_v7_20261019')
    from curved_disk_cell import DiskCell
    with np.load(geometry/'mesh.npz') as mesh:
        assert np.array_equal(ids,mesh['ids'])
        centers=mesh['moments'][:,1:3]/mesh['moments'][:,0,None]
    native=json.loads((geometry/'full_exact_supports.json').read_text())
    result=[]
    assert len(native['cells'])==len(ids)
    for i,record in enumerate(native['cells']):
        assert record['owner']==int(ids[i])
        cell=DiskCell.from_native(record,2.5);points=[]
        for line in cell.lines:points.extend([list(map(float,line.A)),list(map(float,line.B))])
        for arc in cell.arcs:
            theta=np.linspace(arc.lo,arc.hi,max(4,int(np.ceil((arc.hi-arc.lo)*400)))+1)
            points.extend((arc.radius*np.column_stack([np.cos(theta),np.sin(theta)])).tolist())
        xy=np.unique(np.asarray(points),axis=0);delta=xy-centers[i]
        result.append(xy[np.argsort(np.arctan2(delta[:,1],delta[:,0]))])
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('a','b','geometry-a','geometry-b','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--label-a',default='coarse h=0.005')
    p.add_argument('--label-b',default='fine h=0.0025')
    args=p.parse_args()
    a,ra,pa=load(args.a);b,rb,pb=load(args.b)
    assert ra['time']==rb['time'] and pa['physics_sha256']==pb['physics_sha256']
    assert np.array_equal(a['ids'],b['ids'])
    geometry_a=args.geometry_a.resolve();geometry_b=args.geometry_b.resolve()
    # The physics report was generated from exactly these moving geometries.
    for geometry,source in [(geometry_a,args.a),(geometry_b,args.b)]:
        reference=json.loads((source.parent/'reference/report.json').read_text())
        assert Path(reference['actual_state_path']).resolve().parent==geometry.parent
    sets=[polygons(geometry_a,a['ids']),polygons(geometry_b,b['ids'])]
    rho_limits=(min(.01,float(a['actual_density'].min()),float(b['actual_density'].min())),
                max(1.,float(a['actual_density'].max()),float(b['actual_density'].max())))
    error_max=max(float(abs(a['relative_density_reference_error']).max()),float(abs(b['relative_density_reference_error']).max()),1e-15)
    fig,axes=plt.subplots(2,3,figsize=(14,9),constrained_layout=True)
    for row,(data,poly,label) in enumerate(zip((a,b),sets,(args.label_a,args.label_b))):
        for col in range(3):
            ax=axes[row,col]
            values=data['actual_density'] if col==0 else data['relative_density_reference_error']
            norm=LogNorm(*rho_limits) if col==0 else Normalize(-error_max,error_max)
            collection=PolyCollection(poly,array=values,cmap='viridis' if col==0 else 'RdBu_r',norm=norm,edgecolors='none',rasterized=True)
            ax.add_collection(collection)
            extent=.45 if col==2 else 2.55
            ax.set(xlim=(-extent,extent),ylim=(-extent,extent),xlabel='x',ylabel='y')
            ax.set_aspect('equal')
            title='Cell mean density m/A' if col==0 else ('Inner density relative error' if col==2 else 'Density relative error')
            ax.set_title(label+'\n'+title)
            fig.colorbar(collection,ax=ax,shrink=.8)
    fig.suptitle(f"Kepler density, t={ra['time']:g}; common colour scales\nSteady reference on each actual cell; this short checkpoint does not establish long-time stability")
    args.output.mkdir(exist_ok=False)
    image=args.output/'two_dt_density.png';fig.savefig(image,dpi=180);plt.close(fig)
    sha=lambda file:hashlib.sha256(Path(file).read_bytes()).hexdigest()
    manifest=dict(time=ra['time'],cells=len(a['ids']),image=str(image),source_A=pa,source_B=pb,
        source_geometry_sha256={str(file):sha(file) for root in (geometry_a,geometry_b) for file in (root/'mesh.npz',root/'full_exact_supports.json')},
        rendering='Cell averages, no particle scatter; true arcs sampled for display only',
        common_density_limits=rho_limits,common_relative_density_error_limits=[-error_max,error_max],
        physical_providers=0,accepted_time_steps=0,long_time_stability_verified=False)
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(image=str(image),time=ra['time'],cells=len(a['ids']))),flush=True)


if __name__=='__main__':main()

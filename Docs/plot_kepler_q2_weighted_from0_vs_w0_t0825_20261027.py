"""Render same-time native W0 and fresh weighted density fields at t=0.825."""

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np

from plot_kepler_q2_weighted_vs_w0_density_20261027 import case_fields
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS=Path(__file__).resolve().parent
W0=DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T0825_2026-10-27.json'
WEIGHTED=DOCS/'KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/recovery_t0825_audit.json'
OUTPUT=DOCS/'KEPLER_Q2_FRESH_WEIGHTED_VS_W0_DENSITY_T0825_2026-10-27.png'


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    refs=[json.loads(path.read_text()) for path in (W0,WEIGHTED)]
    if refs[0]['time']!=.825 or refs[1]['time']!=.825 or not all(refs[1]['checks'].values()):
        raise ValueError('Wrong time or unaccepted weighted endpoint')
    starts=[Path(x['initial']) for x in refs]
    ends=[Path(x['endpoint']) for x in refs]
    if refs[0]['initial_sha256']!=refs[1]['initial_sha256']:
        raise ValueError('IC differs')
    for path,item in zip(starts,refs):
        if digest(path)!=item['initial_sha256'] or abs(snapshot_time(path,55351))>1e-12:
            raise ValueError('Changed t0 snapshot')
    for path,item in zip(ends,refs):
        if digest(path)!=item['endpoint_sha256'] or abs(snapshot_time(path,55351)-.825)>1e-10:
            raise ValueError('Changed endpoint')
    pixels=560
    domain=(1.,5.)
    axis=np.linspace(*domain,pixels)
    xx,yy=np.meshgrid(axis,axis)
    points=np.column_stack((xx.ravel(),yy.ravel()))
    theta=np.arange(1024)*2*np.pi/1024
    boundary=np.column_stack((3+2.75*np.cos(theta),3+2.75*np.sin(theta)))
    inside=MplPath(boundary).contains_points(points).reshape(xx.shape)
    fig,panels=plt.subplots(2,2,figsize=(10.6,9.2),constrained_layout=True)
    for j,(name,mode) in enumerate((('W=0 Voronoi',False),('Evolving Laguerre',True))):
        density,relative,disk_abs=case_fields(starts[j],ends[j],points,mode)
        density=np.ma.masked_where(~inside,density.reshape(xx.shape))
        relative=np.ma.masked_invalid(np.where(inside,relative.reshape(xx.shape),np.nan))
        top=panels[0,j].imshow(density,origin='lower',extent=(*domain,*domain),
            norm=LogNorm(vmin=.01,vmax=1.1,clip=True),cmap='viridis',interpolation='nearest')
        bottom=panels[1,j].imshow(relative,origin='lower',extent=(*domain,*domain),
            norm=SymLogNorm(linthresh=1e-5,linscale=1.,vmin=-.005,vmax=.005,base=10),
            cmap='coolwarm',interpolation='nearest')
        panels[0,j].set_title(name)
        panels[1,j].set_title(f"disk mean |Σ/Σ₀−1| = {disk_abs:.3g}")
    fig.colorbar(top,ax=panels[0,:].tolist(),label='surface density Σ')
    fig.colorbar(bottom,ax=panels[1,:].tolist(),label='Σ/Σ₀ − 1')
    for ax in panels.flat:
        ax.set_aspect('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xticks((1,2,3,4,5))
        ax.set_yticks((1,2,3,4,5))
    fig.suptitle('Native Kepler density at t=0.825, identical IC and color scales')
    fig.savefig(OUTPUT,dpi=150)
    plt.close(fig)
    print(OUTPUT)


if __name__=='__main__':
    main()

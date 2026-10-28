"""Show accepted aligned, ring-coherent and site-jitter native density fields."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np
from scipy.spatial import cKDTree

from plot_kepler_bounded_aligned_jitter_density_20261027 import accepted_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('aligned_audit', type=Path)
    ap.add_argument('ring_audit', type=Path)
    ap.add_argument('jitter_audit', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cases = [accepted_case('Aligned rings', args.aligned_audit),
             accepted_case('5% ring-coherent', args.ring_audit),
             accepted_case('5% independent site jitter', args.jitter_audit)]
    if any(abs(c['time']-6.5)>1e-10 for c in cases):
        raise ValueError('Expected three accepted native t=6.5 endpoints')
    domain=(1.0,5.0)
    pixels=600
    axis=np.linspace(*domain,pixels)
    xx,yy=np.meshgrid(axis,axis)
    query=np.column_stack((xx.ravel(),yy.ravel()))
    theta=np.arange(1024)*2*np.pi/1024
    boundary=np.column_stack((3+2.75*np.cos(theta),
                              3+2.75*np.sin(theta)))
    inside=MplPath(boundary).contains_points(query).reshape(pixels,pixels)
    fig,panels=plt.subplots(2,3,figsize=(15.2,9.3),constrained_layout=True)
    for col,case in enumerate(cases):
        nearest=cKDTree(case['xy']).query(query,k=1,workers=1)[1]
        density=np.ma.masked_where(~inside,case['density'][nearest].reshape(pixels,pixels))
        relative=np.ma.masked_invalid(np.where(
            inside,case['relative'][nearest].reshape(pixels,pixels),np.nan))
        top=panels[0,col].imshow(density,origin='lower',extent=(*domain,*domain),
            norm=LogNorm(vmin=.01,vmax=1.1,clip=True),cmap='viridis',interpolation='nearest')
        bottom=panels[1,col].imshow(relative,origin='lower',extent=(*domain,*domain),
            norm=SymLogNorm(linthresh=1e-5,linscale=1.,vmin=-.005,vmax=.005,base=10),
            cmap='coolwarm',interpolation='nearest')
        panels[0,col].set_title(case['label'])
        panels[1,col].set_title(f"disk mean |Σ/Σ₀−1| = {case['disk_abs']:.3g}")
    fig.colorbar(top,ax=panels[0,:].tolist(),label='surface density Σ')
    fig.colorbar(bottom,ax=panels[1,:].tolist(),label='Σ/Σ₀ − 1')
    for ax in panels.flat:
        ax.set_aspect('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xticks((1,2,3,4,5))
        ax.set_yticks((1,2,3,4,5))
    fig.suptitle('Native bounded Kepler disk at t=6.5; common density scales and nearest-site display')
    fig.savefig(args.output,dpi=150)
    plt.close(fig)
    print(args.output)


if __name__=='__main__':
    main()

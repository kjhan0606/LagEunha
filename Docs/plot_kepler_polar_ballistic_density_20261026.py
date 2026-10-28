"""Render comparable density fields from exact-orbit Voronoi controls."""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.spatial import cKDTree

from kepler_hybrid_polar_ballistic_20261026 import read_glass, polygon_areas, EPS, N, H

T = 2.221441469079183
DOMAIN = (.75, 5.25)
PIXELS = 720


def maps(xy0):
    dx0,dy0=xy0[:,0]-3.,xy0[:,1]-3.
    r0=np.hypot(dx0,dy0)
    selected=np.flatnonzero(r0<2.35)
    area0=polygon_areas(xy0,selected)
    angle=T*(r0*r0+EPS*EPS)**(-.75)
    xy=np.column_stack((3.+dx0*np.cos(angle)-dy0*np.sin(angle),
                        3.+dx0*np.sin(angle)+dy0*np.cos(angle)))
    area=polygon_areas(xy,selected)
    shape=np.where(r0[selected]<.5,(r0[selected]/.5)**3,
                   np.where(r0[selected]<=2.,1.,(1.+(r0[selected]-2.)/.1)**-3))
    rho0=.01+shape
    density=np.full(N*N,np.nan)
    density[selected]=rho0*area0/area
    error=np.full(N*N,np.nan)
    disk=(r0[selected]>=.5)&(r0[selected]<2.)
    error[selected[disk]]=abs(area0[disk]/area[disk]-1.)
    p=np.linspace(DOMAIN[0],DOMAIN[1],PIXELS)
    gx,gy=np.meshgrid(p,p)
    nearest=cKDTree(xy).query(np.column_stack((gx.ravel(),gy.ravel())),k=1)[1]
    return density[nearest].reshape(PIXELS,PIXELS),error[nearest].reshape(PIXELS,PIXELS)


def main():
    ids=np.arange(N*N)
    cartesian=np.column_stack(((ids%N+.5)*H,(ids//N+.5)*H))
    layouts=[('Cartesian',cartesian),
             ('Polar q=2, no center',read_glass(Path('/tmp/kepler_hybrid_polar_q2c_384.bin')))]
    fields=[maps(xy) for _,xy in layouts]
    fig,ax=plt.subplots(2,2,figsize=(10,9),constrained_layout=True)
    for j,((label,_),(density,error)) in enumerate(zip(layouts,fields)):
        top=ax[0,j].imshow(np.ma.masked_invalid(density),origin='lower',extent=(*DOMAIN,*DOMAIN),
                           norm=LogNorm(vmin=.01,vmax=1.1,clip=True),cmap='viridis',interpolation='nearest')
        bottom=ax[1,j].imshow(np.ma.masked_invalid(np.maximum(error,1e-6)),origin='lower',extent=(*DOMAIN,*DOMAIN),
                              norm=LogNorm(vmin=1e-4,vmax=1,clip=True),cmap='magma',interpolation='nearest')
        ax[0,j].set_title(label)
        for i in range(2):
            ax[i,j].set_aspect('equal')
            ax[i,j].set_xlabel('x')
            ax[i,j].set_ylabel('y')
            ax[i,j].set_xticks([1,2,3,4,5])
            ax[i,j].set_yticks([1,2,3,4,5])
    fig.colorbar(top,ax=ax[0,:],label='surface density')
    fig.colorbar(bottom,ax=ax[1,:],label='|density / initial density - 1|')
    fig.suptitle('Kepler t=2.221: exact generator orbits and Voronoi areas only (no hydro)')
    output=Path(__file__).with_name('KEPLER_POLAR_DENSITY_GEOMETRY_T2221_2026-10-26.png')
    fig.savefig(output,dpi=180)
    print(output)


if __name__=='__main__':
    main()

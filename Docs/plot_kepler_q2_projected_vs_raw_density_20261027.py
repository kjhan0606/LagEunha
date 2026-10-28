"""Show same-t0 native projected/raw q2 density using true power-cell pixels."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np

from plot_kepler_q2_weighted_vs_w0_density_20261027 import accepted, case_fields
from prepare_kepler_sparse_area_restart_20261027 import digest


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("projected_receipt",type=Path)
    ap.add_argument("raw_receipt",type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--pixels",type=int,default=768)
    ap.add_argument("--relative-max",type=float,default=.005)
    ap.add_argument("--domain",type=float,nargs=2,default=(2.4,3.6))
    args=ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if (args.pixels<128 or args.pixels>2048 or args.relative_max<=0
            or not np.all(np.isfinite(args.domain))
            or not 1.<=args.domain[0]<args.domain[1]<=5.):
        raise ValueError("Invalid density raster parameters")
    a,start_a,end_a,ta,failed_a=accepted(args.projected_receipt,())
    b,start_b,end_b,tb,failed_b=accepted(args.raw_receipt,())
    if (failed_a or failed_b or a["initial_sha256"]!=b["initial_sha256"]
            or digest(start_a)!=digest(start_b)
            or abs(ta-6.5)>1e-10 or abs(tb-6.5)>1e-10):
        raise ValueError("Expected fully accepted same-t0 exact t=6.5 native pair")
    axis=np.linspace(*args.domain,args.pixels)
    xx,yy=np.meshgrid(axis,axis)
    points=np.column_stack((xx.ravel(),yy.ravel()))
    theta=np.arange(1024)*2*np.pi/1024
    polygon=np.column_stack((3+2.75*np.cos(theta),3+2.75*np.sin(theta)))
    inside=MplPath(polygon).contains_points(points).reshape(xx.shape)
    fig,panels=plt.subplots(2,2,figsize=(10.5,9),constrained_layout=True)
    for col,(label,start,end) in enumerate((
            ("Pressure-work projected",start_a,end_a),
            ("Raw Jacobi weight rate",start_b,end_b))):
        rho,residual,disk_abs=case_fields(start,end,points,True)
        rho=np.ma.masked_where(~inside,rho.reshape(xx.shape))
        residual=np.ma.masked_invalid(np.where(inside,residual.reshape(xx.shape),np.nan))
        top=panels[0,col].imshow(
            rho,origin="lower",extent=(*args.domain,*args.domain),
            norm=LogNorm(vmin=.01,vmax=1.1,clip=True),cmap="viridis",interpolation="nearest")
        bottom=panels[1,col].imshow(
            residual,origin="lower",extent=(*args.domain,*args.domain),
            norm=SymLogNorm(linthresh=args.relative_max/500,
                            vmin=-args.relative_max,vmax=args.relative_max,base=10),
            cmap="coolwarm",interpolation="nearest")
        panels[0,col].set_title(f"{label} (21/21 gates)")
        panels[1,col].set_title(f"disk mean |ΔΣ/Σ₀| = {disk_abs:.4g}")
    fig.colorbar(top,ax=panels[0,:].tolist(),label="surface density Σ")
    fig.colorbar(bottom,ax=panels[1,:].tolist(),label="Σ/Σ₀ − 1")
    for panel in panels.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        ticks=np.linspace(*args.domain,5)
        panel.set_xticks(ticks)
        panel.set_yticks(ticks)
    fig.suptitle("Native bounded q2 Kepler at t=6.5; same t=0 IC; common density scales")
    fig.savefig(args.output,dpi=150)
    plt.close(fig)
    print(args.output)


if __name__=="__main__":
    main()

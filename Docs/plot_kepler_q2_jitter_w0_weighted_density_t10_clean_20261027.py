"""Compact same-time cell-density comparison for jittered bounded q2 Kepler."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np

from plot_kepler_q2_weighted_vs_w0_density_20261027 import DOMAIN, case_fields
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
OUTPUT = DOCS / "KEPLER_Q2_JITTER_W0_WEIGHTED_DENSITY_T10_CLEAR_2026-10-27.png"
RECEIPTS = (
    ("Voronoi W=0", ROOT / "jitter005_to10_a10_v1/audit_t10.json", False),
    ("Laguerre weights", ROOT / "gpu_energy_neutral_weighted_t10_a10_v1/jitter/audit_t10.json", True),
)


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    axis = np.linspace(*DOMAIN, 512)
    xx,yy = np.meshgrid(axis,axis)
    points = np.column_stack((xx.ravel(),yy.ravel()))
    theta=np.arange(1024)*2*np.pi/1024
    polygon=np.column_stack((3+2.75*np.cos(theta),3+2.75*np.sin(theta)))
    inside=MplPath(polygon).contains_points(points).reshape(xx.shape)
    fig,panels=plt.subplots(2,2,figsize=(12.5,8.5),constrained_layout=True)
    reference=None
    for col,(title,receipt_path,weighted) in enumerate(RECEIPTS):
        receipt=json.loads(receipt_path.read_text())
        start=Path(receipt["initial"])
        end=Path(receipt["endpoint"])
        if (digest(start)!=receipt["initial_sha256"]
                or digest(end)!=receipt["endpoint_sha256"]
                or abs(snapshot_time(end,55351)-10.)>1e-10):
            raise ValueError(f"{title}: changed or wrong-time native state")
        if reference is None:
            reference=receipt["initial_sha256"]
        elif reference!=receipt["initial_sha256"]:
            raise ValueError("Different initial fields")
        density,residual,mean_abs=case_fields(start,end,points,weighted)
        if abs(mean_abs-receipt["disk_density_abs_mass_weighted"])>1e-10:
            raise ValueError(f"{title}: density metric mismatch")
        density=np.ma.masked_where(~inside,density.reshape(xx.shape))
        residual=np.ma.masked_invalid(np.where(inside,residual.reshape(xx.shape),np.nan))
        im0=panels[0,col].imshow(density,origin="lower",extent=(*DOMAIN,*DOMAIN),
                 norm=LogNorm(vmin=.01,vmax=1.1,clip=True),cmap="viridis",interpolation="nearest")
        im1=panels[1,col].imshow(residual,origin="lower",extent=(*DOMAIN,*DOMAIN),
                 norm=SymLogNorm(linthresh=1e-5,vmin=-.01,vmax=.01,base=10),
                 cmap="coolwarm",interpolation="nearest")
        panels[0,col].set_title(title)
        panels[1,col].set_title(f"mean |ΔΣ/Σ₀| = {mean_abs:.4g}")
    fig.colorbar(im0,ax=panels[0,:].tolist(),label="surface density Σ")
    fig.colorbar(im1,ax=panels[1,:].tolist(),label="Σ/Σ₀ − 1")
    for panel in panels.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        panel.set_xticks((1,2,3,4,5))
        panel.set_yticks((1,2,3,4,5))
    fig.suptitle("Native bounded q2 Kepler, 5% site jitter, t=10")
    fig.savefig(OUTPUT,dpi=150)
    plt.close(fig)
    print(OUTPUT)


if __name__=="__main__":
    main()

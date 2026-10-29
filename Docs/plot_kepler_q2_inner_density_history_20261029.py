"""Render saved Q2 inner-disk density at a shared contrast scale."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import mean

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
history=json.loads((docs/'KEPLER_Q2_PRESSURE_HISTORY_2026-10-29.json').read_text())
initial=Path(prior['initial']);middle=Path(prior['cases'][0]['endpoint'])
late=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1/khout.011058.dat')
png=docs/'KEPLER_Q2_INNER_DENSITY_T0_T30_T47P5_2026-10-29.png'
receipt=docs/'KEPLER_Q2_INNER_DENSITY_T0_T30_T47P5_2026-10-29.json'
if png.exists() or receipt.exists():raise FileExistsError('Refusing to overwrite')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial changed')
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
if ma['t']!=0 or len(a)!=55351:raise ValueError('Unexpected initial state')
r0=np.hypot(a['x']-3.,a['y']-3.);disk=(r0>=.5)&(r0<2.)
rho_ref=mean(a['den'][disk].astype(float),a['mass'][disk].astype(float))
fig,axes=plt.subplots(1,3,figsize=(14,4.8),constrained_layout=True)
extent=(2.1,3.9,2.1,3.9);pixels=900
xs=np.linspace(extent[0],extent[1],pixels);ys=np.linspace(extent[2],extent[3],pixels)
xx,yy=np.meshgrid(xs,ys);points=np.column_stack((xx.ravel(),yy.ravel()))
norm=TwoSlopeNorm(vmin=-5.,vcenter=0.,vmax=5.)
cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#20212b')
rows=[]
for ax,path in zip(axes,(initial,middle,late)):
    p,mt=ordered_global(path,ids);digest=sha(path)
    h=next((v for v in history['cases'] if v['time']==mt['t']),None)
    if path==initial:
        if mt['t']!=0:raise ValueError('Initial time changed')
    elif h is None or h['endpoint_sha256']!=digest:raise ValueError('Checkpoint changed')
    if not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Changed mass or weights')
    rho=p['den'].astype(float)
    if np.any(~np.isfinite(rho)) or np.any(rho<=0):raise ValueError('Invalid density')
    xy=np.column_stack((p['x'],p['y']));nearest=cKDTree(xy).query(points,k=1,workers=1)[1]
    contrast=100.*(rho[nearest]/rho_ref-1.);visible=disk[nearest]
    field=np.ma.masked_where(~visible,contrast).reshape(pixels,pixels)
    im=ax.imshow(field,origin='lower',extent=extent,interpolation='nearest',cmap=cmap,norm=norm)
    for radius in (.5,.75):ax.add_patch(Circle((3,3),radius,fill=False,color='black',alpha=.55,linestyle='--',linewidth=.7))
    ax.set(title=f't = {mt["t"]:.3f}',xlabel='x',ylabel='y',aspect='equal')
    rows.append({'time':mt['t'],'source':str(path),'sha256':digest,'visible_pixels':int(visible.sum()),
                 'saturated_visible_pixel_fraction':float(np.mean(np.abs(contrast[visible])>5.))})
fig.colorbar(im,ax=axes,orientation='horizontal',pad=.12,shrink=.83,extend='both',
             label=r'Surface density contrast $100(\Sigma/\Sigma_{ref}-1)$ [%], shared scale')
fig.suptitle('Q2 Kepler inner disk: saved states, identical color scale\n'
             'W=0, fixed masses; nearest-generator raster, ±5% saturation',fontsize=12)
fig.savefig(png,dpi=170);plt.close(fig)
result={'rho_reference':rho_ref,'contrast_limit_percent':5.,'extent':extent,'pixels_per_panel':pixels,
        'cases':rows,'initial_sha256':prior['initial_sha256'],'png_sha256':sha(png),
        'script_sha256':sha(Path(__file__)),
        'scope':'Fixed shared color scale, nearest-generator raster, masked by original disk cohort. t47.5 is before wall impact. Pixel saturation is a rendering statistic, not an error norm or a GIZMO comparison.'}
with receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))

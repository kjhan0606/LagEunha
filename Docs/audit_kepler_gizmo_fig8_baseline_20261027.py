"""Compute Hopkins Fig. 8 annular L1 velocity/density metrics on native states."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS=Path(__file__).resolve().parent
ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
SOURCES={
    'w0_jitter_t0825':DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T0825_2026-10-27.json',
    'weighted_jitter_t0825':DOCS/'KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/recovery_t0825_audit.json',
    'w0_aligned_t20':ROOT/'fast_to20_a10_v1/audit_t20.json',
    'w0_aligned_t30':ROOT/'wallselect_to30_a10_v1/audit_t30.json',
}
OUTPUT=DOCS/'KEPLER_GIZMO_FIG8_METRIC_BASELINE_2026-10-27.json'
EPS=.01175


def annulus_metrics(p,p0,mask):
    x=p['x'].astype(float)-3.
    y=p['y'].astype(float)-3.
    r=np.hypot(x,y)
    vx=p['vx'].astype(float)
    vy=p['vy'].astype(float)
    vphi=(-y*vx+x*vy)/np.maximum(r,1e-30)
    vk=r*(r*r+EPS*EPS)**(-.75)
    rho=p['den'].astype(float)
    rho0=p0['den'].astype(float)
    return {
        'sites':int(np.count_nonzero(mask)),
        'L1_vphi_fractional_paper':float(np.mean(np.abs(vphi[mask]-vk[mask])/vk[mask])),
        'L1_density_literal_from_one_paper':float(np.mean(np.abs(rho[mask]-1.))),
        'L1_density_relative_to_own_t0':float(np.mean(np.abs(rho[mask]/rho0[mask]-1.))),
        'mean_density':float(np.mean(rho[mask])),
        'mean_initial_density':float(np.mean(rho0[mask])),
    }


def main():
    result={}
    for name,source in SOURCES.items():
        audit=json.loads(source.read_text())
        if 'checks' in audit and not all(audit['checks'].values()):
            raise ValueError(f'{name} did not pass its native audit')
        initial=Path(audit['initial'])
        endpoint=Path(audit['endpoint'])
        if digest(initial)!=audit['initial_sha256'] or digest(endpoint)!=audit['endpoint_sha256']:
            raise ValueError(f'{name} snapshot changed')
        p0,m0=ordered(initial,55351)
        p,m=ordered(endpoint,55351)
        t=float(audit['time'])
        if abs(float(m0['t']))>1e-12 or abs(float(m['t'])-t)>1e-10:
            raise ValueError(f'{name} snapshot time mismatch')
        r=np.hypot(p['x'].astype(float)-3.,p['y'].astype(float)-3.)
        r0=np.hypot(p0['x'].astype(float)-3.,p0['y'].astype(float)-3.)
        current=(r>.8)&(r<1.2)
        original=(r0>.8)&(r0<1.2)
        result[name]={
            'receipt':str(source),'receipt_sha256':digest(source),
            'initial_sha256':audit['initial_sha256'],
            'endpoint_sha256':audit['endpoint_sha256'],'time':t,
            'orbits_at_r1_unsoftened':t/(2*np.pi),
            'current_radius_annulus':annulus_metrics(p,p0,current),
            'initial_radius_annulus':annulus_metrics(p,p0,original),
        }
    output={
        'paper':'Hopkins 2015, arXiv:1409.7395, section 4.2.4, Figures 7-8',
        'paper_url':'https://arxiv.org/pdf/1409.7395',
        'paper_annulus':'0.8 < R < 1.2 at current particle positions',
        'paper_velocity_metric':'arithmetic mean |v_phi-v_K(R)|/v_K(R)',
        'paper_density_metric':'arithmetic mean |Sigma-1|',
        'paper_fixed_time':float(20*np.pi),
        'paper_thresholds':{'L1_velocity':.01,'L1_density':.25},
        'paper_longest_figure_time':600.,
        'cases':result,
        'scope':'Same algebraic annular metrics on available native endpoints. Current bounded polar IC differs from the paper uniform effective-256-squared IC, so these values are not a direct reproduction. The own-t0 density metric removes the paper literal 0.01 floor offset; both are reported.'}
    with OUTPUT.open('x') as stream:
        json.dump(output,stream,indent=2,allow_nan=False)
        stream.write('\n')
    for key,row in result.items():
        print(key,row['time'],row['current_radius_annulus'])


if __name__=='__main__':
    main()

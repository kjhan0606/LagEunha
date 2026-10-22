"""Mass/volume weighted comparisons on retained initial IDs and native RHS."""
from pathlib import Path
import argparse,json,re,sys
import numpy as np
root=Path(__file__).resolve().parent;sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
args=argparse.ArgumentParser();args.add_argument('--run-dir',default='native_v3');args=args.parse_args()
base=root/'core_refinement';native=base/args.run_dir;output=base/('results_'+args.run_dir+'.json')
out={'scope':'fixed-eps discretization diagnosis; not a substitute for long Kepler validation','cases':{}}
manifest=native/'manifest.json'
if not manifest.exists():
    out['status']='native preflight did not finish; do not interpret missing cases as passes'
    output.write_text(json.dumps(out,indent=2)+'\n')
    raise SystemExit(0)
record=json.loads(manifest.read_text());out['preflight']=record.get('preflight')
gamma=float(np.float32(1.666667));eps=.047
def arrange(p):return p[np.argsort(p['u4if']&np.uint64((1<<56)-1))]
for case in record['cases']:
    d=native/case['name'];p0,m0=read_snapshot(d/'khout.000001.dat');p0=arrange(p0)
    r0=np.hypot(p0['x']-12,p0['y']-12)
    row={'exit_code':case['exit_code'],'classification':case['classification'],'initial_data':case['initial_data'],'radial':[]}
    snapshots=[p for p in d.glob('khout.*.dat') if p.name!='khout.000001.dat']
    if snapshots:
        path=max(snapshots,key=lambda p:int(p.name.split('.')[1]));p,m=read_snapshot(path);p=arrange(p)
        assert np.array_equal(p0['u4if'],p['u4if'])
        omega=(r0*r0+eps*eps)**(-.75);phase=omega*m['t'];c,s=np.cos(phase),np.sin(phase)
        x=p0['x']-12;y=p0['y']-12;xe=x*c-y*s;ye=x*s+y*c
        poserr=np.hypot(p['x']-12-xe,p['y']-12-ye)/r0
        verr=np.hypot(p['vx']+omega*ye,p['vy']-omega*xe)/(omega*r0)
        # Individual masses stay fixed, so entropy ratio is U/U0*(V/V0)^(gamma-1).
        entropy=p['ie'].astype(float)/p0['ie']*(p['volume'].astype(float)/p0['volume'])**(gamma-1)
        row.update(time=m['t'],checkpoint=str(path),mass_max_relative_change=float(np.max(abs(p['mass']-p0['mass'])/p0['mass'])),
                   total_internal_energy=float(p['ie'].astype(float).sum()),minimum_internal_energy=float(p['ie'].min()))
        radius=np.hypot(p['x']-12,p['y']-12)
        E0=float(np.sum(p0['ie'].astype(float)+p0['mass']*(.5*(p0['vx']**2+p0['vy']**2)-1/np.sqrt(r0*r0+eps*eps))))
        E1=float(np.sum(p['ie'].astype(float)+p['mass']*(.5*(p['vx']**2+p['vy']**2)-1/np.sqrt(radius*radius+eps*eps))))
        L0=float(np.sum(p0['mass']*((p0['x']-12)*p0['vy']-(p0['y']-12)*p0['vx'])))
        L1=float(np.sum(p['mass']*((p['x']-12)*p['vy']-(p['y']-12)*p['vx'])))
        row.update(energy_relative_change_from_IC=(E1-E0)/abs(E0),angular_momentum_relative_change_from_IC=(L1-L0)/abs(L0))
        for lo,hi in ((0,.1),(.1,.2),(.2,.5),(.5,.65),(.65,1),(1,2),(2,4),(4,8)):
            sel=(r0>=lo)&(r0<hi)
            if not sel.any():continue
            mass=p0['mass'][sel];U0=p0['ie'][sel].astype(float);V=p['volume'][sel].astype(float)
            Uad=U0*(p0['volume'][sel].astype(float)/V)**(gamma-1)
            row['radial'].append({'initial_Rlo':lo,'initial_Rhi':hi,'n':int(sel.sum()),'mass':float(mass.sum()),
                'entropy_error_mass_weighted':float(np.sum(mass*abs(entropy[sel]-1))/mass.sum()),
                'entropy_ratio_max':float(entropy[sel].max()),
                'orbital_position_error_mass_weighted':float(np.sum(mass*poserr[sel])/mass.sum()),
                'orbital_velocity_error_mass_weighted':float(np.sum(mass*verr[sel])/mass.sum()),
                'orbital_position_error_max':float(poserr[sel].max()),
                'internal_energy':float(p['ie'][sel].astype(float).sum()),
                'thermal_energy_excess_above_adiabatic':float(np.sum(p['ie'][sel].astype(float)-Uad)),
                'pressure_max':float(p['pressure'][sel].max())})
    rows=np.genfromtxt(d/'rk4_stage_r0.tsv',names=True,dtype=None,encoding=None,delimiter='\t')
    assert np.all(rows['stage']=='k1_rhs')
    Udot=np.full(len(p0),np.nan);V=np.full(len(p0),np.nan)
    Udot[rows['id']]=rows['Udot'];V[rows['id']]=rows['V']
    if np.isnan(Udot).any():raise ValueError('missing first RHS cell')
    row['first_RHS']=[]
    for lo,hi in ((0,.1),(.1,.2),(.2,.5),(.5,.65),(.65,1),(1,2),(2,4),(4,8)):
        sel=(r0>=lo)&(r0<hi)
        if not sel.any():continue
        row['first_RHS'].append({'initial_Rlo':lo,'initial_Rhi':hi,'n':int(sel.sum()),
            'sum_Udot':float(Udot[sel].sum()),'positive_Udot_sum':float(np.maximum(Udot[sel],0).sum()),
            'max_Udot_over_U':float(np.max(Udot[sel]/p0['ie'][sel])),
            'native_volume_relative_error_max':float(np.max(abs(V[sel]-p0['volume'][sel])/p0['volume'][sel]))})
    log=(d/'log').read_text()
    errors=re.findall(r'dEtot/\|Etot0\|=([0-9.eE+-]+)',log)
    row['global_energy_error_max']=max([abs(float(x)) for x in errors],default=None)
    floor_values=[float(v) for v in re.findall(r'floor_cum=([0-9.eE+-]+)',log)]
    row['floor_present']=any(v!=0 for v in floor_values)
    failed=sorted(d.glob('khdiag.step*.failed.dat'))
    if failed:
        p,m=read_snapshot(failed[0]);row['invalid_stage_snapshot']=str(failed[0]);row['invalid_stage_min_U']=float(p['ie'].min())
    out['cases'][case['name']]=row
    print(case['name'],'exit',case['exit_code'],'time',row.get('time'),'first_core_heat',sum(x['positive_Udot_sum'] for x in row['first_RHS'] if x['initial_Rhi']<=.5),flush=True)
output.write_text(json.dumps(out,indent=2)+'\n')

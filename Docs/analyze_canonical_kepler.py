"""Read explicit production snapshots using a compiler-verified layout.
Uses the run's actual Pakmor/spline IC, not a synthetic particle lattice.
Reports raw field errors; no long-time pass/fail is inferred from conservation.
"""
import argparse
import hashlib
import json
import re
import struct
from pathlib import Path
import numpy as np


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def circular_speed(r, rc):
    u = r / rc
    enclosed = np.ones_like(r)
    lo = u <= .5
    mid = (u > .5) & (u < 1)
    q = u[lo]
    enclosed[lo] = q**3 * (32/3 - 192/5*q*q + 32*q**3)
    q = u[mid]
    enclosed[mid] = -1/15 + q**3*(64/3 - 48*q + 192/5*q*q - 32/3*q**3)
    return np.sqrt(np.maximum(enclosed, 0)/np.maximum(r, 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--layout', required=True)
    ap.add_argument('--initial', required=True)
    ap.add_argument('--environment', required=True)
    ap.add_argument('--params', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('snapshots', nargs='+')
    args = ap.parse_args()
    layout = json.loads(Path(args.layout).read_text())
    env = json.loads(Path(args.environment).read_text())
    assert env['EUNHA_KEPLER_PROFILE'] == 'pakmor'
    assert env['EUNHA_KEPLER_SOFT'] == 'spline'
    assert env['EUNHA_KEPLER_BG_ROTATE'] == '1'
    assert layout['postype_size'] == layout['id_size'] == 8
    for key in ('EUNHA_KEPLER_FORCE', 'EUNHA_HILL_OMEGA'):
        assert key not in env, ('unsupported reference', key)
    with open(env['LAGMFM_GLASS_FILE'], 'rb') as f:
        n = struct.unpack('<i', f.read(4))[0]
        lx, ly = struct.unpack('<dd', f.read(16))
        xy = np.frombuffer(f.read(), dtype='<f8')
    assert len(xy) == 2*n
    volume0 = np.fromfile(env['EUNHA_IC_CELLVOL'], dtype='<f8')
    assert len(volume0) == n and np.all(volume0 > 0)
    cx, cy = float(env['EUNHA_KEPLER_CX']), float(env['EUNHA_KEPLER_CY'])
    x0, y0 = xy[:n]-cx, xy[n:]-cy
    r0 = np.hypot(x0, y0)
    rin, rout = float(env['EUNHA_KEPLER_RIN']), float(env['EUNHA_KEPLER_ROUT'])
    rc = float(env['EUNHA_KEPLER_RC'])
    rho0 = np.where((r0 >= rin) & (r0 <= rout), 1., float(env['EUNHA_KEPLER_FLOOR']))
    mass0 = rho0*volume0
    gamma_text = re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)', Path(args.params).read_text())[1]
    gamma = float(np.float32(float(gamma_text)))  # GAS_GAMMA is stored as float.
    p0 = float(env.get('EUNHA_KEPLER_P', (gamma-1)*(5*gamma/2)*1e-5))
    vc0 = circular_speed(r0, rc)
    vx0, vy0 = -y0/np.maximum(r0,1e-30)*vc0, x0/np.maximum(r0,1e-30)*vc0
    names = ['u4if','x','y','mass','vx','vy','den','pressure','ie','volume']
    dtype = np.dtype({'names':names, 'formats':['<u8']+['<f8']*(len(names)-1),
                      'offsets':[layout[k] for k in names], 'itemsize':layout['size']})
    def load_snapshot(path):
        with open(path,'rb') as f:
            magic,count,size,stress=struct.unpack('<iiii',f.read(16))
            assert magic==-1 and count==n and size==layout['size'] and stress==1
            a=np.frombuffer(f.read(count*size),dtype=dtype)
            trailer=f.read()
        assert len(a)==n and len(trailer)==52, 'incomplete snapshot'
        t,dt=struct.unpack_from('<dd',trailer)
        nx,ny=struct.unpack_from('<ii',trailer,20)
        bx,by=struct.unpack_from('<dd',trailer,28)
        assert nx*ny==n and (bx,by)==(lx,ly)
        ids=(a['u4if'] & layout['id_mask']).astype(np.int64)
        assert np.array_equal(np.sort(ids),np.arange(n)), 'bad particle IDs'
        a=a[np.argsort(ids)]
        assert all(np.isfinite(a[k]).all() for k in names[1:])
        return a,t,dt
    initial,t0,_=load_snapshot(args.initial)
    assert t0==0
    assert np.array_equal(initial['x']-cx,x0) and np.array_equal(initial['y']-cy,y0)
    # Exact C initial classification matters at rings on r=0.5 and r=2.
    # Python hypot and optimized C sqrt can classify boundary points differently.
    rho0=np.where(initial['den']>.5,1.,float(env['EUNHA_KEPLER_FLOOR'])); mass0=initial['mass']; vx0=initial['vx']; vy0=initial['vy']
    p0=float(initial['pressure'][0]); assert np.all(initial['pressure']==p0)
    vc0=np.hypot(vx0,vy0)
    freeze_r=float(env['EUNHA_KEPLER_EVOLVE_R'])
    ambiguous_freeze=np.abs(r0-freeze_r)<=1e-12
    output = {'reference':'Actual IC files, Pakmor uniform pressure, co-rotation, spline force',
              'gamma':gamma,'pressure0':p0,'initial_snapshot':args.initial,'initial_sha256':digest(args.initial),'layout':layout,'ic_sha256':{
              env[k]:digest(env[k]) for k in ('LAGMFM_GLASS_FILE','EUNHA_IC_CELLVOL')},
              'snapshots':[], 'passed_long_term_validation':False}
    for path in args.snapshots:
        a,t,dt=load_snapshot(path)
        x,y = a['x']-cx,a['y']-cy
        r = np.hypot(x,y)
        vc = circular_speed(r,rc)
        vr = (x*a['vx']+y*a['vy'])/np.maximum(r,1e-30)
        vp = (x*a['vy']-y*a['vx'])/np.maximum(r,1e-30)
        theta_err = np.angle(np.exp(1j*(np.arctan2(y,x)-np.arctan2(y0,x0)-vc0/np.maximum(r0,1e-30)*t)))
        regions = {'initial_inner_background':(rho0<1)&(r0<=rin),
                   'initial_disk':rho0==1,
                   'initial_disk_inner_edge':(rho0==1)&(r0<rin+.1),
                   'initial_disk_core':(rho0==1)&(r0>=rin+.1)&(r0<=rout-.1),
                   'initial_disk_outer_edge':(rho0==1)&(r0>rout-.1),
                   'initial_outer_active':(rho0<1)&(r0>=rout)&(r0<freeze_r)&~ambiguous_freeze,
                   'freeze_boundary_ambiguous':ambiguous_freeze,
                   'initial_frozen':r0>freeze_r+1e-12}
        for low,high in ((0,.125),(.125,.25),(.25,.375),(.375,.5)):
            regions[f'inner_r0_{low:g}_{high:g}']=(rho0<1)&(r0>=low)&(r0<high)
        disk_partition=sum(regions[k].astype(int) for k in
                           ('initial_disk_inner_edge','initial_disk_core','initial_disk_outer_edge'))
        assert np.array_equal(disk_partition,regions['initial_disk'].astype(int))
        assert np.all(a['pressure']>0) and np.all(a['den']>0), 'nonpositive state: entropy diagnostic invalid'
        ratio=a['den']/initial['den']
        pad=initial['pressure']*ratio**gamma
        pressure_error=a['pressure']/p0-1
        adiabatic_error=pad/p0-1
        nonadiabatic_error=(a['pressure']-pad)/p0
        log_entropy=np.log(a['pressure']/initial['pressure'])-gamma*np.log(ratio)
        # Endpoint entropy proxy only: not an attribution to HG, TCF, floors,
        # or physical shocks. That requires integrating their separate budgets.
        result = {'path':str(path),'t':t,'dt':dt,'particle_count':n,
                  'mass_max_abs_change':float(np.max(np.abs(a['mass']-mass0))),
                  'mass_total':float(a['mass'].sum()),'pressure_min':float(a['pressure'].min()),
                  'density_min':float(a['den'].min()),'regions':{}}
        for name,sel in regions.items():
            if not sel.any():
                result['regions'][name]={'n':0}
                continue
            weights=mass0[sel]
            def rms(v): return float(np.sqrt(np.average(v[sel]**2,weights=weights)))
            pressure_bias=float(np.average(pressure_error[sel],weights=weights))
            pressure_variance=float(np.average((pressure_error[sel]-pressure_bias)**2,weights=weights))
            pressure_mse=float(np.average(pressure_error[sel]**2,weights=weights))
            assert np.isclose(pressure_mse,pressure_bias**2+pressure_variance,rtol=1e-12,atol=1e-28)
            result['regions'][name]={'n':int(sel.sum()),
                'initial_mass':float(weights.sum()),
                'rho_relative_rms_material':rms(a['den']/rho0-1),
                'rho_change_rms_from_actual_initial':rms(ratio-1),
                'rho_relative_rms_material_reference':'ideal material density 1 or floor; includes initial tessellation density error',
                'pressure_relative_rms':rms(a['pressure']/p0-1),
                'pressure_relative_mean_bias':pressure_bias,
                'pressure_relative_std_about_mean':float(np.sqrt(pressure_variance)),
                'pressure_mse_mean_bias_fraction':pressure_bias**2/pressure_mse if pressure_mse>0 else None,
                'pressure_error_mass_weighted_sum_squares':float(np.sum(weights*pressure_error[sel]**2)),
                'radial_displacement_rms':rms(r-r0),
                'vr_over_vc_rms':rms(vr/np.maximum(vc,1e-30)),
                'vphi_relative_rms':rms(vp/np.maximum(vc,1e-30)-1),
                'Lz_change':float(np.sum(a['mass'][sel]*(x[sel]*a['vy'][sel]-y[sel]*a['vx'][sel])-mass0[sel]*(x0[sel]*vy0[sel]-y0[sel]*vx0[sel])))}
            result['regions'][name]['thermal_decomposition']={
                'adiabatic_pressure_error_rms':rms(adiabatic_error),
                'nonadiabatic_pressure_error_rms':rms(nonadiabatic_error),
                'cross_term_2mean_ad_nonad':float(2*np.average((adiabatic_error*nonadiabatic_error)[sel],weights=weights)),
                'log_entropy_rms':rms(log_entropy),
                'log_entropy_mean':float(np.average(log_entropy[sel],weights=weights)),
                'log_entropy_min':float(log_entropy[sel].min()),
                'log_entropy_max':float(log_entropy[sel].max()),
                'negative_entropy_particle_fraction':float(np.mean(log_entropy[sel]<-1e-10)),
                'pressure_min_over_initial':float(np.min(a['pressure'][sel]/p0)),
                'pressure_max_over_initial':float(np.max(a['pressure'][sel]/p0))}
            if name not in ('initial_frozen','freeze_boundary_ambiguous'): result['regions'][name]['orbital_phase_error_rms_rad']=rms(theta_err)
            else:
                result['regions'][name]['velocity_change_max']=float(np.max(np.hypot(a['vx'][sel]-vx0[sel],a['vy'][sel]-vy0[sel])))
                result['regions'][name]['position_change_max']=float(np.max(np.hypot(x[sel]-x0[sel],y[sel]-y0[sel])))
        inner_ids=np.flatnonzero(regions['initial_inner_background'])
        worst=inner_ids[np.argsort(np.abs(log_entropy[inner_ids]))[-10:][::-1]]
        result['inner_top_entropy_errors']=[{'id':int(i),'r0':float(r0[i]),'r':float(r[i]),
            'log_entropy':float(log_entropy[i]),'rho_over_initial':float(ratio[i]),
            'pressure_over_initial':float(a['pressure'][i]/p0)} for i in worst]
        result['thermal_scope']='Endpoint adiabatic reconstruction and entropy proxy; no attribution to a particular source term.'
        output['snapshots'].append(result)
    Path(args.output).write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'output':args.output,'times':[x['t'] for x in output['snapshots']],
                      'mass_errors':[x['mass_max_abs_change'] for x in output['snapshots']]}))

if __name__=='__main__': main()

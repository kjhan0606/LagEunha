"""Read native canonical CUDA flux packets; no time integration."""
from pathlib import Path
import argparse
import json

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DTYPE=np.dtype([('i','<i8'),('j','<i8'),('a','<f8',(29,)),('flags','<i4',(4,))])
FOCUS=(67935,68001,68526,68592)


def load(d):
    data=[];cells=[];metas=[]
    for path in sorted(d.glob('ale_shadow_r*_call1.bin')):
        with path.open('rb') as fh:
            header=np.fromfile(fh,dtype='<u8',count=4)
            assert tuple(header[[0,1,3]])==(0x474653414c453031,1,29)
            p=np.fromfile(fh,dtype=DTYPE,count=int(header[2]));assert not fh.read(1)
        assert len(p)==header[2] and np.all(p['flags'][:,0]==1)
        data.append(p)
    for path in sorted(d.glob('ale_cells_r*_call1.tsv')):
        cells.append(np.genfromtxt(path,names=True,delimiter='\t',ndmin=1))
    for path in sorted(d.glob('ale_meta_r*_call1.json')):metas.append(json.loads(path.read_text()))
    p=np.concatenate(data);c=np.concatenate(cells)
    order=np.lexsort((p['j'],p['i']));p=p[order];c=c[np.argsort(c['id'])]
    assert np.all(p['i']<p['j']) and np.all(np.diff(c['id'])>0)
    assert not np.any((p['i'][1:]==p['i'][:-1]) & (p['j'][1:]==p['j'][:-1]))
    a=np.searchsorted(c['id'],p['i']);b=np.searchsorted(c['id'],p['j'])
    assert np.array_equal(c['id'][a],p['i']) and np.array_equal(c['id'][b],p['j'])
    return p,c,a,b,metas


def accumulate(values,a,b,n):
    shape=(n,)+values.shape[1:];out=np.zeros(shape)
    np.add.at(out,a,-values);np.add.at(out,b,values)
    return out


def thermal(rates,c):
    return rates[:,3]-c['vx']*rates[:,1]-c['vy']*rates[:,2]+.5*(c['vx']**2+c['vy']**2)*rates[:,0]


def summarize(d):
    p,c,a,b,metas=load(d);n=len(c);z=p['a']
    rates={key:accumulate(z[:,start:start+4],a,b,n) for key,start in
           [('full',5),('pressure',9),('advection',13),('riemann_residual',17),('native_canonical',21),('affine_reference',25)]}
    u={key:thermal(value,c) for key,value in rates.items()}
    native=c['Edot_native']-c['vx']*c['pxdot_native']-c['vy']*c['pydot_native']
    dv=-accumulate(z[:,4],a,b,n)
    sx=z[:,3]-z[:,1];sy=-(z[:,2]-z[:,0])
    closure=np.column_stack([-accumulate(sx,a,b,n),-accumulate(sy,a,b,n)])
    defect=rates['full']-rates['pressure']-rates['advection']-rates['riemann_residual']
    out=dict(directory=str(d),n_cells=n,n_unique_pairs=len(p),
        clamped_pairs=int((p['flags'][:,1]!=0).sum()),
        HLL_fallback_pairs=int(((p['flags'][:,2] & (1<<5))!=0).sum()),
        pair_antisymmetry_exact=bool(np.all(z[:,5:9]+(-z[:,5:9])==0)),
        max_surface_closure=float(np.linalg.norm(closure,axis=1).max()),
        full_pressure_advection_residual_cell_max=float(abs(defect).max()),
        global_conservative_rate=rates['full'].sum(axis=0).tolist(),
        sum_absolute_conservative_rates=abs(rates['full']).sum(axis=0).tolist(),focus=[],
        native_reconstruction_note='Native clamps are retained and counted; pressure-only and full conservative solvers differ. Full difference is decomposed, never assigned solely to a missing mass term.')
    # Separate changed pressure closure/quadrature from the native work-speed
    # change at the actual pressure-ratio switch.
    area=np.hypot(sx,sy)
    pn=(z[:,22]*sx+z[:,23]*sy)/(area*area)
    native_geom=z[:,21:25].copy();native_geom[:,3]=pn*z[:,4]
    rg=accumulate(native_geom,a,b,n)
    pressure_closure=u['pressure']-thermal(rg,c)
    work_speed=thermal(rg,c)-u['native_canonical']
    for ident in FOCUS:
        hit=np.flatnonzero(c['id']==ident)
        if not len(hit):continue
        i=int(hit[0]);terms=dict(pressure_closure=float(pressure_closure[i]),
             native_work_speed=float(work_speed[i]),advection=float(u['advection'][i]),
             riemann_residual=float(u['riemann_residual'][i]),
             native_canonical_vs_actual=float(u['native_canonical'][i]-native[i]))
        out['focus'].append(dict(id=ident,U=float(c['U'][i]),rho=float(c['rho'][i]),P=float(c['P'][i]),
             native_Udot=float(native[i]),full_ALE_Udot=float(u['full'][i]),
             mdot=float(rates['full'][i,0]),Vdot_geom=float(dv[i]),
             rho_dot_full=float((rates['full'][i,0]-c['rho'][i]*dv[i])/c['V'][i]),
             rho_dot_fixed_mass=float(-c['rho'][i]*dv[i]/c['V'][i]),
             signed_full_minus_native_terms=terms,
             difference_identity_residual=float(u['full'][i]-native[i]-sum(terms.values()))))
    if metas[0]['affine_active']:
        coeff=np.array(metas[0]['affine']);rho,P,cx,cy,vx0,vy0,bxx,bxy,byx,byy=coeff
        interior=(abs(c['x']-.5)<.25)&(abs(c['y']-.5)<.25)
        ref=rates['affine_reference']
        mass_err=ref[:,0]-rho*dv
        scale=max(float(abs(ref[interior,0]).max()),float(abs(rho*dv[interior]).max()),1e-30)
        # Independent degree-3 polynomial surface integration (three Gauss
        # points), using analytic affine fluid states and native mesh motion.
        dx=(c['x'][b]-c['x'][a]+.5)%1-.5;dy=(c['y'][b]-c['y'][a]+.5)%1-.5
        dist=np.hypot(dx,dy);ex=dx/dist;ey=dy/dist;nx=sx/area;ny=sy/area
        dux=c['vx'][b]-c['vx'][a];duy=c['vy'][b]-c['vy'][a]
        expected=np.zeros((len(p),4));gamma=metas[0]['Gamma']
        for xi,weight in ((-np.sqrt(3/5),5/9),(0.,8/9),(np.sqrt(3/5),5/9)):
            ox=.5*(z[:,0]+z[:,2])+.5*xi*(z[:,2]-z[:,0])
            oy=.5*(z[:,1]+z[:,3])+.5*xi*(z[:,3]-z[:,1])
            delta=-(dux*(ox-.5*dx)+duy*(oy-.5*dy))/dist
            wx=.5*(c['vx'][a]+c['vx'][b])+delta*ex
            wy=.5*(c['vy'][a]+c['vy'][b])+delta*ey
            xx=c['x'][a]+ox-cx;yy=c['y'][a]+oy-cy
            vx=vx0+bxx*xx+bxy*yy;vy=vy0+byx*xx+byy*yy
            vn=vx*nx+vy*ny;relative=vn-wx*nx-wy*ny
            fm=rho*relative;eden=P/(gamma-1)+.5*rho*(vx*vx+vy*vy)
            physical=np.column_stack([fm,fm*vx+P*nx,fm*vy+P*ny,eden*relative+P*vn])
            expected+=.5*weight*area[:,None]*physical
        quadrature_diff=z[:,25:29]-expected
        qscale=np.maximum(np.max(abs(expected),axis=0),1e-30)
        speed=max(float(np.hypot(c['vx'],c['vy']).max()),float(np.hypot(vx0,vy0)),1.e-20)
        sound=np.sqrt(gamma*P/rho)
        natural_mass=rho*(2*speed+sound)*float(area.max())
        natural=np.array([natural_mass,natural_mass*speed+P*area.max(),
                          natural_mass*speed+P*area.max(),
                          natural_mass*(speed*speed+gamma*P/((gamma-1)*rho))])
        qerror=np.max(abs(quadrature_diff),axis=0)
        out['affine']=dict(coefficients=coeff.tolist(),certified_interior_cells=int(interior.sum()),
            actual_native_rho_range=[float(c['rho'].min()),float(c['rho'].max())],
            actual_native_P_range=[float(c['P'].min()),float(c['P'].max())],
            GCL_reference_mdot_minus_rho_Vdot_abs_max=float(abs(mass_err[interior]).max()),
            GCL_reference_relative_to_peak_mass_rate=float(abs(mass_err[interior]).max()/scale),
            Vdot_abs_max_interior=float(abs(dv[interior]).max()),
            fixed_mass_density_rate_abs_max=float(abs(c['rho'][interior]*dv[interior]/c['V'][interior]).max()),
            actual_native_reconstruction_GCL_error_abs_max=float(abs(rates['full'][interior,0]-rho*dv[interior]).max()),
            polynomial_quadrature_max_abs=np.max(abs(quadrature_diff),axis=0).tolist(),
            polynomial_quadrature_relative_to_peak=(qerror/qscale).tolist(),
            polynomial_quadrature_relative_to_natural_flux=(qerror/natural).tolist(),
            zero_rate_note='Peak-relative errors are uninformative when the expected component vanishes; gate uses 1e-9*peak +128*machine-epsilon*natural flux scale.')
        assert np.all(qerror<=1e-9*qscale+128*np.finfo(float).eps*natural),out['affine']
        if bxx==bxy==byx==byy==0:
            assert abs(ref[interior,0]).max()<1e-15
        else:
            assert abs(mass_err[interior]).max()<1e-13
    assert out['pair_antisymmetry_exact'] and out['max_surface_closure']<1e-12
    assert out['full_pressure_advection_residual_cell_max']<1e-10
    return out,p,c,rates


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=('patches','real'),required=True)
    args=parser.parse_args();base=ROOT/'rk4_diagnostics'/('ale_shadow_'+args.phase+'_v1')
    target=base/'analysis.json'
    if target.exists():raise SystemExit('Refuse to overwrite analysis')
    manifest=json.loads((base/'manifest.json').read_text())
    result=dict(scope='read-only conservative ALE shadow; no fixed-mass update or long stability proof',cases={})
    loaded={}
    for row in manifest['cases']:
        out,p,c,r=summarize(base/row['name']);result['cases'][row['name']]=out;loaded[row['name']]=(p,c,r)
    if args.phase=='patches':
        p,c,_=loaded['shear'];q,_,_=loaded['shear_boost']
        assert np.array_equal(p['i'],q['i']) and np.array_equal(p['j'],q['j'])
        boost=np.array([3.,-2.]);full=p['a'][:,5:9];expected=full.copy()
        expected[:,1:3]+=full[:,0,None]*boost
        expected[:,3]+=full[:,1:3]@boost+.5*np.dot(boost,boost)*full[:,0]
        delta=q['a'][:,5:9]-expected;scale=np.maximum(abs(expected).max(axis=0),1e-30)
        result['GPU_boost']=dict(velocity_boost=boost.tolist(),max_abs=abs(delta).max(axis=0).tolist(),
            relative_to_peak=(abs(delta).max(axis=0)/scale).tolist())
        q,_,_=loaded['shear_mpi4']
        assert np.array_equal(p['i'],q['i']) and np.array_equal(p['j'],q['j'])
        delta=q['a'][:,5:9]-p['a'][:,5:9];scale=np.maximum(abs(p['a'][:,5:9]).max(axis=0),1e-30)
        result['MPI1_MPI4']=dict(max_abs=abs(delta).max(axis=0).tolist(),relative_to_peak=(abs(delta).max(axis=0)/scale).tolist())
        assert max(result['GPU_boost']['relative_to_peak'])<1e-8
        assert max(result['MPI1_MPI4']['relative_to_peak'])<1e-8
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(target),cases={k:dict(n_cells=v['n_cells'],pairs=v['n_unique_pairs'],
        affine=v.get('affine'),focus=v['focus']) for k,v in result['cases'].items()},**{k:v for k,v in result.items() if k in ('GPU_boost','MPI1_MPI4')}),indent=2),flush=True)


if __name__=='__main__':main()

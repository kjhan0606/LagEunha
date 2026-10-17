"""Read-only field sensitivity on common cells plus a geometric strip bound.

Separate Eulerian probes can miss very narrow regions swept by a displaced
cell face. Integrate where cell labels agree, and bound the complementary
area using analytic envelopes of the numerical fields. Quadrature and
floating geometry remain numerical, not rigorous interval enclosures.
"""
import hashlib,json,os
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from physical_field_accuracy_audit import environment
from simple_wave_evolution import install_domain
from compare_physical_budget import compatible,scales_and_goals,endpoint_matches
from simple_wave_exact import polygon_rule
import moment_ale_evolve as model
from kinetic_factor_rhs import pressure_gradients
from quadratic_velocity_reconstruction import prepare as prepare_velocity,correction
from quadratic_pressure_reconstruction import prepare as prepare_pressure,value as pressure_value


def density_envelope(density,poly):
    x,y=np.max(abs((np.asarray(poly)-density.origin)/density.scale),axis=0);envelope=np.array([x,y,x*x,x*y,y*y]);variation=float(abs(density.theta)@envelope)
    logvalue=np.log(density.mass)-density.logZ+variation;slack=128*np.finfo(float).eps*(abs(np.log(density.mass))+abs(density.logZ)+variation+1.)
    return float(np.exp(logvalue+slack))


def velocity_envelope(state,poly):
    extent=np.max(abs(np.asarray(poly)-state['c']),axis=0);bound=abs(state['v'])+abs(state['B'])@extent;meta=state.get('velocity_quadratic')
    if meta is not None:
        x,y=np.max(abs((np.asarray(poly)-meta['center'])/meta['scale']),axis=0);psi=np.array([x*x,x*y,y*y])+np.array([1.,x,y])@abs(meta['projection']);bound+=abs(meta['coefficient'])@psi
    return float(np.linalg.norm(bound)*(1+128*np.finfo(float).eps))


def pressure_envelope(meta,poly):
    x,y=np.max(abs((np.asarray(poly)-meta['center'])/meta['scale']),axis=0)
    if meta.get('degree',2)>2:envelope=np.array([x**a*y**b for a,b in meta['powers']])+abs(meta['basis_means'])
    else:
        S=meta['second'];envelope=np.array([x,y,x*x+abs(S[0,0]),x*y+abs(S[0,1]),y*y+abs(S[1,1])])
    return float((abs(meta['mean'])+abs(meta['coefficient'])@envelope)*(1+128*np.finfo(float).eps))


def load(root):
    root=Path(root);cfg=json.loads((root/'simple_wave_config.json').read_text());cp=json.loads((root/'checkpoint.json').read_text());initial=json.loads((root/'initial_state.json').read_text());r=json.loads((root/'simple_wave_report.json').read_text());identity=hashlib.sha256((root/'simple_wave_config.json').read_bytes()).hexdigest()
    if cp['config_sha256']!=identity or initial['config_sha256']!=identity or r['config']!=cfg or not r['completed'] or not r['numerical_passed'] or not endpoint_matches(cp['t'],cfg['end']):raise ValueError('overlap audit source incomplete or inconsistent')
    np.testing.assert_array_equal(cp['mass'],initial['mass']);os.environ.update(environment(cfg));a=np.asarray(cp['state'])[:-4].reshape(-1,20);mass=np.asarray(cp['mass']);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t']);states=[model.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)]
    if os.environ['LAGFORCE_QUADRATIC_VELOCITY']=='1':prepare_velocity(states,cells,edges,mass,[np.asarray(v) for v in cp['factors']])
    centers,grad=pressure_gradients(states,edges);polynomials=None
    if os.environ['LAGFORCE_QUADRATIC_PRESSURE']=='1':polynomials,_=prepare_pressure(states,cells)
    caps=dict(density=0.,velocity=0.,pressure=0.)
    for i,(s,p) in enumerate(zip(states,cells)):
        caps['density']=max(caps['density'],density_envelope(s['density'],p));caps['velocity']=max(caps['velocity'],velocity_envelope(s,p))
        cap=pressure_envelope(polynomials[i],p) if polynomials is not None else abs(s['P'])+abs(grad[i])@np.max(abs(p-centers[i]),axis=0)
        caps['pressure']=max(caps['pressure'],float(cap))
    return dict(config=cfg,cp=cp,initial=initial,cells=cells,polygons=[Polygon(p) for p in cells],states=states,pressure=polynomials,centers=centers,grad=grad,caps=caps)


def fields(data,i,points):
    s=data['states'][i];velocity=s['v']+(points-s['c'])@s['B'].T+correction(s,points)
    pressure=s['P']+(points-data['centers'][i])@data['grad'][i] if data['pressure'] is None else pressure_value(data['pressure'][i],points)
    return dict(density=s['density'].value(points),velocity=velocity,pressure=pressure)


def intersections(left,right):
    common=[a.intersection(b) for a,b in zip(left,right)];strip=.5*sum(a.symmetric_difference(b).area for a,b in zip(left,right));return common,float(strip)


def run():
    install_domain();cfg=json.loads(Path('physical_budget_comparison_config.json').read_text());data=[load(r) for r in cfg['roots']]
    if len(data)!=2:raise ValueError('two overlap endpoints required')
    compatible([d['config'] for d in data],[d['initial'] for d in data]);scales,goals=scales_and_goals(data[0]['config']);common,strip=intersections(*[d['polygons'] for d in data]);rows=[]
    for order in [12,20]:
        totals={k:0. for k in scales}
        for i,p in enumerate(common):
            if p.is_empty or p.area==0:continue
            points,w=polygon_rule(np.asarray(p.exterior.coords[:-1]),order);left,right=[fields(d,i,points) for d in data]
            for k in scales:
                delta=left[k]-right[k];norm=np.linalg.norm(delta,axis=1) if delta.ndim==2 else abs(delta);totals[k]+=float(w@norm/scales[k])
        rows.append(dict(order=order,common_cell_L1_over_scale=totals))
    strip_bound={k:strip*(data[0]['caps'][k]+data[1]['caps'][k])/scales[k] for k in scales};upper={k:rows[-1]['common_cell_L1_over_scale'][k]+strip_bound[k]+abs(rows[-1]['common_cell_L1_over_scale'][k]-rows[0]['common_cell_L1_over_scale'][k]) for k in scales};fractions={k:upper[k]/goals[k] for k in scales}
    result=dict(completed=True,roots=cfg['roots'],quadratures=rows,different_cell_label_area=strip,field_envelopes=[d['caps'] for d in data],strip_contribution_upper_bound_over_scale=strip_bound,conditional_total_L1_sensitivity_upper_estimate=upper,fraction_of_physical_goal=fractions,fraction_of_goal_limit=cfg.get('fraction_of_goal_limit',.05),conditional_sensitivity_passed=bool(max(fractions.values())<cfg.get('fraction_of_goal_limit',.05)),partition_area_balance_error=float(abs(sum(p.area for p in data[0]['polygons'])-sum(p.area for p in common)-strip)),note='Same-label intersections integrated at12/20 positive Gaussian orders. On the different-label region, field difference bounded by the sum of analytic field envelopes times its area. Observed quadrature difference added as a numerical screen. Bounds are conditional on numerical geometry/quadrature; no rigorous interval enclosure is claimed. No field/state altered or evolution repeated.')
    atomic('physical_overlap_sensitivity_report.json',result);print(result,flush=True)


if __name__=='__main__':run()

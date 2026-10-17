"""Postprocess nine small 2D results; no evolution."""
import json
import sys
from pathlib import Path
import numpy as np


def analyze(directory):
    d=Path(directory); rows=[json.loads(s) for s in (d/'evolution.jsonl').read_text().splitlines()]
    assert len(rows)==9,len(rows)
    waves=sorted([r for r in rows if r['kind']=='wave' and r['boost']==[0.,0.] and r['refine']==1],key=lambda r:r['n'])
    errors=[r['cell_average_l1'][0] for r in waves]
    orders=[float(np.log2(a/b)) for a,b in zip(errors,errors[1:])]
    coarse=[r for r in rows if r['kind']=='wave' and r['n']==16 and r['refine']==1][0]
    fine=[r for r in rows if r['refine']==2][0]
    still=[r for r in waves if r['n']==8][0]
    boost=[r for r in rows if r['kind']=='wave' and r['boost']!=[0.,0.]][0]
    summary=dict(wave_errors=errors,orders=orders,amplitude_ratios=[r['amplitude_ratio'] for r in waves],
        phase_errors=[r['phase_error'] for r in waves],split_faces=[r['split_faces'] for r in waves],
        topology_changes=[r['topology_changes'] for r in waves],
        time_density_difference=float(np.max(abs(np.array(coarse['final_density'])-fine['final_density']))),
        time_l1_error_difference=abs(coarse['cell_average_l1'][0]-fine['cell_average_l1'][0]),
        boosted_density_difference=float(np.max(abs(np.array(still['final_density'])-boost['final_density']))),
        max_mass_flux=max(r['mass_flux'] for r in rows),max_energy_relative=max(abs(r['energy_relative']) for r in rows),
        max_momentum_error=max(r['momentum_error'] for r in rows),min_pressure=min(r['pressure_min'] for r in rows),
        min_area=min(r['area_min'] for r in rows),max_iterations=max(r['iterations'] for r in rows),
        max_edge_mass_flux=max(r['edge_mass_flux'] for r in rows),
        contact_density_error=max(r['cell_average_l1'][0] for r in rows if r['kind']=='contact'),
        contact_translation_error=max(r['translated_generator_error'] for r in rows if r['kind']=='contact'),
        contact_weight_change=max(r['weight_change'] for r in rows if r['kind']=='contact'))
    assert all(a>b for a,b in zip(errors,errors[1:]))
    assert summary['max_mass_flux']<3e-13 and summary['max_energy_relative']<1e-10
    (d/'analysis.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':analyze(sys.argv[1])

# Additional static diagnostic can be called explicitly after analyze().
def midpoint_diagnostic(directory):
    from periodic2d_hllc import initial,mesh,rhs,fan,integrate_face,sample
    d=Path(directory); rows=[json.loads(s) for s in (d/'evolution.jsonl').read_text().splitlines()]
    r=[r for r in rows if r['kind']=='wave' and r['n']==16 and r['refine']==1][0]
    mass,_=initial(16,'wave',np.array([0.,0.])); y=np.array(r['final_state']); dy,_=rhs(y,mass)
    _,area,edges=mesh(y[:,:2],y[:,2]); v=y[:,3:5]/mass[:,None]
    u=np.column_stack([mass/area,y[:,3:5]/area[:,None],y[:,5]/area])
    exact_change=np.zeros((len(y),4)); mid_change=np.zeros_like(exact_change); midpoint_max_face_error=0.
    for i,j,a,b,delta,length in edges:
        distance=np.linalg.norm(delta); normal=delta/distance; dv=v[j]-v[i]
        b0=(delta@v[j]-dv@a)/distance; b1=(delta@v[j]-dv@b)/distance
        offset=(dy[i,2]-dy[j,2])/(2*distance); f=fan(u[i],u[j],normal)
        exact_flux,_,_=integrate_face(f,b0,b1,offset,length)
        mid_flux=sample(f,(b0+b1)/2+offset)[0]*length
        midpoint_max_face_error=max(midpoint_max_face_error,float(np.max(abs(exact_flux-mid_flux))))
        exact_change[i]-=exact_flux;exact_change[j]+=exact_flux
        mid_change[i]-=mid_flux;mid_change[j]+=mid_flux
    result=dict(exact_mass_residual=float(np.max(abs(exact_change[:,0]))),
        midpoint_mass_residual_same_weights=float(np.max(abs(mid_change[:,0]))),
        midpoint_rhs_error=np.max(abs(exact_change-mid_change),axis=0).tolist(),max_face_flux_difference=midpoint_max_face_error)
    (d/'midpoint_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))

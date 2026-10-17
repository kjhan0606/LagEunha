"""Focused solver-tolerance check on the AF N4 high-range interior face."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
import face_density_range_audit as audit
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from adaptive_coupled_mass import atomic


def main():
    audit.OPTIONS.update(primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10)
    name='n4_r100';cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    rebuilt,_=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)]);low=float(means.min());high=float(means.max())
    previous=json.loads(Path('predecessor_completion.json').read_text())['cases'][name];rows=[]
    for level in (2,3):
        pred=next(r for r in previous['comparisons'] if r['level']==level and r['angle_degrees']==0.)
        for i,j in ((11,14),(14,11)):
            edge=next(e for e in edges if (e[0]==i and e[1]==j) or (e[0]==j and e[1]==i))
            p=audit.problem(cells[i],np.r_[mass[i],a[i,3:9]],level,pred['cells'][i]['prior'],means[i],low,high)
            weights,_=audit.trace_weights(p['parts'],edge[2],edge[3]);coefficient=weights*p['mass']/p['raw'][:,0]/high
            for budget in (0.,1e-10,1e-8):
                extrema=[audit.extremum(p,coefficient,budget,sign) for sign in (1.,-1.)]
                good=all(x.get('gate',False) for x in extrema)
                rows.append(dict(level=level,cell=i,neighbor=j,budget=budget,all_gates=good,extrema=extrema,
                    density_width=(extrema[1]['value']-extrema[0]['value'])*high if good else None))
    atomic('strict_report.json',dict(completed=True,all_gates=all(x['all_gates'] for x in rows),rows=rows,
        note='Posthoc check on AF worst-case face; stricter HiGHS 1e-10 tolerances. Zero budget remains a numerical, not exact-arithmetic, test.'))


if __name__=='__main__':main()

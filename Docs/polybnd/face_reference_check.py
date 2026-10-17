"""Independent inward-shift Shapely check of analytic one-sided face traces."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon,LineString
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from adaptive_coupled_mass import atomic


def main():
    report=json.loads(Path('completion_report.json').read_text());cases={}
    for name,result in report['cases'].items():
        cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
        _,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);a0=y0[:-4].reshape(-1,20)
        initial,_=model.geometry(a0[:,:2],a0[:,2],0.);refs=[]
        for v in initial:
            points=v.copy();points[:,0]+=cp['t']*(.2+1.1*points[:,1]);refs.append(Polygon(points))
        phases=np.where(a0[:,0]<.5,1.,cfg['amplitude'])
        a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
        rows=[r for r in result['records'] if r['level']==3 and r['budget']==1e-5];checks=[]
        for epsilon in (1e-10,1e-12):
            errors=[]
            for row in rows:
                i=row['cell'];offset=epsilon*np.sqrt(Polygon(cells[i]).area)
                for edge_id,exact in zip(row['face_ids'],row['exact_face_density']):
                    ei,ej,va,vb,normal,length=edges[edge_id]
                    if i==ej:normal=-normal
                    segment=LineString([va-offset*normal,vb-offset*normal])
                    value=sum(rho*segment.intersection(ref).length for rho,ref in zip(phases,refs))/length
                    errors.append(abs(value-exact)/max(1.,abs(exact)))
            checks.append(dict(epsilon=epsilon,max_relative_difference=float(max(errors))))
        cases[name]=checks
    atomic('reference_check.json',dict(completed=True,cases=cases,note='Inward offsets approximate analytic traces; errors should decrease with epsilon until floating-point limitations. Independent Shapely line intersections.'))


if __name__=='__main__':main()

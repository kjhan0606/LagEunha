"""Read-only diagnosis of phase extrema near stiff accepted path intervals."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
from adaptive_coupled_mass import atomic
import moment_ale_evolve as model


def run():
    source=json.loads(Path('checkpoint.initial.json').read_text());cfg=source['config'];mass=np.asarray(source['mass'])
    setup('shear_contact',cfg['n'],cfg['amplitude']);manifest=json.loads(Path('steps_manifest.json').read_text());rows=[]
    for entry in manifest:
        saved=json.loads(Path(entry['file']).read_text());record=saved['record'];a=np.asarray(record['y'])[:-4].reshape(-1,20)
        polys,_=model.geometry(a[:,:2],a[:,2],record['t']);means=np.array([m/Polygon(v).area for m,v in zip(mass,polys)])
        order=np.argsort(means);low,high=means[order[0]],means[order[-1]]
        d0=np.asarray(record['d0'])[:-4].reshape(-1,20);dm=np.asarray(record['dm'])[:-4].reshape(-1,20)
        accel=2*(dm[:,2]-d0[:,2])/record['dt'];local=np.array(saved['audit']['absolute'])/mass
        row=dict(step=entry['step'],t=record['t'],dt=record['dt'],means=means,low_indices=order[:3],high_indices=order[-3:],
            low_gap=means[order[1]]-low,high_gap=high-means[order[-2]],low=low,high=high,
            nearest_low_fraction=(means[order[1]]-low)/(high-low),nearest_high_fraction=(high-means[order[-2]])/(high-low),
            max_weight_acceleration=float(np.max(abs(accel))),weight_acceleration=accel,
            worst_mass_cell=int(np.argmax(local)),relative_step_absolute_mass_defect=float(local.max()))
        rows.append(row);print({k:row[k] for k in ('step','t','dt','low_gap','high_gap','max_weight_acceleration','worst_mass_cell')},flush=True)
    atomic('phase_extrema_report.json',dict(completed=True,rows=rows,
        note='Read-only accepted states. Extrema gaps are diagnostics, not proof of the source of stiffness. No density thresholds or evolution equations changed.'))


if __name__=='__main__':run()

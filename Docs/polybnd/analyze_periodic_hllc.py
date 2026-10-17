"""Analyze existing small JSON outputs; performs no evolution."""
import json
import sys
from pathlib import Path
import numpy as np
from periodic_hllc_convergence import initial,geometry


def analyze(directory):
    directory=Path(directory)
    rows=[json.loads(line) for line in (directory/'evolution.jsonl').read_text().splitlines()]
    assert len(rows)==18, len(rows)
    summary={}
    for kind,boost in [('contact',.7),('smooth_contact',.7),('wave',0.),('wave',.7)]:
        selected=sorted([r for r in rows if r['kind']==kind and r['boost']==boost and r['cfl']==.15],key=lambda r:r['n'])
        assert [r['n'] for r in selected]==[16,32,64,128]
        errors=[r['reconstructed_density_l1'] if kind=='smooth_contact' else r['cell_average_l1'][0] for r in selected]
        orders=None if kind=='contact' else [float(np.log2(a/b)) for a,b in zip(errors,errors[1:])]
        if kind!='contact': assert all(a>b for a,b in zip(errors,errors[1:]))
        summary[f'{kind}_{boost}']={'errors':errors,'orders':orders,'amplitude_ratios':[r['mode_amplitude_ratio'] for r in selected],
            'phase_errors':[r['mode_phase_error'] for r in selected],'speed_errors':[r['mode_speed_error'] for r in selected]}
    def density(row):
        m,_=initial(row['n'],row['kind'],row['boost']); y=np.array(row['final_state'])
        return m/geometry(y[:,0],y[:,1])[2]
    summary['time_refinement']={}
    for boost in (0.,.7):
        r=[r for r in rows if r['kind']=='wave' and r['boost']==boost and r['n']==128]
        summary['time_refinement'][str(boost)]={'max_density_difference':float(np.max(abs(density(r[0])-density(r[1])))),
            'density_l1_error_difference':abs(r[0]['cell_average_l1'][0]-r[1]['cell_average_l1'][0])}
    summary['galilean_density_difference']={}
    for n in (16,32,64,128):
        r=[r for r in rows if r['kind']=='wave' and r['cfl']==.15 and r['n']==n]
        summary['galilean_density_difference'][str(n)]=float(np.max(abs(density(r[0])-density(r[1]))))
    summary['gates']={
        'max_mass_flux':max(r['mass_flux'] for r in rows),
        'max_energy_relative':max(abs(r['energy_relative']) for r in rows),
        'max_momentum_error':max(abs(r['momentum_error']) for r in rows),
        'minimum_pressure':min(r['pressure_min'] for r in rows),
        'minimum_width':min(r['width_min'] for r in rows),
        'max_newton_evaluations':max(r['iterations'] for r in rows),
        'max_edge_mass_flux':max(r['edge_flux'] for r in rows)}
    (directory/'analysis.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__': analyze(sys.argv[1])

"""Static full-RHS regression of the optional continuous phase-floor map."""
import json
from pathlib import Path
import numpy as np
from predictive_coupling_batch import saved_case,metrics
from predictive_boundary_density import BoundaryDensity
from adaptive_coupled_mass import atomic
import kinetic_factor_rhs as kinetic


def run():
    rows=[]
    for name in ('n3_r10','n3_r100','n4_r100'):
        cp,mass,y,L=saved_case(name)
        for mode in ('neighbor','moment'):
            values=[];diags=[]
            for smooth in (False,True):
                provider=BoundaryDensity(mode,vector_geometry=True,smooth_phase_floor=smooth)
                dy,d=kinetic.rhs(cp['t'],y,mass,covariance_factors=L,face_density_provider=provider,face_segments=64,face_order=2,exact_bin_transport=True)
                values.append(dy);diags.append(metrics(dy,d,mass))
            difference=float(np.max(abs(values[0]-values[1])))
            good=difference<1e-11 and all(m['closure_moment_error']<1e-8 and m['mass_flux']<1e-10 and m['density_min']>0 and m['closure_volume_cap_excess']<1e-8 and m['pressure_rate']<1e-8 and max(m['affine_rhs_errors'].values())<1e-8 and m['ledger_rate_error']<1e-10 for m in diags)
            row=dict(case=name,mode=mode,passed=bool(good),max_rhs_difference=difference,metrics=diags)
            rows.append(row);atomic('smoothing_progress.json',dict(rows=rows));print(name,mode,difference,good,flush=True)
    probe=json.loads(Path('threshold_probe_report.json').read_text());near=next(r for r in probe['rows'] if r['multiplier']==1.01)
    probe_gate=probe['smooth_phase_floor'] and near['max_profile_difference']<2e-6 and near['normalized_frozen_weight_mass_rate']<1e-11
    report=dict(passed=all(r['passed'] for r in rows) and bool(probe_gate),rows=rows,threshold_probe_gate=bool(probe_gate),
        note='Optional C1 blend only within1..4 roundoff thresholds; full static old/smooth comparison and prior-map probe. No time evolution or universal regularity claim. Active AQ uses the original unsmoothed map.')
    atomic('smoothing_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()

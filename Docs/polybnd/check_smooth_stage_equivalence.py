"""Confirm the shared pressure-gradient refactor against an immutable RHS."""
import json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
import kinetic_factor_rhs as kinetic
from smooth_mass_path import initialize_worker


def run():
    initialize_worker();initial=json.loads(Path('initial_state.json').read_text())
    r=restore(json.loads(Path('step_1_trial.json').read_text())['record']);mass=np.asarray(initial['mass'])
    dy,diag=kinetic.rhs(r['t'],r['y'],mass,covariance_factors=r['L'],face_order=4)
    difference=float(np.max(abs(dy-r['d0'])))
    result=dict(passed=difference==0.,maximum_rhs_difference=difference,
        note='One immutable stage recomputed only to verify the pressure-gradient code extraction. Approved dynamical steps are reused.')
    atomic('stage_equivalence.json',result);print(result,flush=True)
    raise SystemExit(0 if result['passed'] else 2)


if __name__=='__main__':run()

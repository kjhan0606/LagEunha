"""One shared new RHS at a completed endpoint for independent dt comparisons."""
import hashlib,json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from simple_wave_evolution import preserved_initial,install_domain
import kinetic_factor_rhs as kinetic


def run():
    cfg=json.loads(Path('time_source_config.json').read_text());install_domain()
    source,mass,state,L=preserved_initial(cfg)
    dy,diag=kinetic.rhs(source['t'],state,mass,covariance_factors=L,face_order=4)
    passed=diag['mass_flux']<1e-10 and diag['min_pressure']>0 and np.all(np.isfinite(dy))
    atomic('shared_source_rhs.json',dict(t=source['t'],state_sha256=hashlib.sha256(state.tobytes()).hexdigest(),
        mass_sha256=hashlib.sha256(mass.tobytes()).hexdigest(),dy=dy,diag=diag))
    atomic('time_source_gate.json',dict(passed=bool(passed),source_time=source['t'],mass_flux=diag['mass_flux'],min_pressure=diag['min_pressure'],
        note='First new RHS at the preserved final state, shared by all dt branches. No accepted integration repeated.'))
    print('time source',passed,flush=True)
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':run()

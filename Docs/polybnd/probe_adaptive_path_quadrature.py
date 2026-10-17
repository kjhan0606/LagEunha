"""Refine a frozen rejected path; no fluid step or weight root is repeated."""
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import initialize_case,density_options,initialize_path_worker
from predictive_boundary_density import BoundaryDensity
from predictive_path_refinement import assess_adaptive


def run():
    cp=json.loads(Path('source_initial_config.json').read_text());mass,_=initialize_case(cp['config'])
    np.testing.assert_array_equal(cp['mass'],mass)
    rejected=json.loads(Path('frozen_rejection.json').read_text());options=density_options(cp['config'],True)
    seed=json.loads(Path('seed_cache.json').read_text())
    if not Path('adaptive_nodes.json').exists():atomic('adaptive_nodes.json',seed)
    original={str(i):seed['samples'][str(i*256)] for i in range(65)}
    with ProcessPoolExecutor(max_workers=24,mp_context=mp.get_context('spawn'),initializer=initialize_path_worker,
            initargs=(cp['config']['n'],cp['config']['amplitude'],True,True,options)) as executor:
        checked=assess_adaptive(rejected['record'],mass,1e-4,BoundaryDensity(**options),original,
            'adaptive_nodes.json','frozen rejected trial',executor,min_nodes=513)
    report=dict(completed=True,passed=bool(checked['ratio']<=1 and checked['pressure']>0 and checked['factor_R']>=-1e-20),
        source_time=rejected['t'],dt=rejected['dt'],old_ratio=rejected['ratio'],audit=checked,
        note='Same frozen path, mass/rate/quadrature budgets unchanged. Refine quadrature before reducing evolution dt. Finest nested difference after focus replacement; also control absolute-integral uncertainty. Old coarse corrections are diagnostics, not the error estimate of the finer replaced integral.')
    atomic('adaptive_quadrature_report.json',report);raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()

"""Explicit final numerical gates and duration-aware continuation diagnostics."""
import copy
import numpy as np


def budget_duration(cfg,source_time):
    """Allocate each accepted path against the declared full target.

    A shorter execution endpoint may be used for inspection, but cannot
    enlarge the per-step mass/quadrature allocation or relax the rate gate.
    Configuration hashes still protect continuation identity.
    """
    end=float(cfg['end']);target=float(cfg.get('validation_end',end))
    source=float(source_time)
    if not all(np.isfinite(x) for x in (end,target,source)) or end<=source or target<end:
        raise ValueError('finite validation endpoint must cover the execution interval')
    return target-source


def retarget_peak(peak,old_duration,new_duration):
    if not np.isfinite(old_duration+new_duration) or old_duration<=0 or new_duration<=0:raise ValueError('positive finite durations required')
    result=copy.deepcopy(peak);result['rate']*=new_duration/old_duration
    return result


def numerical_gates(accum,mass,peak,cfg):
    mass=np.asarray(mass)
    passed=bool(np.max(np.asarray(accum['absolute'])/mass)<=cfg['path_mass_budget'] and np.max(abs(np.asarray(accum['quad_delta']))/mass)<=cfg['path_quadrature_budget'])
    if cfg.get('adaptive_path_quadrature'):
        passed=passed and all(np.max(np.asarray(accum[k])/mass)<=cfg['path_quadrature_budget'] for k in ('quad_error_bound','absolute_quad_error_bound'))
    checks=dict(integrated_mass_and_quadrature=bool(passed),path_rate=peak['rate']<=cfg['path_rate_budget'],generator_path=peak['generator']<=cfg['generator_path_budget'],ledger=peak['ledger']<1e-10,stage_mass=peak['stage_mass']<1e-10,positive_pressure=peak['min_pressure']>0,positive_factor=peak['min_factor_R']>=-1e-20)
    return {k:bool(v) for k,v in checks.items()}

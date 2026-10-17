"""Nested time quadrature on an exact integer lattice; old samples are retained."""
import numpy as np
from adaptive_coupled_mass import simpson


def requires_refinement(checked,cfg,adaptive=True,integrated_geometry=False):
    """Resolve quadrature before treating its absolute integral as a failure.

    A sampled rate or generator violation is already a genuine pointwise
    failure. An unresolved quadrature estimate of the mass integral is not.
    No budget is changed by this decision.
    """
    if not adaptive or checked['ratio']<=1:return False
    if checked['rate']>cfg['path_rate_budget'] or checked['generator']>cfg['generator_path_budget']:return False
    if not integrated_geometry and not cfg.get('resolve_mass_quadrature_before_reject',False) and checked['physical_ratio']>1:return False
    return checked['quadrature_ratio']>1


def node_indices(count,lattice=128):
    if lattice not in (128,1024) or count not in (9,17,33,65,129,257,513,1025) or count-1>lattice:raise ValueError('unsupported smooth path node count')
    return list(range(0,lattice+1,lattice//(count-1)))


def upgrade_nine_samples(samples):
    if any(str(int(k))!=k or not 0<=int(k)<=8 for k in samples):raise ValueError('not a legacy nine-node cache')
    return {str(int(k)*16):v for k,v in samples.items()}


def ensure_sample_lattice(trial,pending_path,sample_path,target):
    """Crash-recoverable exact integer reindexing; no sample recomputation."""
    import json,hashlib
    from pathlib import Path
    from adaptive_coupled_mass import atomic
    pending_path=Path(pending_path);sample_path=Path(sample_path)
    old=trial.get('sample_lattice')
    if target not in (128,1024) or old not in (128,1024) or target<old:raise ValueError('unsupported path lattice migration')
    if old==target and 'sample_lattice_upgrade' not in trial:
        return trial,json.loads(sample_path.read_text()) if sample_path.exists() else {}
    if 'sample_lattice_upgrade' not in trial:
        samples=json.loads(sample_path.read_text()) if sample_path.exists() else {}
        from predictive_evolution_batch import record_digest
        archive=sample_path.parent/('samples_before_lattice_'+record_digest(trial['record'])+'.json')
        atomic(archive,samples)
        trial['sample_lattice_upgrade']=dict(source=old,target=target,archive=archive.name,sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
        atomic(pending_path,trial)
    migration=trial['sample_lattice_upgrade'];archive=sample_path.parent/migration['archive']
    if migration['target']!=target or hashlib.sha256(archive.read_bytes()).hexdigest()!=migration['sha256']:raise ValueError('path lattice archive mismatch')
    samples=json.loads(archive.read_text());factor=target//migration['source']
    if any(str(int(k))!=k or not 0<=int(k)<=migration['source'] for k in samples):raise ValueError('invalid cached path index')
    updated={str(int(k)*factor):value for k,value in samples.items()}
    atomic(sample_path,updated);trial['sample_lattice']=target;trial.pop('sample_lattice_upgrade');atomic(pending_path,trial)
    return trial,updated


def assess_values(values,mass,dt,duration,cfg,check_absolute_quadrature=False):
    flux=np.asarray([v['flux'] for v in values]);mass=np.asarray(mass)
    signed=simpson(flux,dt);absolute=simpson(abs(flux),dt)
    quad=signed-simpson(flux[::2],dt);absolute_quad=absolute-simpson(abs(flux[::2]),dt)
    rate=duration*float(np.max(abs(flux)/mass));gen=max(v['diag']['generator_defect'] for v in values)
    pressure=min(v['diag']['min_pressure'] for v in values);gram=min(v['diag']['min_factor_R'] for v in values)
    mass_ratio=float(np.max(absolute/mass))/(cfg['path_mass_budget']*dt/duration)
    signed_quad_ratio=float(np.max(abs(quad)/mass))/(cfg['path_quadrature_budget']*dt/duration)
    absolute_quad_ratio=float(np.max(abs(absolute_quad)/mass))/(cfg['path_quadrature_budget']*dt/duration)
    physical_ratio=max(mass_ratio,rate/cfg['path_rate_budget'],gen/cfg['generator_path_budget'])
    quadrature_ratio=max(signed_quad_ratio,absolute_quad_ratio) if check_absolute_quadrature else signed_quad_ratio
    return dict(signed=signed,absolute=absolute,quad_delta=quad,absolute_quad_delta=absolute_quad,
        rate=rate,generator=gen,pressure=pressure,gram=gram,mass_ratio=mass_ratio,
        physical_ratio=physical_ratio,signed_quadrature_ratio=signed_quad_ratio,
        absolute_quadrature_ratio=absolute_quad_ratio,quadrature_ratio=quadrature_ratio,
        ratio=max(physical_ratio,quadrature_ratio))

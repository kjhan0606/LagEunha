"""Scan each accepted lineage prefix once; never splice in an unaccepted suffix."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029')
BASE = ROOT/'native_u64_kdk_individual_20261104_v1/kdk_block2_wallkdk_v1'
ADAPTER = BASE/'analysis_tools/full_reference_audit_20261104_v1/audit.py'
spec = importlib.util.spec_from_file_location('reference_audit', ADAPTER)
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


def scan(log, start_step, end_step, start_time, end_time, rk4):
    digest = hashlib.sha256()
    records = dict(time_rows=0, energy_rows=0, wall_rows=0, impacts=0,
                   rk4_pressure_torque=0., collision_torque=0.,
                   collision_torque_absolute_by_record=0., kinetic_roundoff=0.)
    bad = []
    previous_step = start_step
    last_exact_time = None
    cutoff_bytes = 0
    with log.open('rb') as stream:
        for raw in stream:
            digest.update(raw)
            cutoff_bytes += len(raw)
            line = raw.decode(errors='replace')
            is_rk4_energy = '[RK4E]' in line
            is_kdk_energy = '[U64_KDK_AUDIT]' in line
            if is_rk4_energy or is_kdk_energy or '[PHASE1E]' in line:
                fields = dict(re.findall(r'\b([A-Za-z][A-Za-z0-9_]*)=([0-9.eE+-]+)', line))
                guards = ('floor_inj','floor_cum','de_inj','de_cum','sfl_inj','sfl_cum',
                          'cum_inj','cum_floor','clip_inj','clip_ev','n_clip_all',
                          'n_ie_le0','invalid','clips','floors')
                for name in guards:
                    if name in fields and float(fields[name]) != 0:
                        bad.append(dict(field=name,value=fields[name],line=line.strip()))
                if is_kdk_energy:
                    last_exact_time = float(fields['t'])
                    for name in ('Umin','Vmin'):
                        if not (math.isfinite(float(fields[name])) and float(fields[name])>0):
                            bad.append(dict(field=name,value=fields[name],line=line.strip()))
                if is_rk4_energy or is_kdk_energy:
                    records['energy_rows'] += 1
            if '[WALL_ACCEPT]' in line:
                records['wall_rows'] += 1
                values = dict(re.findall(r'\b(t_start|h|impacts|kinetic_roundoff)=([0-9.eE+-]+)',line))
                last_exact_time = float(values['t_start']) + float(values['h'])
                records['impacts'] += int(values['impacts'])
                records['kinetic_roundoff'] += float(values['kinetic_roundoff'])
                for tag, key in [('wall_pressure_J','rk4_pressure_torque'),
                                 ('wall_collision_J','collision_torque')]:
                    match = re.search(tag+r'=\(([^)]+)\)',line)
                    if not match:
                        raise ValueError('Incomplete accepted RK4 wall ledger')
                    torque = float(match[1].split(',')[-1])
                    records[key] += torque
                    if key == 'collision_torque':
                        records['collision_torque_absolute_by_record'] += abs(torque)
            if '[KDK_WALL_LEDGER]' in line or '[KDK_LTS_WALL_LEDGER]' in line:
                records['wall_rows'] += 1
                values = dict(re.findall(r'\b(hits|wall_torque|kinetic_roundoff)=([0-9.eE+-]+)',line))
                records['impacts'] += int(values['hits'])
                torque = float(values['wall_torque'])
                records['collision_torque'] += torque
                records['collision_torque_absolute_by_record'] += abs(torque)
                records['kinetic_roundoff'] += float(values['kinetic_roundoff'])
            if line.startswith('Time='):
                match = re.search(r'\bstep=\s*(\d+)',line)
                if not match:
                    raise ValueError('Missing absolute accepted step')
                step = int(match[1])
                if step != previous_step + 1:
                    raise ValueError(f'Accepted-step discontinuity: {previous_step}->{step}')
                records['time_rows'] += 1
                previous_step = step
                if step == end_step:
                    break
        else:
            raise ValueError(f'No accepted cutoff step {end_step} in {log}')
    if records['energy_rows'] != records['time_rows']:
        raise ValueError('Energy ledger does not cover every accepted step')
    if rk4 and records['wall_rows'] != records['time_rows']:
        raise ValueError('RK4 wall ledger does not cover every accepted step')
    if last_exact_time is None or abs(last_exact_time-end_time)>1e-10:
        raise ValueError(f'Exact ledger/footer endpoint mismatch: {last_exact_time}, {end_time}')
    return dict(log=str(log),start_step=start_step,end_step=end_step,
                start_time=start_time,end_time=end_time,last_exact_time=last_exact_time,
                prefix_bytes=cutoff_bytes,prefix_sha256=digest.hexdigest(),
                records=records,nonzero_or_invalid_guards=bad,
                pressure_wall_torque_complete=rk4,
                accepted_cell_ID_wall_history_complete=rk4)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('endpoint',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    t = adapter.gate.snapshot_time(args.endpoint)
    endpoint_step = int(args.endpoint.name.split('.')[1])
    cases = [
        (ROOT/'native_u64_w0_to20pi_20261101_v1/run/log',1,164092,0.,40.00009144754977,True),
        (BASE/'run_wall_global_t4002/log',164092,164175,40.00009144754977,40.02,False),
        (BASE/'run_wall_global_t4010_adiabatic/log',164175,164511,40.02,40.1,False),
        (BASE/'run_wall_global_t4020_roundoff/log',164511,164947,40.1,40.2,False),
        (BASE/'run_wall_lts_t4030_eta01/log',164947,165209,40.2,40.3,False),
        (args.endpoint.parent/'log',165209,endpoint_step,40.3,t,False)]
    rows = [scan(*case) for case in cases]
    result = dict(endpoint=str(args.endpoint),time=t,segments=rows,
                  all_accepted_steps_guarded=all(not row['nonzero_or_invalid_guards'] for row in rows),
                  accepted_steps=sum(row['records']['time_rows'] for row in rows),
                  impacts=sum(row['records']['impacts'] for row in rows),
                  known_collision_torque=sum(row['records']['collision_torque'] for row in rows),
                  known_RK4_pressure_torque=sum(row['records']['rk4_pressure_torque'] for row in rows),
                  KDK_pressure_wall_torque_complete=False,
                  KDK_unique_impacted_cell_mass_complete=False,
                  qualification='Exact producer/consumer checkpoints are separately SHA-bound. '
                                'No overlapping suffix or duplicate warmup is included. Guard checks '
                                'cover accepted steps; incomplete pressure/ID wall accounting is '
                                'explicit and is not assigned a fictitious zero value.')
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='segments'},indent=2))


if __name__ == '__main__':
    main()

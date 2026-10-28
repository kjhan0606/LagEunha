"""Provenance + corruption checks on a real completed native checkpoint."""
import argparse
import copy
import json
from pathlib import Path
from full_moment_restart_provenance import seal, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint', type=Path, required=True)
    p.add_argument('--initial', type=Path, required=True)
    p.add_argument('--physics', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(exist_ok=False)
    manifest = seal(args.checkpoint, args.initial, args.physics, args.output/'sealed')
    record = json.loads(manifest.read_text())
    header, _ = verify(args.checkpoint, manifest, args.initial, args.physics, record['h'])
    assert header == (args.checkpoint/'header.bin').read_bytes()
    changes = dict(physics_sha256='0'*64, initial_state_sha256='0'*64,
                   operator_sha256='0'*64, integrator='RK4', fixed_mass=False,
                   Re_infinity=False, steering=True, time=record['time']+.001,
                   original_run_token=record['original_run_token']+1)
    rejected = {}
    for key, value in changes.items():
        bad = copy.deepcopy(record)
        bad[key] = value
        path = args.output/(key+'.json')
        path.write_text(json.dumps(bad))
        try:
            verify(args.checkpoint, path, args.initial, args.physics, record['h'])
        except AssertionError:
            rejected[key] = True
        else:
            raise AssertionError('Invalid restart accepted: '+key)
    for key in ('state_files', 'evidence'):
        bad = copy.deepcopy(record)
        first = next(iter(bad[key]))
        bad[key][first] = '0'*64
        path = args.output/(key+'.json')
        path.write_text(json.dumps(bad))
        try:
            verify(args.checkpoint, path, args.initial, args.physics, record['h'])
        except AssertionError:
            rejected[key] = True
        else:
            raise AssertionError('Corrupt restart digest accepted: '+key)
    try:
        verify(args.checkpoint, manifest, args.initial, args.physics, record['h']/2)
    except AssertionError:
        rejected['silent_timestep_change'] = True
    else:
        raise AssertionError('Restart timestep change accepted')
    report = dict(status='passed_native_restart_provenance', manifest=str(manifest),
        header_exact=True, rejected=rejected, physical_providers=0,
        accepted_time_steps=0, runtime_physical_restart_verified=False)
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()

"""Compare cached physical audits only at matching actual evolved times."""
import json,argparse
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,default=Path('/gpfs/kjhan/LagForce/Codex_smooth_qr_physical_to_t10_20261018'))
parser.add_argument('--case',default='smooth_from_t2')
args=parser.parse_args();root=args.root
manifest = json.loads((root / 'manifest.json').read_text())
start = manifest['cases'][args.case]['start_time']
cache = root / 'physical_frames'
legacy = []
for name in ('legacy_history_to_t4', 'legacy_continue'):
    for path in sorted((cache / name).glob('khout.*.dat.json')):
        legacy.append(json.loads(path.read_text()))
dest = root / 'physical_comparisons'
dest.mkdir(exist_ok=True)
summary = {'matched': [], 'awaiting_matching_legacy': [], 'source_frames_excluded': []}
for path in sorted((cache / args.case).glob('khout.*.dat.json')):
    smooth = json.loads(path.read_text())
    if smooth['time'] <= start + 1e-10:
        summary['source_frames_excluded'].append(smooth['time'])
        continue
    matches = [row for row in legacy if abs(row['time'] - smooth['time']) < 1e-10]
    if not matches:
        summary['awaiting_matching_legacy'].append(smooth['time'])
        continue
    # A copied continuation source may duplicate an original history frame.
    assert len({row['sha256'] for row in matches}) == 1, 'Ambiguous same-time legacy states'
    old = matches[0]
    assert smooth['mass_unchanged'] and old['mass_unchanged']
    assert smooth['weights_zero'] and old['weights_zero']
    assert smooth['regions'].keys() == old['regions'].keys()
    regions = {}
    for name, values in smooth['regions'].items():
        reference = old['regions'][name]
        assert values['cells'] == reference['cells']
        regions[name] = {
            key: {'legacy': reference[key], 'smooth': value,
                  'smooth_over_legacy': value / reference[key] if reference[key] else None}
            for key, value in values.items() if key != 'cells'
        }
    result = {'smooth_time': smooth['time'], 'legacy_time': old['time'],
              'smooth_snapshot': smooth['snapshot'], 'legacy_snapshot': old['snapshot'],
              'smooth_sha256': smooth['sha256'], 'legacy_sha256': old['sha256'],
              'regions': regions,
              'scope': 'Same actual time, original-IC material masks and errors. Smooth evolution starts at t2 with inherited error. Ratios are physical error comparisons, not stability certification. Phase is wrapped; inspect position and velocity as well.'}
    output = dest / (path.stem + '.comparison.json')
    if output.exists():
        assert json.loads(output.read_text()) == result, 'Cached comparison differs'
    else:
        with output.open('x') as stream:
            json.dump(result, stream, indent=2)
    summary['matched'].append({'time': smooth['time'], 'output': str(output), 'regions': regions})
print(json.dumps(summary, indent=2))

#!/usr/bin/env python3
"""Copy current source/small evidence; leave native datasets on shared GPFS."""
import hashlib, json, os, shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

CONTROL = Path(__file__).resolve().parent
CFG = json.loads((CONTROL / 'config.json').read_text())
ROOT = Path(CFG['root'])
PUB = Path(CFG['publish'])
ARC = PUB / 'Research/Kepler'
CODE = {'.c', '.cpp', '.h', '.cu', '.inc', '.py', '.sh', '.slurm',
        '.sbatch', '.md', '.json', '.csv', '.sha256', '.patch'}
TEXT = {'.json', '.sha256', '.txt', '.log', '.tsv', '.sh', '.sbatch', '.py'}


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def copy_record(pair):
    src, dst = pair
    for attempt in range(3):
        before = src.stat()
        data = src.read_bytes()
        after = src.stat()
        if before.st_size == after.st_size == len(data) and before.st_mtime_ns == after.st_mtime_ns:
            break
    else:
        raise RuntimeError('Source kept changing during shutdown snapshot: ' + str(src))
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists() or dst.read_bytes() != data:
        dst.write_bytes(data)
    dst.chmod(0o755 if src.suffix in {'.sh', '.slurm', '.sbatch'} else 0o644)
    return {'source': str(src), 'repository_path': str(dst.relative_to(PUB)),
            'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def make_snapshot():
    assert PUB.is_dir() and (PUB / '.git').exists()
    manifest = ARC / 'SNAPSHOT_SOURCE_MANIFEST.json'
    previous = json.loads(manifest.read_text())
    history = ARC / 'history/pre_poweroutage_518f357'
    history.mkdir(parents=True, exist_ok=True)
    for name in ['SNAPSHOT_SOURCE_MANIFEST.json', 'SNAPSHOT_VERIFICATION.json',
                 'SNAPSHOT_CONTENTS.sha256']:
        if (ARC / name).exists() and not (history / name).exists():
            shutil.copyfile(ARC / name, history / name)
    pairs = {}

    def take(src, dst, limit=20 * 2**20):
        if src.is_file() and not src.is_symlink() and src.stat().st_size <= limit:
            pairs[str(dst)] = (src, dst)

    # The previously published original solver files retain their exact mapping.
    for rec in previous['files']:
        if rec['repository_path'].startswith('Exam/'):
            take(Path(rec['source']), PUB / rec['repository_path'])
    for srcdir, destdir in [(ROOT, ARC / 'rk4_diagnostics'),
                            (ROOT.parent, ARC / 'stream_candidates')]:
        for src in srcdir.iterdir():
            if src.suffix in CODE or src.name == 'AGENTS.md':
                take(src, destdir / src.name)

    bundles = ['material_cell_reference_lifetime_v1_ada_mpi4',
               'material_cell_reference_lifetime_v2_ada_mpi4',
               'material_cell_reference_continuation_v1_h200_mpi4',
               'material_cell_reference_lifetime_build_evidence_v1',
               'material_reference_resource_block_audit_20261023_v1',
               'material_step6_reference_continuation.claim',
               'material_cell_condensation_capacity_v1_h200']
    for name in bundles:
        base = ROOT / name
        if not base.exists():
            continue
        for src in base.iterdir():
            if src.suffix in TEXT | CODE:
                take(src, ARC / 'evidence' / name / src.name, 2**20)
        for source_name in ['src', 'src_revision2']:
            srcdir = base / source_name
            if srcdir.is_dir():
                for src in srcdir.iterdir():
                    if src.suffix in CODE:
                        take(src, ARC / 'frozen' / name / source_name / src.name)
        if (base / 'run').is_dir():
            for src in (base / 'run').rglob('*'):
                if src.suffix in TEXT:
                    take(src, ARC / 'evidence' / name / src.relative_to(base), 2**20)

    # Latest geometry programs/scripts and compact result summaries. Large edge
    # tables, binary128 packets and repeated frozen source snapshots stay on GPFS.
    geometry = ROOT / 'material_compact_actual_source_v1'
    if geometry.is_dir():
        for src in geometry.iterdir():
            if src.is_file() and src.suffix in CODE | {'.log', '.txt'}:
                take(src, ARC / 'geometry' / src.name, 2**20)
            elif src.is_dir():
                for child in src.iterdir():
                    if child.is_file() and child.suffix in CODE | {'.log', '.txt'}:
                        take(child, ARC / 'geometry' / child.relative_to(geometry), 2**20)

    for src in CONTROL.iterdir():
        if src.suffix in {'.py', '.sh'}:
            take(src, ARC / 'operations/poweroutage_20261023' / src.name)
    config_public = {k: v for k, v in CFG.items() if k != 'ssh_auth_sock'}
    (ARC / 'operations/poweroutage_20261023').mkdir(parents=True, exist_ok=True)
    (ARC / 'operations/poweroutage_20261023/config.json').write_text(
        json.dumps(config_public, ensure_ascii=False, indent=2) + '\n')

    with ThreadPoolExecutor(max_workers=16) as pool:
        fresh = list(pool.map(copy_record, pairs.values()))
    records = {x['repository_path']: x for x in previous['files']}
    records.update({x['repository_path']: x for x in fresh})
    latest = dict(previous)
    latest.update(snapshot_date='2026-10-23T10:50:00+09:00',
                  snapshot_kind='user_requested_pre_poweroutage_stop',
                  note='Current source and small evidence snapshot. Large native datasets remain on GPFS. Historical status/pending fields describe their original observation time; read the latest handover.',
                  files=sorted(records.values(), key=lambda x: x['repository_path']))
    outcome = ROOT / 'material_cell_reference_continuation_v1_h200_mpi4/run/continuation_result.json'
    if outcome.exists():
        result = json.loads(outcome.read_text())
        latest.update(accepted_time=result['accepted_time'],
                      new_accepted_steps_after_step5=result['new_accepted_steps'])
    manifest.write_text(json.dumps(latest, ensure_ascii=False, indent=2) + '\n')
    for rec in latest['files']:
        assert digest(PUB / rec['repository_path']) == rec['sha256'], rec['repository_path']
    verification = {'scope': 'copied source/archive byte verification, no new physical test',
                    'records_checked': len(records), 'fresh_files_copied': len(fresh),
                    'repository_bytes_match_manifest': True,
                    'long_time_Kepler_stability_proven': False}
    (ARC / 'SNAPSHOT_VERIFICATION.json').write_text(json.dumps(verification, indent=2) + '\n')
    # Rebuild the old index because current development files have changed.
    content_files = [p for p in ARC.rglob('*') if p.is_file()
                     and p.name != 'SNAPSHOT_CONTENTS.sha256'
                     and 'history' not in p.relative_to(ARC).parts]
    content_files += [PUB / 'Docs/HANDOVER_2026-10-23.md',
                      PUB / 'Docs/POWER_OUTAGE_HANDOVER_2026-10-23.md']
    with ThreadPoolExecutor(max_workers=16) as pool:
        lines = list(pool.map(lambda p: digest(p) + '  ' + str(p.relative_to(PUB)) + '\n',
                              sorted(set(content_files))))
    (ARC / 'SNAPSHOT_CONTENTS.sha256').write_text(''.join(lines))
    return verification

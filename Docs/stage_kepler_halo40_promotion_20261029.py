"""Prepare an unselected candidate chain and finite A100 validation inputs."""
import ast
import datetime
import json
from pathlib import Path
import shutil

from kepler_halo40_candidate_chain_20261029 import stage_chain, CANDIDATE_SHA, sha
from run_kepler_halo40_promotion_20261029 import make_case

DOCS = Path(__file__).resolve().parent
Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
MAIN = Q1/'halo40_from_t2654918_a100_20261028'
SOURCE = Q1/'vertex_halo_pair_t2654918_a100_20261029'
ROOT = Q1/'halo40_promotion_validation_t3_a100_20261029'
CHAIN = Q1/'halo40_reuse_alias_from_t3_a100_20261029'


def main():
    candidate = Q1/'vertex_halo_native_build_20261029/eunha2.vertex.exe'
    canonical = (MAIN/'to3p0/eunha2.exe').resolve()
    canonical_sha = '6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0'
    if sha(canonical) != canonical_sha or sha(candidate) != CANDIDATE_SHA:
        raise ValueError('Native executable identity changed')
    ada = json.loads((Q1/'halo40_net_wall_ada_t2654918_20261029/net_wall_comparison.json').read_text())
    if not ada['numerical_passed'] or not ada['performance']['clear_improvement']:
        raise ValueError('Qualified Ada numerical/performance prerequisite missing')
    template = stage_chain(MAIN, CHAIN, candidate)
    ROOT.mkdir()
    helpers = ['run_kepler_halo40_promotion_20261029.py',
               'kepler_halo40_candidate_chain_20261029.py',
               'audit_kepler_halo40_net_wall_20261029.py',
               'audit_kepler_restart_bytes_20261029.py',
               'audit_kepler_padding_reuse_pair_20261028.py',
               'audit_kepler_geometry_force_dump_20261029.py',
               'audit_kepler_polar_long_gate_20261027.py',
               'audit_kepler_native_radial_20261026.py']
    for name in helpers:
        ast.parse((DOCS/name).read_text())
        shutil.copy2(DOCS/name, ROOT/name)
    guard = 'run_kepler_checkpoint_validation_20261029.py'
    shutil.copy2(MAIN/guard, ROOT/guard)
    (ROOT/'canonical.exe').symlink_to(canonical)
    (ROOT/'candidate.exe').symlink_to(CHAIN/'eunha2.exe')
    (ROOT/'initial.dat').symlink_to((SOURCE/'control40/khout.000000.dat').resolve())
    shutil.copy2(SOURCE/'control40/params.dat', ROOT/'params.dat')
    (ROOT/'templates').mkdir()
    for kind, case in [('canonical', 'qualifier_old'), ('control', 'control40')]:
        shutil.copy2(SOURCE/case/'run.sh', ROOT/'templates'/f'{kind}.sh')
    checkpoint = (SOURCE/'control40/khout.023008.dat').resolve()
    static = [('fresh_reference40', 'off', 4), ('combined40', 'combo', 4),
              ('old_a', 'canonical', 0), ('combined_b', 'combo', 0),
              ('combined_b2', 'combo', 0), ('old_a2', 'canonical', 0)]
    for name, kind, dumps in static:
        make_case(ROOT, name, kind, checkpoint, 23008, 2.66, dumps)
    files = helpers+[guard, 'canonical.exe', 'candidate.exe', 'initial.dat', 'params.dat',
                     'templates/canonical.sh', 'templates/control.sh']
    for name, _, _ in static:
        files += [name+'/run.sh', name+'/inputs.sha256']
    config = {'parent_job': '438957', 'stage': str(MAIN/'to3p0'), 'root': str(ROOT),
              'candidate_root': str(CHAIN), 'candidate_template': template,
              'expected_next_job': str(MAIN/'to4p0/job.slurm'),
              'cached_control': str(SOURCE/'control40'),
              'candidate_binary_sha256': CANDIDATE_SHA, 'canonical_binary_sha256': canonical_sha,
              'static_cases': [name for name, _, _ in static],
              'frozen_hashes': {name: sha(ROOT/name) for name in files}}
    config_path = ROOT/'config.json'
    config_path.write_text(json.dumps(config, indent=2)+'\n')
    receipt = {'observed_at': datetime.datetime.now().astimezone().isoformat(),
               'root': str(ROOT), 'candidate_root': str(CHAIN),
               'config': str(config_path), 'config_sha256': sha(config_path),
               'parent_job': '438957', 'production_selected': False,
               'original_main_checker_or_next_job_modified': False,
               'new_Slurm_jobs': 0, 'maximum_auxiliary_seconds': 3000,
               'required_allocation_seconds_remaining': 3600,
               'native_cases': '5 full cases if qualified same-allocation reference reusable, otherwise 6; plus 3 natural-step bridges',
               'output_estimated_GB': 1.4,
               'bridge': 'EXACT_STOP=0 and no diagnostic dt cap; nextafter controls loop exit, not RK4 h; checkpoint first actual accepted production step',
               'promotion': 'Only after actual A100 halo40 operator/full-state/performance/natural-step bridge and verified cleanup; atomic next_job CAS, no successor submission here',
               'fallback': 'Original next4 on ordinary failure after cleanup; operational stop if residual steps or interrupted transaction are unverified',
               'stops': [4,5,6,6.5,7,8,9,10]}
    (ROOT/'staging_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (DOCS/'KEPLER_HALO40_PROMOTION_STAGE_2026-10-29.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()

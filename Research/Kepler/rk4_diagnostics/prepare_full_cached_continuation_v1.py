"""Freeze native continuation only; never reinitialize mass or evolve in Python."""
import hashlib
import json
import pathlib
import shutil

root = pathlib.Path(__file__).resolve().parent
donor = root / 'material_full_GPU_mass_v1_h200_mpi4'
dest = root / 'material_full_cached_continuation_v1_h200_mpi4'
assert json.loads((donor / 'run/full_GPU_mass_continuation_result.json').read_text())['first_GPU_mass_continuation_pass']
assert (donor / 'native_exit_code.txt').read_text().strip() == '0'
dest.mkdir()
shutil.copytree(donor / 'Exam', dest / 'Exam')
for name in ('material_full_initial_driver.cpp', 'material_full_initial_batch.h',
             'material_full_initial_cuda.cu', 'material_full_evolved_cuda.inc',
             'material_full_cached_cuda.inc', 'material_full_cached_host.h',
             'material_full_evolution_host.h', 'material_full_gpu_mass_host.h',
             'material_full_gpu_mass.h', 'material_full_gpu_mass.inc'):
    shutil.copy2(root / name, dest / 'Exam' / name)
for name in ('run_env.sh', 'immutable_payload.sha256', 'retained_build_objects.sha256'):
    shutil.copy2(donor / name, dest / name)
shutil.copy2(root / 'run_full_cached_continuation_v1.sh', dest / 'run.sh')
shutil.copy2(donor / 'link_command.txt', dest / 'donor_link_command.txt')
retained = list((donor / 'run/first_GPU_mass_step').glob('accepted_r*.bin'))
retained += [donor / 'run/first_GPU_mass_step/step.json', donor / 'run/full_GPU_mass_continuation_result.json',
             donor / 'native_exit_code.txt', donor / 'predecessor_adoption.json', donor / 'binary.sha256']
assert len(retained) == 9
def manifest(paths):
    return ''.join(hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() + '  ' + str(p) + '\n' for p in paths)
(dest / 'accepted_input.sha256').write_text(manifest(retained))
sources = sorted(p for p in (dest / 'Exam').iterdir() if p.suffix in {'.c', '.cpp', '.cu', '.h', '.inc', '.sh'})
sources += [dest / 'run.sh', dest / 'run_env.sh', dest / 'immutable_payload.sha256',
            dest / 'accepted_input.sha256', dest / 'donor_link_command.txt']
(dest / 'source.sha256').write_text(manifest(sources))
(dest / 'execution_contract.json').write_text(json.dumps({
    'scope': 'native full trajectory continuation with bounded same-state point cache',
    'resume_time': 5.8133562848030723e-5,
    'resume_prefix': str(donor / 'run/first_GPU_mass_step/accepted_r'),
    'target': 1, 'wall_seconds': 86400,
    'stop_on_measured_cost_exceeding_wall_to_target': True,
    'point_cache_capacity_per_rank': 1048576,
    'point_cache_size_bytes_per_point': 568,
    'GPU_cards': 1, 'MPI_ranks': 4, 'OMP_threads_per_rank': 3,
    'new_independent_IC_or_mass_integral': False,
    'first_new_RHS_one_tile_per_rank_compares_original_stream': True,
    'full_map_injectivity_guard_complete': False,
    'long_Kepler_validated': False
}, indent=2) + '\n')
(dest / 'status.txt').write_text('FROZEN_BUILD_REQUIRED_NO_NEW_PHYSICAL_EXECUTION\n')
print(dest)

from pathlib import Path
import hashlib,json,os,subprocess,time
root=Path(__file__).parent
build=Path('/gpfs/kjhan/LagForce/Codex_full_disk_mpi_bridge_build_v2_20261019')
binary=build/'eunha2.exe'
expected=(build/'binary.sha256').read_text().strip()
assert hashlib.sha256(binary.read_bytes()).hexdigest()==expected
manifest=dict(scope='Existing native eunha2 MPI entrypoint + owned CUDA stage, immutable captured-input integration gate; no accepted time evolution',binary=str(binary),binary_sha256=expected,Slurm_job=os.environ.get('SLURM_JOB_ID'),nodes=os.environ.get('SLURM_JOB_NODELIST'),allocated_GPUs=2,allocated_CPUs=16,memory_GiB=64,input_cases=[1,2],physical_providers=0,evolution_steps=0,live_Ada_runs_modified=False)
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
env=dict(os.environ,EUNHA_MOMENT_BINARY_SHA=expected);cases=[]
for ranks in (1,2):
 begin=time.perf_counter()
 with (root/f'rank{ranks}.log').open('x') as log:
  process=subprocess.run(['mpirun','-np',str(ranks),'-ppn',str(ranks),str(binary),str(root/f'params_rank{ranks}.dat')],env=env,stdout=log,stderr=subprocess.STDOUT)
 row=dict(ranks=ranks,exit=process.returncode,wall_seconds=time.perf_counter()-begin)
 (root/f'rank{ranks}'/'exit.json').write_text(json.dumps(row,indent=2));cases.append(row)
 if process.returncode:raise RuntimeError(row)
 row.update(json.loads((root/f'rank{ranks}'/'native_result.json').read_text()))
 assert row['status']=='passed_native_owned_MPI_component'
assert (root/'rank1/stage.bin').read_bytes()==(root/'rank2/stage.bin').read_bytes()
assert (root/'rank1/Q.bin').read_bytes()==(root/'rank2/Q.bin').read_bytes()
assert (root/'rank1/hydro.bin').read_bytes()==(root/'rank2/hydro.bin').read_bytes()
t1,t2=[row['MPI_stage_wall_seconds'] for row in cases]
report=dict(status='passed_native_executable_MPI_component',cases=cases,rank1_rank2_outputs_bitwise_equal=True,same_H200_node_speedup=t1/t2,stage_parallel_efficiency=t1/t2/2,allocated_GPU_hours=2*sum(row['wall_seconds'] for row in cases)/3600,physical_providers=0,evolution_steps=0,MPI_evolution_verified=False,limits='Captured stage only, includes canonical flux distribution/owned CUDA/gather. Global flux oracle outside timed stage is additional. No density/velocity service or actual time integration in this gate; not a full solver speedup.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)

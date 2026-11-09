from pathlib import Path
import os,json,shutil,hashlib,struct,subprocess,numpy as np,h5py
R=Path('/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1');B=R/'Research/MaterialVolume_20261105';D=B/'material_pressure_acceleration_20261107_v1';O=D/'kh_t10_approved_20261109_v1';O.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
Q=D/'common_closed_volume_20261109_v1/full_kh128';G=B/'material_kh_128x256_20261106_v1/gizmo_t135';S=B/'gizmo_matched_kh_20261106_v1/source'
(O/'bin').mkdir();shutil.copy2(Q/'native.exe',O/'bin/native128.exe');shutil.copy2(D/'common_portable_20261109_v1/native/native.exe',O/'bin/native256.exe');shutil.copy2(S/'GIZMO',O/'bin/GIZMO');shutil.copy2(S/'Config.sh',O/'bin/GIZMO_Config.sh')
assert sha(O/'bin/native128.exe')=='c3554c696ec0b2044ac34ba20f54c4dfe7c920bf142eb6fd3e67a91f799a9ffb'
assert sha(O/'bin/native256.exe')=='11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e'
extract=B/'material_kh_stop_20261106_v1/extract_ids'
for nx,start in [(128,1.35),(256,0.)]:
 p=O/f'n{nx}';p.mkdir();(p/'input').mkdir();ic=p/'input/native.dat'
 if nx==128:
  source=Q/'native_t135/khout.001808.dat';shutil.copy2(source,ic);shutil.copy2(Q/'initial_ids.bin',p/'input/t0_ids.bin')
 else:
  source=B/'kh_unseeded_20261106_v1/initial'
  with (p/'input/ic_generation.json').open('x') as f:subprocess.run([str(source),str(nx),str(2*nx),'1','2','kh',str(5/3),str(ic),'dynamic','.01'],check=True,stdout=f)
 subprocess.run([str(extract),str(ic),str(p/'input/native_ids.bin')],check=True,stdout=subprocess.DEVNULL)
 a=np.fromfile(p/'input/native_ids.bin').reshape(-1,13);a=a[a[:,12].argsort()]
 assert len(a)==nx*2*nx and np.isfinite(a).all() and (a[:,7:9]>0).all()
 if nx==256:shutil.copy2(p/'input/native_ids.bin',p/'input/t0_ids.bin')
 if nx==128:
  with h5py.File(G/'output/snap_001.hdf5') as f:
   attrs=dict(f['Header'].attrs);data={k:f['PartType0'][k][:] for k in ['Coordinates','Velocities','Masses','InternalEnergy','ParticleIDs']}
  z=np.fromfile(D/'gizmo_sync128_t135.bin').reshape(-1,10);z=z[z[:,7].argsort()];order=np.argsort(data['ParticleIDs']);data={k:v[order] for k,v in data.items()}
  assert np.array_equal(data['ParticleIDs'],z[:,7]) and np.max(abs(data['Coordinates'][:,:2]-z[:,:2]))<1e-12
  data['Velocities'][:,:2]=z[:,2:4];data['InternalEnergy']=z[:,5];assert np.array_equal(data['Masses'],z[:,4])
  synced_source=D/'gizmo_sync128_t135.bin'
 else:
  n=len(a);num=np.array([n,0,0,0,0,0],dtype=np.uint32)
  attrs=dict(NumPart_ThisFile=num,NumPart_Total=num,NumPart_Total_HighWord=np.zeros(6,dtype=np.uint32),MassTable=np.zeros(6),Time=0.,Redshift=0.,BoxSize=1.,NumFilesPerSnapshot=1,Omega0=0.,OmegaLambda=0.,HubbleParam=1.,Flag_Sfr=0,Flag_Cooling=0,Flag_StellarAge=0,Flag_Metals=0,Flag_Feedback=0,Flag_DoublePrecision=1,Flag_Entropy_ICs=0)
  x=np.zeros((n,3));x[:,:2]=a[:,:2];v=np.zeros((n,3));v[:,:2]=a[:,2:4]
  data=dict(Coordinates=x,Velocities=v,Masses=a[:,6],InternalEnergy=a[:,7]/a[:,6],ParticleIDs=a[:,12].astype(np.uint64)+1);synced_source=ic
 with h5py.File(p/'input/snap_001.hdf5','x') as f:
  h=f.create_group('Header')
  for k,v in attrs.items():h.attrs[k]=v
  g=f.create_group('PartType0')
  for k,v in data.items():g[k]=v
 assert np.array_equal(a[:,12]+1,data['ParticleIDs']) and np.array_equal(a[:,6],data['Masses'])
 params=(Q/'native_t135/params.dat').read_text()
 import re
 for key,value in [('Nx',nx),('Ny',2*nx)]:params,n=re.subn(r'(?m)^(define\s+'+key+r'\s*=\s*)\d+',lambda m:m[1]+str(value),params);assert n==1
 (p/'input/params.dat').write_text(params)
 e=json.loads((Q/'native_t135/environment.json').read_text());e.update(GFS_MATERIAL_DT=str(.001*128/nx),EUNHA_DUMP_DT='20')
 (p/'input/environment.json').write_text(json.dumps(e,indent=2))
 plan=dict(nx=nx,ny=2*nx,cells=2*nx*nx,start=start,target=10.,outputs=[2.1,4.7,9.2,10.],source=str(source),source_sha256=sha(source),native_input_sha256=sha(ic),gizmo_input_sha256=sha(p/'input/snap_001.hdf5'),gizmo_export_source=str(synced_source),gizmo_export_sha256=sha(synced_source),native_mpi_ranks=1,native_omp_threads=8,gizmo_mpi_ranks=4,gizmo_omp_threads_per_rank=2,physical_method='Existing fixed-mass common endpoint-volume SSPRK2; C=0.2; no physical/model changes',native_binary=str(O/f'bin/native{nx}.exe'),native_binary_sha256=sha(O/f'bin/native{nx}.exe'),gizmo_binary_sha256=sha(O/'bin/GIZMO'),restart_scope='128: both resume t=1.35; GIZMO rank change uses synchronized-state snapshot restart (flag 2), resets time-bin hierarchy. 256: authorized spatial refinement from same continuous t=0 IC.',native_checkpoint_bytes=ic.stat().st_size,estimated_max_storage_bytes=4*1024**3 if nx==128 else 12*1024**3,free_bytes=shutil.disk_usage(O).free)
 (p/'plan.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(plan),flush=True)
(O/'authorization.json').write_text(json.dumps(dict(date='2026-11-09',scope='User approved KH and matching GIZMO to t=10, intermediate density/conservation, increased MPI and a larger simultaneous Slurm calculation. Stop after this stage; no other tests or method development.',mpi_limitation='Native material driver explicitly requires NID=1; do not bypass guard. GIZMO increased to4. Native OMP increased to8; larger native uses one H200.'),indent=2))
shutil.copy2(__file__,O/'prepare.py')

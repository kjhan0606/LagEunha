from pathlib import Path
import numpy as np,json,sys,re
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_moving_geometry_20261017")

def read(n,stem,cols):
 a=np.concatenate([np.loadtxt(root/f'np{n}'/f'{stem}.rank{r:04d}') for r in range(n)])
 assert a.shape==(36864,cols) and np.isfinite(a).all();a=a[np.argsort(a[:,0])];assert np.array_equal(a[:,0],np.arange(36864));return a
out={'scope':'HG, diagnostic velocity drift, actual periodic_postStage migration, rebuild, HG; not a full RK step','ranks':{}}
base_after=read(1,'raw_faces.dat.repeat02.hgafter',13)
for n in [1,2,4]:
 assert (root/f'np{n}/exit.code').read_text().strip()=='0'
 pre=read(n,'raw_faces.dat.repeat01.motion_before',37);post=read(n,'raw_faces.dat.repeat01.motion_after',37)
 assert np.array_equal(pre[:,4:],post[:,4:]),'migration changed non-position state or RK fields'
 dt=float(json.loads((root/f'np{n}/environment.json').read_text())['GFS_MPI_PROBE_DRIFT']);frozen=np.hypot(pre[:,2]-3,pre[:,3]-3)>=2.5
 expected=pre[:,2:4].copy();expected[~frozen]+=dt*pre[~frozen,5:7];expected=np.fmod(expected+6,6)
 position_error=float(np.max(np.abs(post[:,2:4]-expected)));assert position_error<2e-14
 changed=int(np.count_nonzero(pre[:,1]!=post[:,1]));assert n==1 or changed>0,'no actual ownership crossing'
 a=read(n,'raw_faces.dat.repeat02.hgafter',13);diff=np.abs(a[:,3:]-base_after[:,3:]);bits=np.ascontiguousarray(a[:,3:]).view(np.uint64)!=np.ascontiguousarray(base_after[:,3:]).view(np.uint64)
 assert np.all(a[:,2]==2)
 cache=read(n,'raw_faces.dat.repeat02.cache',25);reference=read(1,'raw_faces.dat.repeat02.cache',25)
 out['ranks'][str(n)]={'ownership_changed_particles':changed,'all_25_RK4_scalars_preserved':True,'non_position_state_preserved_in_transport':True,'position_update_max_abs_error':position_error,'post_rebuild_HG_bitwise_different_particles':int(np.count_nonzero(np.any(bits,axis=1))),'post_rebuild_HG_max_abs_difference':dict(zip(['x','y','mass','vx','vy','ie','volume','density','pressure','csound'],np.max(diff,axis=0).tolist())),'reconstruction_max_abs_difference':float(np.max(np.abs(cache[:,3:]-reference[:,3:])))}
(root/'migration_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

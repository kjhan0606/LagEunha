from pathlib import Path
import json,sys,struct
import numpy as np
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_hg_modes_20261017")

def read(path):
 with path.open('rb') as f:h=struct.unpack('<4q',f.read(32))
 magic,nm,ns,ncol=h;assert magic==20261017 and ncol==9 and nm>0 and ns>0
 assert path.stat().st_size==32+nm*9*8+ns*3*8
 m=np.memmap(path,mode='r',dtype='<f8',offset=32,shape=(nm,9));c=np.memmap(path,mode='r',dtype='<f8',offset=32+nm*9*8,shape=(ns,3))
 assert np.isfinite(m).all() and np.isfinite(c).all() and np.all(m[1:,1]>m[:-1,1])
 assert np.all(m[:,4:6]==0) and np.all(m[:,3]>0)
 assert np.array_equal(m[:,2],np.r_[0,np.cumsum(m[:-1,3])]) and int(m[:,3].sum())==ns
 return m,c
base,bc=read(root/'np1/raw_faces.dat.hgm.rank0000');out={'scope':'HG mode assembly with actual two-hop ghost state; no mode application, global ordering or evolution','reference_modes':len(base),'reference_support_entries':len(bc),'ranks':{}}
for n in [1,2,4]:
 assert (root/f'np{n}/exit.code').read_text().strip()=='0'
 keys=[];reports=[]
 for rank in range(n):
  m,c=read(root/f'np{n}/raw_faces.dat.hgm.rank{rank:04d}');keys.append(np.array(m[:,1]));ix=np.searchsorted(base[:,1],m[:,1]);assert np.all(ix<len(base)) and np.array_equal(base[ix,1],m[:,1]);assert np.array_equal(base[ix,3],m[:,3])
  owned=np.loadtxt(root/f'np{n}/raw_faces.dat.cache.rank{rank:04d}')[:,0];assert np.isin(m[:,0],owned).all()
  cols=[0,1,3,4,5,6,7,8];md=np.abs(m[:,cols]-base[ix][:,cols]);offset=np.repeat((base[ix,2]-m[:,2]).astype(np.int64),m[:,3].astype(np.int64));ridx=np.arange(len(c))+offset;ref=bc[ridx]
  iddiff=int(np.count_nonzero(c[:,0]!=ref[:,0]));cd=np.abs(c[:,1:]-ref[:,1:])
  reports.append({'rank':rank,'modes':len(m),'support_entries':len(c),'max_metadata_abs_difference':np.max(md,axis=0).tolist(),'support_id_order_mismatches':iddiff,'coefficient_different_rows':int(np.count_nonzero(np.any(cd!=0,axis=1))),'max_coefficient_abs_difference':np.max(cd,axis=0).tolist(),'dump_bytes':(root/f'np{n}/raw_faces.dat.hgm.rank{rank:04d}').stat().st_size})
  del ref,offset,ridx,md,cd,m,c
 keys=np.sort(np.concatenate(keys));assert np.array_equal(keys,base[:,1]),'mode union incomplete/duplicate'
 out['ranks'][str(n)]=reports
(root/'mode_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

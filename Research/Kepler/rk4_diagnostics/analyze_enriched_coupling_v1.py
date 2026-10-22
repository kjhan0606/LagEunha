from pathlib import Path
import argparse,json,struct
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--version',default='v5');args=p.parse_args()
dest=ROOT/f'rk4_diagnostics/enriched_{args.version}';manifest=json.loads((dest/'manifest.json').read_text())
projection=np.dtype(dict(names=['id','geometry','scale','L','mean'],
    formats=['<u8',('<f8',6),('<f8',4),('<f8',(10,10)),('<f8',10)],offsets=[0,24,88,120,920],itemsize=1528))
coupled=np.dtype(dict(names=['id','xy','mu','mass','z','rhs','geo'],
    formats=['<u8',('<f8',2),'<f8','<f8','<f8',('<f8',(4,10)),('<f8',(4,10))],offsets=[0,8,24,32,40,48,368],itemsize=688))
pairs=np.dtype(dict(names=['id','generator','corners','ledger'],
    formats=[('<u8',2),('<f8',2),('<f8',4),('<f8',10)],offsets=[0,16,32,64],itemsize=144))
out=dict(scope='actual initial native moving-cell Euler RHS and weight-rate audit; no evolution',cases=[]);arrays={};projections={}
def stable_sum(a):return float(np.sum(a,dtype=np.longdouble))
for case in manifest['cases']:
    name=case['name'];d=dest/name
    if not (d/'enriched_coupled.bin').is_file():continue
    report=json.loads((d/'enriched_coupled.json').read_text())
    with (d/'enriched_coupled.bin').open('rb') as f:h=struct.unpack('<8Q',f.read(64));perturb=struct.unpack('<2d',f.read(16))
    assert h[:2]==(0x4746534543523031,1) and h[3:6]==(10,4,688)
    a=np.fromfile(d/'enriched_coupled.bin',dtype=coupled,offset=80,count=h[2]);assert len(np.unique(a['id']))==len(a)
    pr=[];face=[]
    for rank in range(case['ranks']):
        path=d/f'enriched_project_r{rank}.bin'
        with path.open('rb') as f:hh=struct.unpack('<8Q',f.read(64))
        assert hh[5]==projection.itemsize
        pr.append(np.fromfile(path,dtype=projection,offset=64,count=hh[2]))
        path=d/f'enriched_pairs_r{rank}.bin'
        with path.open('rb') as f:hh=struct.unpack('<4Q',f.read(32))
        assert hh[:2]==(0x4746534543503031,1) and hh[3]==144
        face.append(np.fromfile(path,dtype=pairs,offset=32,count=hh[2]))
    b=np.concatenate(pr);b=b[np.argsort(b['id'])];face=np.concatenate(face)
    assert np.array_equal(a['id'],b['id'])
    keys=face['id'];assert len(np.unique(keys,axis=0))==len(keys)
    ii=np.searchsorted(a['id'],keys[:,0]);jj=np.searchsorted(a['id'],keys[:,1]);assert np.array_equal(a['id'][ii],keys[:,0])
    mean_face=np.zeros((len(a),4));vol=np.zeros(len(a));mass_at_w0=np.zeros(len(a))
    for c in range(4):np.add.at(mean_face[:,c],ii,-face['ledger'][:,c]);np.add.at(mean_face[:,c],jj,face['ledger'][:,c])
    np.add.at(vol,ii,face['ledger'][:,4]);np.add.at(vol,jj,-face['ledger'][:,4])
    np.add.at(mass_at_w0,ii,-face['ledger'][:,6]);np.add.at(mass_at_w0,jj,face['ledger'][:,6])
    g=b['geometry'];L=b['L'];mean=b['mean'];scale=b['scale'];area=g[:,0]
    rate=np.sqrt(2*scale[:,3]/scale[:,0])/np.minimum(g[:,3],g[:,5]);normalization=area[:,None]*scale*rate[:,None]
    defect=a['rhs']-a['geo'];normalized=np.max(abs(defect)/normalization[:,:,None],axis=(1,2))
    ax=g[:,3]*L[:,1,1];ay=g[:,4]*L[:,1,1]+g[:,5]*L[:,2,1];by=g[:,5]*L[:,2,2]
    x0=a['xy'][:,0]+g[:,1]+g[:,3]*L[:,1,0]+ax*mean[:,1]
    y0=a['xy'][:,1]+g[:,2]+g[:,4]*L[:,1,0]+g[:,5]*L[:,2,0]+ay*mean[:,1]+by*mean[:,2]
    center=.5 if name in ('uniform','rigid','shear_boost','varying_rho_isobaric') else 12.
    r=a['rhs'];Jdot=(x0-center)*r[:,2,0]+ax*r[:,2,1]-(y0-center)*r[:,1,0]-ay*r[:,1,1]-by*r[:,1,2]
    total_J=stable_sum(Jdot);face_J=stable_sum(face['ledger'][:,9]-face['ledger'][:,8])
    result=dict(name=name,report=report,mean_rates=[stable_sum(a['rhs'][:,c,0]) for c in range(4)],
        mass_ledger_max_error=float(np.max(abs(a['rhs'][:,0,0]-mean_face[:,0])/normalization[:,0])),
        max_mass_weight_rate_error=float(np.max(abs(a['rhs'][:,0,0]+a['mu']*a['z'])/a['mu'])),
        initial_W0_max_abs_mass_rate=float(np.max(abs(mass_at_w0))),
        max_abs_weight_rate=float(np.max(abs(a['z']))),
        angular_rate=total_J,periodic_face_angular_transport=face_J,
        angular_ledger_absolute_error=abs(total_J-face_J),
        max_normalized_rhs_minus_geometry=float(normalized.max()),focus=[])
    if name=='uniform':
        result['uniform_GCL_mass_error']=float(np.max(abs(a['rhs'][:,0,0]-a['mu']*vol)/normalization[:,0]))
    if name.startswith('fine'):
        radius=np.linalg.norm(a['xy']-12,axis=1);result['radial']=[]
        for lo,hi in [(0,.047),(.047,.1),(.1,.2),(.2,.5),(.5,1),(1,2),(2,8),(8,10),(10,18)]:
            z=(radius>=lo)&(radius<hi)
            if np.any(z):result['radial'].append(dict(lo=lo,hi=hi,cells=int(z.sum()),max=float(normalized[z].max()),median=float(np.median(normalized[z]))))
        for pid in (67943,68000,68527,68584,67935,68001,68526,68592):
            k=int(np.searchsorted(a['id'],pid));result['focus'].append(dict(id=pid,weight_rate=float(a['z'][k]),
                mass_rate=float(a['rhs'][k,0,0]),max_normalized_rhs_minus_geometry=float(normalized[k])))
    arrays[name]=a;projections[name]=b;out['cases'].append(result)
if 'fine_fresh' in arrays and 'fine_mpi4' in arrays:
    a,c=arrays['fine_fresh'],arrays['fine_mpi4'];b=projections['fine_fresh'];assert np.array_equal(a['id'],c['id'])
    rate=np.sqrt(2*b['scale'][:,3]/b['scale'][:,0])/np.minimum(b['geometry'][:,3],b['geometry'][:,5])
    normal=b['geometry'][:,0,None]*b['scale']*rate[:,None]
    out['MPI1_4']=dict(max_weight_rate_difference=float(np.max(abs(a['z']-c['z']))),
        max_normalized_RHS_difference=float(np.max(abs(a['rhs']-c['rhs'])/normal[:,:,None])),
        max_normalized_geometry_difference=float(np.max(abs(a['geo']-c['geo'])/normal[:,:,None])))
if 'fine_fresh' in arrays and 'fine_perturb' in arrays:
    a,c=arrays['fine_fresh'],arrays['fine_perturb'];b=projections['fine_fresh'];bx,by=.013,-.007
    transformed=a['rhs'][:,:,0].copy();old=transformed.copy()
    transformed[:,1]+=bx*old[:,0];transformed[:,2]+=by*old[:,0]
    transformed[:,3]+=bx*old[:,1]+by*old[:,2]+.5*(bx*bx+by*by)*old[:,0]
    rate=np.sqrt(2*b['scale'][:,3]/b['scale'][:,0])/np.minimum(b['geometry'][:,3],b['geometry'][:,5])
    normal=b['geometry'][:,0,None]*b['scale']*rate[:,None]
    out['uniform_velocity_perturbation']=dict(max_weight_rate_difference=float(np.max(abs(a['z']-c['z']))),
        max_normalized_mean_rate_transform_error=float(np.max(abs(c['rhs'][:,:,0]-transformed)/normal)))
(dest/'coupled_analysis.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

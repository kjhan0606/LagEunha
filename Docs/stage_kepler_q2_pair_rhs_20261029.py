"""Stage a zero-physical-time native RHS probe at a saved Q2 checkpoint."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_t3_support_20261029 import orbit

docs=Path(__file__).resolve().parent
root=Path('/gpfs/kjhan/LagForce/Kepler_Q2_PairRHS_Ada_20261029')
initial=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/t001_v2/khout.000000.dat')
endpoint=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1/khout.011058.dat')
binary=Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_v2/eunha2.exe')
params=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/params.dat')
template=docs/'KEPLER_Q1_BOUNDED_GLOBAL_2026-10-29/rhs_t10_probe_run.sh'
expected={initial:'a4d557ae5837e933c8c2c8c335f91e5c97d36f9a3bdb8acd5568b0499140b665',endpoint:'ceb49b0881577489ec51a93e4259c458256e4f16d8f4adc2b342c3656bee304b',binary:'f26816122df101ea71f72b4f3fb9233e965d06784ed94e88cfb46245fa251444'}
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def replace_once(s,a,b):
    if s.count(a)!=1:raise ValueError('Unexpected template: '+a)
    return s.replace(a,b)
def main():
    if root.exists():raise FileExistsError(root)
    for p,h in expected.items():
        if sha(p)!=h:raise ValueError('Changed source '+str(p))
    a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK;p,meta=ordered_global(endpoint,ids)
    if ma['t']!=0 or meta['t']!=47.50031653561342 or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Unexpected native state')
    xy=np.column_stack((p['x'],p['y']));xy0=np.column_stack((a['x'],a['y']))
    r=np.hypot(xy0[:,0]-3.,xy0[:,1]-3.);sel=np.flatnonzero((r>=.5)&(r<.75))
    d,j=cKDTree(xy).query(xy[sel],k=2,workers=1);exact=orbit(xy0,meta['t'])
    ratio=d[:,1]/np.linalg.norm(exact[sel]-exact[j[:,1]],axis=1)
    chosen=[]
    for k in np.argsort(ratio):
        pair=(int(ids[sel[k]]),int(ids[j[k,1]]))
        if pair[0] in chosen or pair[1] in chosen:continue
        chosen.extend(pair)
        if len(chosen)==6:break
    if chosen!=[31088,30249,29637,23949,25319,23161]:raise ValueError('Unexpected selected pairs')
    root.mkdir();(root/endpoint.name).symlink_to(endpoint);(root/'eunha2.exe').symlink_to(binary);shutil.copy2(params,root/'params.dat')
    selected=[{'id':uid,'purpose':'close_pair','initial_radius':float(r[np.searchsorted(ids,uid)])} for uid in chosen]
    cfg={'initial':str(initial),'endpoint':str(endpoint),'endpoint_sha256':expected[endpoint],
         'binary':str(binary),'binary_sha256':expected[binary],'restart_step':11058,'time':meta['t'],'selected':selected,
         'scope':'Single selected K1 native RHS evaluation; no physical evolution. Uses tested wall-event v2 binary, not historical Q2 executable.'}
    (root/'probe.json').write_text(json.dumps(cfg,indent=2)+'\n')
    s=template.read_text();s=replace_once(s,'cd /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/rhs_t10_selected_v1','cd '+str(root))
    s=replace_once(s,'EUNHA_KEPLER_DTETA=0','EUNHA_KEPLER_DTETA=0.06')
    s=replace_once(s,'HYDRO_TSTOP=10.0001','HYDRO_TSTOP=47.5004')
    s=replace_once(s,'GFS_RK4_AUDIT_STEP_MIN=20790 GFS_RK4_AUDIT_STEP_MAX=20790','GFS_RK4_AUDIT_STEP_MIN=11058 GFS_RK4_AUDIT_STEP_MAX=11058')
    s=replace_once(s,'GFS_RK4_AUDIT_IDS=73920,73540,73532,74301,76610,72774,65495,64215','GFS_RK4_AUDIT_IDS='+','.join(map(str,chosen)))
    s=replace_once(s,'params.dat 20790','params.dat 11058')
    lines=[v for v in s.splitlines() if not ('probe.py validate' in v)]
    s='\n'.join(lines)+'\n';s=replace_once(s,'set +e\n','sha256sum -c inputs.sha256\nset +e\n')
    (root/'run.sh').write_text(s);subprocess.run(['bash','-n',str(root/'run.sh')],check=True)
    deps=[*expected,params,root/'params.dat',root/'probe.json',root/'run.sh']
    (root/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in deps))
    print(json.dumps(cfg,indent=2))
if __name__=='__main__':main()

"""Owned-CSR reassembly versus the existing global compiled-stage oracle.

Stored capture only. This is not an MPI evolution or an executable integration
claim. Canonical fluxes are computed once, including cross-partition faces.
"""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_compiled_midpoint_modules_20261018'
sys.path.insert(0,str(modules))
sys.path.insert(0,str(base/'Codex_full_disk_mpi_bridge_build_v1_20261019'))
from moment_owned_backend import OwnedBackend,partition
from moment_bridge_audit import audit_stage
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
capture=base/'Codex_kinetic_face_capture_20261018';prep=base/'Codex_hydro_cell_backends_20261018';stage=base/'Codex_moment_stage_capture_20261018'
def binary(root,name,dtype,width=None):
    a=np.fromfile(root/name,dtype=dtype);return a.reshape(-1,width) if width else a
data=dict(node_inputs=binary(capture,'inputs.bin','<f8',33),cells=binary(prep,'cells.bin','<f8',19),owners=binary(prep,'owners.bin','<i8',2),
    offsets=binary(prep,'offsets.bin','<i8'),incidence=binary(prep,'incidence.bin','<i8'),frozen=binary(prep,'frozen.bin','<i8'),
    extra=binary(stage,'node_extra.bin','<f8',20),factors=binary(stage,'factors.bin','<f8',4).reshape(-1,2,2),
    gravity_nodes=binary(stage,'gravity_nodes.bin','<f8',7),gravity_offsets=binary(stage,'gravity_offsets.bin','<i8'))
st=np.load(stage/'reference.npz');expect=st['expected'];hy=np.load(prep/'reference.npz')['expected']
reference=dict(raw=expect[:,:18],unmasked_hydro=hy[:,:18],hydro=hy[:,18:36],pressure=hy[:,36:54],source=hy[:,78:96],
    F=hy[:,72:74],H=hy[:,74:78].reshape(-1,2,2),K=expect[:,22:26].reshape(-1,2,2),K_unmasked=expect[:,63:67].reshape(-1,2,2),
    Q=st['Q'],gravity=expect[:,30:48],B=expect[:,18:22].reshape(-1,2,2),G=expect[:,48:52].reshape(-1,2,2),
    Ag=expect[:,52:56].reshape(-1,2,2),cov_A=expect[:,56:60].reshape(-1,2,2))
n=len(data['cells']);nn=len(data['node_inputs']);ni=len(data['incidence']);rows=[]
library=base/'Codex_full_disk_mpi_bridge_build_v1_20261019/liblagmoment_owned.so'

def ledgers(output,nf):
    result=np.zeros(12)
    for k in np.flatnonzero(data['owners'][:,1]<0):
        f=nf[k];result[0]-=f[7];result[1]-=f[8];result[2]-=f[17];result[3]-=f[10]-f[11]
    for i in np.flatnonzero(data['frozen']):
        r=output['hydro'][i];result[4]-=r[7];result[5]-=r[8];result[6]-=r[17];result[7]-=r[10]-r[11]
    for i in np.flatnonzero(~data['frozen'].astype(bool)):
        g=output['stage'][i,30:48];result[8]+=g[7];result[9]+=g[8];result[10]+=g[17];result[11]+=g[10]-g[11]
    return result

for kind,threads in [('cpu',1),('cpu',4),('cuda',1)]:
    backend=OwnedBackend(library,kind,threads);nf=backend.nodes(data['node_inputs']);baseline=backend.stage(data)
    for label,partitions in [('one',[np.arange(n)]),('two_contiguous',np.array_split(np.arange(n),2)),('two_interleaved',[np.arange(0,n,2),np.arange(1,n,2)])]:
        output=dict(stage=np.empty((n,94)),Q=np.empty((ni,16)),hydro=np.empty((n,106)))
        counts=[];cross_faces=0
        membership=np.empty(n,int)
        for rank,owned in enumerate(partitions):membership[owned]=rank
        internal=data['owners'][:,1]>=0
        cross_faces=int(np.sum(membership[data['owners'][internal,0]]!=membership[data['owners'][internal,1]]))
        for owned in partitions:
            local=partition(data,nf,owned);out=backend.owned_stage(n,local)
            output['stage'][owned]=out['stage'];output['hydro'][owned]=out['hydro']
            for li,gi in enumerate(owned):output['Q'][data['offsets'][gi]:data['offsets'][gi+1]]=out['Q'][local['offsets'][li]:local['offsets'][li+1]]
            counts.append(dict(cells=len(owned),receiver_nodes=len(local['node_inputs']),incidences=len(local['incidence'])))
        output['ledgers']=ledgers(output,nf)
        exact={k:bool(np.array_equal(output[k],baseline[k])) for k in ('stage','Q','hydro','ledgers')}
        assert all(exact.values()),(kind,threads,label,exact)
        rows.append(dict(backend=kind,threads=threads,partition=label,owned_counts=counts,cross_partition_nodes=cross_faces,
            canonical_nodes_computed_once=nn,exact_global_reassembly=exact))
        print(json.dumps(rows[-1]),flush=True)
        if label=='two_contiguous':
            # One independent floating point oracle uses the existing bounds.
            audit=audit_stage(data,reference,{kind+str(threads):output})
            (args.output/f'oracle_{kind}{threads}.json').write_text(json.dumps(audit,indent=2))
    local=partition(data,nf,np.arange(n//2));failures={}
    for name in ('stale_epoch','wrong_owner','missing_incidence','bad_spatial_moment'):
        bad={k:v.copy() for k,v in local.items()};expected_epoch=1
        if name=='stale_epoch':expected_epoch=2
        if name=='wrong_owner':bad['owned_registry_index'][0]=bad['owned_registry_index'][1]
        if name=='missing_incidence':bad['incidence'][0]=-bad['incidence'][0]
        if name=='bad_spatial_moment':bad['cells'][0,3]=-1.
        status,out=backend.owned_stage(n,bad,expected_epoch=expected_epoch,sentinel=-947.)
        assert status!=0 and all(np.all(v==-947.) for v in out.values()),(kind,name,status)
        failures[name]=dict(status=status,caller_outputs_unchanged=True)
    rows.append(dict(backend=kind,threads=threads,transactional_failures=failures))
result=dict(status='passed_owned_moment_bridge_component',cells=n,nodes=nn,partitions=rows,
    library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),seconds=time.perf_counter()-started,
    physical_providers=0,evolution_steps=0,MPI_evolution_verified=False,
    limits='Captured-stage owned ABI only. Canonical flux and physical frozen flags preserved, no rank walls or local massroots. Native executable/runtime/service and actual MPI transactions still pending.')
(args.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(status=result['status'],seconds=result['seconds'])),flush=True)

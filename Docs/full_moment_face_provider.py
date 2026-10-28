"""Original global face/root/audit path with the GPU stage owned by native MPI.

The service yields typed data to the native controller at backend.stage, then
continues the unchanged wall/gravity/conservation audits on native outputs.
No service integrator, MPI, CUDA initialization, or state acceptance.
"""
import ctypes as C,json,os,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from full_moment_provider_wire import *
from full_moment_state_io import load as load_state

def stage_module():
    provider_dir = Path(__file__).resolve().parent
    fallback_dir = '/gpfs/kjhan/LagForce/Codex_full_disk_fast_modules_v1_20261019'
    if os.environ.get('LAGFORCE_NATIVE_FACE_STAGE')=='1':
        local_stage = provider_dir / 'native_full_disk_stage.py'
        if local_stage.is_file():
            while fallback_dir in sys.path:
                sys.path.remove(fallback_dir)
            if str(provider_dir) in sys.path:
                sys.path.remove(str(provider_dir))
            sys.path.insert(0, str(provider_dir))
            sys.path.append(fallback_dir)
        else:
            sys.path.insert(0, fallback_dir)
        import native_full_disk_stage
        return native_full_disk_stage
    sys.path.insert(0, fallback_dir)
    import full_disk_stage_memory
    return full_disk_stage_memory

class GravityShape:
    """Counts only; the real quadrature values remain in owned rank buffers."""
    def __init__(self,count):self.count=count;self.nbytes=count*7*8
    def __len__(self):return self.count

def evaluate(service,request,response,read_exact,write_exact,write_json):
    stage_api=stage_module()
    ctx=service.ctx;root=service.root
    assert request.rank==0 and VELOCITY in service.completed
    fits=json.loads(Path(bytes(request.fits_path).decode()).read_text())
    rows=json.loads(Path(bytes(request.aux_path).decode()).read_text())
    initial=load_state(ctx['state'])
    with np.load(Path(ctx['geometry'])/'mesh.npz') as f:mesh={k:f[k] for k in f.files}
    native=json.loads((Path(ctx['geometry'])/'full_exact_supports.json').read_text())
    ids=initial['ids'];n=len(ids)
    assert [r['id'] for r in rows]==ids.tolist() and [r['id'] for r in fits]==ids.tolist()
    assert all(r['status']=='passed' for r in rows+fits)
    prepared={k:initial[k] for k in ('ids','state','factors','mass','Psi')}
    prepared.update({key:np.array([r[key] for r in rows]) for key in ('generator','gravity','gravity_absolute_sum','moment_defect')})
    prepared.update({key:np.array([r['metadata'][key] for r in rows]) for key in rows[0]['metadata']})
    prepared['gravity_offsets']=np.r_[0,np.cumsum([r['nodes24'] for r in fits])].astype(np.int64)
    ng=int(prepared['gravity_offsets'][-1]);native_root=root/'native_stage';native_root.mkdir(exist_ok=False)
    out_root=native_root/'outputs';out_root.mkdir()
    class NativeBackend:
        kind=1
        def stage(self,data):
            assert len(data['cells'])==n and len(data['gravity_nodes'])==ng
            header=np.array([request.epoch,request.run_token,n,len(data['node_inputs']),len(data['incidence']),ng],dtype='<u8')
            arrays={k:v for k,v in data.items() if k!='gravity_nodes'}
            arrays['header']=header
            for key,value in arrays.items():
                value=np.ascontiguousarray(value)
                assert value.dtype in (np.dtype('<f8'),np.dtype('<i8'),np.dtype('<u8')),(key,value.dtype)
                with (native_root/(key+'.bin')).open('xb') as f:value.tofile(f)
            event=Response.from_buffer_copy(bytes(response));event.status=NATIVE_STAGE_REQUIRED
            event.output_path=str(native_root).encode();event.gravity_nodes=ng
            event.message=b'Canonical face/root ready; native MPI owned GPU stage required'
            write_exact(service.protocol,bytes(event))
            raw=read_exact(0,C.sizeof(Request));assert raw is not None,'Native controller closed before GPU result'
            resume=Request.from_buffer_copy(raw)
            assert resume.magic==MAGIC and resume.version==ABI and resume.command==RESUME_STAGE
            assert resume.request_id==service.sequence+1 and resume.epoch==service.epoch and resume.run_token==service.run_token
            assert (resume.rank,resume.nranks,resume.ncells,bytes(resume.physics_sha))==service.signature
            assert resume.ld_bytes==np.dtype(np.longdouble).itemsize and resume.ld_mantissa==np.finfo(np.longdouble).nmant+1
            service.sequence=resume.request_id
            assert Path(bytes(resume.output_dir).decode()).resolve()==out_root
            result={}
            for key,shape in [('stage',(n,94)),('Q',(len(data['incidence']),16)),('hydro',(n,106)),('ledgers',(12,))]:
                path=out_root/(key+'.bin');assert path.stat().st_size==int(np.prod(shape))*8
                result[key]=np.memmap(path,mode='r',dtype=np.float64,shape=shape)
            return result
    def event(row):
        print(json.dumps(dict(full_moment_face_event=row)),flush=True)
    result=stage_api.evaluate(initial,native,mesh,fits,prepared,GravityShape(ng),ctx['stage_key'],NativeBackend(),root/'actual_stage_audit',callback=event)
    rates=root/'rates';rates.mkdir(exist_ok=False)
    for key in ('generator','weight_rates','Psi_rates','ledgers'):
        value=np.ascontiguousarray(result.saved[key]);assert value.dtype==np.dtype('<f8')
        with (rates/(key+'.bin')).open('xb') as f:value.tofile(f)
    # Compact rates plus incoming columns are node-local, consumed by the
    # native integrator. NPZ here is serialization only, not state evolution.
    with (root/'stage_rates.npz').open('xb') as f:np.savez(f,**result.saved)
    report=dict(result.report,native_MPI_GPU_stage=True,service_MPI_calls=0,service_GPU_calls=0)
    if hasattr(result,'native_provenance'):report['native_preparation']=result.native_provenance
    if os.environ.get('LAGFORCE_COUPLED_JOURNAL')=='1':
        journal=root/'actual_stage_audit/cell_transport_journal.json'
        record=json.loads(journal.read_text())
        assert record['stage_key']==ctx['stage_key'] and record['time']==float(initial['time'])
        report['cell_transport_journal']=dict(report_path=str(journal),**record)
    write_json(root/'actual_stage_report.json',report)
    response.request_id=service.sequence;response.output_path=str(root/'actual_stage_report.json').encode()
    response.gravity_nodes=ng;response.message=b'Native GPU outputs passed original global physical stage audits'
    service.completed.add(FACE_STAGE)
    return response

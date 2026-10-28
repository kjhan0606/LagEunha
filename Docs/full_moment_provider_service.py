"""Persistent typed provider for the native FullMomentKepler MPI controller.

CPU cell preparation only. No MPI calls, CUDA, integrator, restart, retry or
state acceptance. The native controller owns global assembly and commitment.
"""
import argparse,ctypes as C,hashlib,json,multiprocessing as mp,os,sys,time,traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_PROC_BIND='false')
import numpy as np
from full_moment_provider_wire import *
from full_moment_provider_workers import initialize,task
from full_moment_state_io import load as load_state,dump as dump_state,digest as state_digest,HEADER as STATE_HEADER

def text_path(field):return bytes(field).decode('utf8')
def read_exact(fd,n):
    chunks=[];left=n
    while left:
        data=os.read(fd,left)
        if not data:
            if left==n:return None
            raise EOFError('Truncated provider request')
        chunks.append(data);left-=len(data)
    return b''.join(chunks)
def write_exact(fd,data):
    view=memoryview(data)
    while view:
        n=os.write(fd,view)
        if n<=0:raise BrokenPipeError('Provider response pipe closed')
        view=view[n:]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,data):
    with Path(path).open('x') as f:json.dump(data,f,separators=(',',':'))

class Service:
    def __init__(self,workers):
        assert 1<=workers<=len(os.sched_getaffinity(0))
        assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1'
        assert os.environ['OMP_PROC_BIND']=='false'
        # Resolve all global-face dependencies before any expensive owned fit.
        from full_moment_face_provider import stage_module
        stage_module()
        self.pool=ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize)
        self.workers=workers;self.run_token=None;self.sequence=0;self.epoch=0
        self.signature=None;self.ctx=None;self.root=None;self.live=[]
        self.completed=set();self.stopped=False;self.pending_geometry=None
    def handle(self,r):
        assert r.magic==MAGIC and r.version==ABI
        assert r.ld_bytes==np.dtype(np.longdouble).itemsize and r.ld_mantissa==np.finfo(np.longdouble).nmant+1
        assert r.run_token and r.request_id==self.sequence+1
        if self.run_token is None:self.run_token=r.run_token
        assert self.run_token==r.run_token
        self.sequence=r.request_id
        response=Response(magic=MAGIC,epoch=r.epoch,run_token=r.run_token,
            request_id=r.request_id,version=ABI,ncells=r.ncells,
            ld_bytes=r.ld_bytes,ld_mantissa=r.ld_mantissa)
        if r.command==STOP:self.stopped=True;return response
        signature=(r.rank,r.nranks,r.ncells,bytes(r.physics_sha))
        if self.signature is None:
            assert r.command in (CONFIGURE,VERIFY_RESTART,VERIFY_FORK) and r.nranks>0 and r.rank<r.nranks and r.ncells>=r.nranks
            assert len(bytes(r.physics_sha))==64
            int(bytes(r.physics_sha),16)
            self.signature=signature
        assert self.signature==signature,'Rank/physics registry changed'
        n=r.ncells;lo=n*r.rank//r.nranks;hi=n*(r.rank+1)//r.nranks
        response.owned_start=lo;response.owned_count=hi-lo
        if r.command==VERIFY_FORK:
            assert r.rank==0 and self.epoch==0 and self.ctx is None and not self.live
            assert not hasattr(self,'fork_verified') and not hasattr(self,'restart_verified')
            import struct
            from full_moment_coupled_fork import verify_serialize
            output=Path(text_path(r.output_dir)).resolve()
            assert output==Path(f'/dev/shm/eunha_run_{r.run_token}/fork_verified')
            header,record=verify_serialize(text_path(r.state_path),text_path(r.bounds_path),
                text_path(r.fits_path),text_path(r.aux_path),struct.unpack('<d',bytes(r.reserved)[:8])[0],output)
            assert record['physics_sha256']==bytes(r.physics_sha).decode() and record['cells']==n
            self.fork_verified=True;response.output_path=str(output).encode()
            response.message=b'Explicit-to-coupled fork verified and serialized; native token import required'
            return response
        if r.command==MIDPOINT_OBSERVE and r.outer==1:
            assert r.rank==0 and self.ctx is None and not self.live and self.pending_geometry is None and r.epoch>self.epoch
            source=Path(text_path(r.state_path))
            header=STATE_HEADER.unpack((source/'header.bin').read_bytes())
            assert header[2]==n and header[3]==r.epoch and header[4]==r.run_token
            from full_moment_coupled_observation import summary,agreement
            physics_file=Path(text_path(r.aux_path));physics=json.loads(physics_file.read_text())['physics']
            assert hashlib.sha256(json.dumps(physics,sort_keys=True,separators=(',',':')).encode()).hexdigest()==bytes(r.physics_sha).decode()
            proof=summary(source,text_path(r.geometry_path),text_path(r.fits_path),physics['gamma'])
            with np.load(Path(text_path(r.bounds_path))/'summary.npz') as f:provider={k:f[k] for k in f.files}
            output=agreement(provider,proof,text_path(r.output_dir))
            response.output_path=str(output).encode();response.message=b'Actual provider and accepted arithmetic midpoint density/pressure agree'
            return response
        if r.command==VERIFY_RESTART:
            assert r.rank==0 and self.epoch==0 and self.ctx is None and not self.live
            assert not hasattr(self,'restart_verified'),'Restart verification is a single bootstrap transaction'
            import struct
            from full_moment_restart_provenance import verify
            requested_h=struct.unpack('<d',bytes(r.reserved)[:8])[0]
            header,record=verify(text_path(r.state_path),text_path(r.bounds_path),
                text_path(r.fits_path),text_path(r.aux_path),requested_h)
            assert record['physics_sha256']==bytes(r.physics_sha).decode() and record['cells']==n
            output=Path(text_path(r.output_dir)).resolve()
            assert output==Path(f'/dev/shm/eunha_run_{r.run_token}/restart_verified')
            output.mkdir(exist_ok=False)
            with (output/'expected_header.bin').open('xb') as f:f.write(header)
            write_json(output/'receipt.json',dict(status='verified_restart_provenance_no_state_update',
                source=str(Path(text_path(r.state_path)).resolve()),
                manifest_sha256=sha(text_path(r.bounds_path)),time=record['time'],
                old_run_token=record['original_run_token'],new_run_token=int(r.run_token),
                old_epoch=record['original_epoch'],new_epoch=1,physics_sha256=record['physics_sha256'],
                requested_h=requested_h,all_source_files_verified=True,numerical_updates=0))
            self.restart_verified=True;response.output_path=str(output).encode()
            response.message=b'Checkpoint model, baseline, payload and audit hashes verified; native import required'
            return response
        if r.command==DISPOSE_GEOMETRY:
            # Native collective RELEASE acknowledgement must precede this
            # separate rank0-only command. No rank is still hashing this mesh.
            assert r.rank==0 and r.epoch==self.epoch and self.ctx is None and not self.live
            if self.pending_geometry is not None:
                import shutil
                geometry=Path(self.pending_geometry)
                assert geometry==Path(text_path(r.geometry_path)).resolve()
                prefix=Path(f'/dev/shm/eunha_run_{r.run_token}/geometry')
                assert geometry.name=='geometry' and geometry.parent.parent==prefix
                shutil.rmtree(geometry.parent);self.pending_geometry=None
            response.message=b'Shared geometry disposed after collective release acknowledgement'
            return response
        if r.command==GEOMETRY:
            assert r.rank==0 and r.epoch>self.epoch and not self.live
            from full_moment_global_audit import geometry
            path=geometry(self,r)
            response.output_path=str(path).encode();response.message=b'Global true-circle geometry passed original gates'
            return response
        if r.command==CONFIGURE:
            assert r.epoch>self.epoch and not self.live,'Unreleased or stale epoch'
            assert self.pending_geometry is None,'Shared geometry disposal acknowledgement missing'
            state=Path(text_path(r.state_path)).resolve();geometry=Path(text_path(r.geometry_path)).resolve()
            self.root=Path(text_path(r.output_dir)).resolve();self.root.mkdir(exist_ok=False)
            assert str(self.root).startswith('/dev/shm/'),'Provider bulk buffers must be node-local shared memory'
            initial=load_state(state)
            assert len(initial['ids'])==n and initial['Psi'].dtype==np.dtype(np.longdouble)
            assert np.all(np.diff(initial['ids'])>0)
            if state.is_dir():
                header=STATE_HEADER.unpack((state/'header.bin').read_bytes())
                assert header[3]==r.epoch and header[4]==r.run_token,'Foreign state epoch/token'
            physics_file=Path(text_path(r.aux_path)).resolve()
            physics=json.loads(physics_file.read_text())['physics']
            assert hashlib.sha256(json.dumps(physics,sort_keys=True,separators=(',',':')).encode()).hexdigest()==bytes(r.physics_sha).decode(), 'Physics manifest SHA mismatch'
            self.ctx=dict(state=str(state),geometry=str(geometry),ncells=n,lo=lo,hi=hi,
                epoch=int(r.epoch),run_token=int(r.run_token),rank=int(r.rank),nranks=int(r.nranks),
                stage_key=state_digest(state)+':'+sha(geometry/'full_exact_supports.json'),
                physics_sha=bytes(r.physics_sha).decode(),q48=bool(r.outer),bounds=None,fits=None,
                state_sha=state_digest(state),geometry_sha=sha(geometry/'full_exact_supports.json'),mesh_sha=sha(geometry/'mesh.npz'),physics_file=str(physics_file))
            self.epoch=r.epoch;self.completed=set();write_json(self.root/'context.json',self.ctx)
            response.message=b'Configured global immutable state; no evolution performed'
            return response
        assert self.ctx is not None and r.epoch==self.epoch,'Stale or unconfigured epoch'
        assert state_digest(self.ctx['state'])==self.ctx['state_sha'],'State modified within epoch'
        assert sha(Path(self.ctx['geometry'])/'mesh.npz')==self.ctx['mesh_sha'],'Mesh modified within epoch'
        if r.command==MIDPOINT_OBSERVE:
            assert r.rank==0 and r.outer==0 and FACE_STAGE in self.completed and MIDPOINT_OBSERVE not in self.completed
            from full_moment_coupled_observation import capture_provider
            output=capture_provider(self,r);self.completed.add(MIDPOINT_OBSERVE)
            response.output_path=str(output).encode();return response
        if r.command==COUPLED_JOURNAL:
            assert r.rank==0 and DENSITY in self.completed and COUPLED_JOURNAL not in self.completed
            from full_moment_coupled_observation import final_journal
            output=final_journal(self,r);self.completed.add(COUPLED_JOURNAL)
            response.output_path=str(output).encode();return response
        if r.command==COUPLED_METADATA:
            assert r.rank==0 and DENSITY in self.completed and COUPLED_METADATA not in self.completed
            from full_moment_coupled_metadata import write_metadata
            output=self.root/'coupled_metadata'
            saved=load_state(self.ctx['state'])
            write_metadata(text_path(r.fits_path),self.ctx['geometry'],saved['ids'],output,self.ctx.get('physics_file'))
            response.output_path=str(output).encode()
            response.message=b'Actual nonlinear residual scales exported; no physics evaluation'
            self.completed.add(COUPLED_METADATA);return response
        if r.command==EXPORT_STATE:
            assert EXPORT_STATE not in self.completed
            path=dump_state(self.root/'state_flat',load_state(self.ctx['state']),r.epoch,r.run_token)
            response.output_path=str(path).encode();response.message=b'Exact typed native state serialization'
            self.completed.add(EXPORT_STATE);return response
        if r.command==FACTOR_AUDIT:
            from full_moment_factor_audit import evaluate
            path=evaluate(self,r)
            response.output_path=str(path).encode();response.message=b'Original factor conditioning and positivity gates passed'
            return response
        if r.command==ADMISSIBILITY:
            assert r.rank==0 and DENSITY in self.completed
            from full_moment_global_audit import admissibility
            path=admissibility(self,r)
            response.output_path=str(path).encode();response.message=b'Tentative state passed requested physical gates; native commit required'
            return response
        if r.command==RELEASE:
            # Explicit native acknowledgement of consumption; no recursive cleanup.
            for path in self.live:Path(path).unlink()
            self.live=[];response.message=b'Released owned node buffers'
            if r.outer:
                # Explicit acknowledgement from the native transaction owner.
                # Failed and historical runs never enter this branch.
                import shutil
                assert self.root.parent==Path('/dev/shm') and self.root.name.startswith('eunha_fmp_')
                geometry=Path(self.ctx['geometry']);prefix=Path(f'/dev/shm/eunha_run_{r.run_token}/geometry')
                if r.rank==0 and geometry.name=='geometry' and geometry.parent.parent==prefix:
                    assert self.pending_geometry is None
                    self.pending_geometry=str(geometry)
                shutil.rmtree(self.root);self.root=None;self.ctx=None;self.completed=set()
            return response
        if r.command==FACE_STAGE:
            assert FACE_STAGE not in self.completed
            assert state_digest(self.ctx['state'])==self.ctx['state_sha']
            assert sha(Path(self.ctx['geometry'])/'full_exact_supports.json')==self.ctx['geometry_sha']
            from full_moment_face_provider import evaluate
            return evaluate(self,r,response,read_exact,write_exact,write_json)
        assert r.command in (BOUNDS,DENSITY,VELOCITY) and r.command not in self.completed
        assert state_digest(self.ctx['state'])==self.ctx['state_sha'],'State modified within epoch'
        assert sha(Path(self.ctx['geometry'])/'full_exact_supports.json')==self.ctx['geometry_sha']
        ctx=dict(self.ctx)
        if r.command>=DENSITY:
            ctx['bounds']=str(Path(text_path(r.bounds_path)).resolve())
            assert Path(ctx['bounds']).is_file()
        gravity=None
        if r.command==VELOCITY:
            ctx['fits']=str(Path(text_path(r.fits_path)).resolve())
            fits=json.loads(Path(ctx['fits']).read_text())
            assert len(fits)==n and all(row['status']=='passed' for row in fits)
            ng=sum(row['nodes24'] for row in fits[lo:hi]);assert ng>0
            gravity=self.root/'gravity.bin';fd=os.open(gravity,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
            try:os.ftruncate(fd,ng*7*8)
            finally:os.close(fd)
            self.live.append(str(gravity));ctx['gravity']=str(gravity);response.gravity_nodes=ng
        descriptor=self.root/f'context_{r.command}.json';write_json(descriptor,ctx)
        started=time.perf_counter();rows=[]
        for i,row in zip(range(lo,hi),self.pool.map(task,((str(descriptor),r.command,i) for i in range(lo,hi)),chunksize=4)):
            if row['status']!='passed':
                write_json(self.root/f'failure_cell_{i}.json',row);raise RuntimeError(f'Cell {i} preparation failed: {row}')
            rows.append(row)
        path=self.root/f'rows_{r.command}.json';write_json(path,rows)
        result=dict(status='passed_owned_cpu_preparation',command=int(r.command),
            epoch=int(r.epoch),rank=int(r.rank),ranks=int(r.nranks),lo=lo,hi=hi,
            workers=self.workers,seconds=time.perf_counter()-started,rows=str(path),
            cells=hi-lo,gravity_path=str(gravity) if gravity else None,
            gravity_nodes=int(response.gravity_nodes),
            ownership='Contiguous rows in global stable-ID registry; all stencils global',
            provider_calls=0,evolution_steps=0)
        manifest=self.root/f'result_{r.command}.json';write_json(manifest,result)
        response.output_path=str(manifest).encode();response.payload_bytes=path.stat().st_size
        response.token=r.request_id;response.message=b'Owned preparation complete; native acceptance required'
        self.completed.add(r.command)
        return response

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workers',type=int,required=True);args=p.parse_args()
    # Reserve protocol FD before spawned workers inherit stdout. All incidental
    # Python/native output goes to the per-rank diagnostic log, never the pipe.
    protocol=os.dup(1);os.set_inheritable(protocol,False);os.dup2(2,1)
    service=Service(args.workers);service.protocol=protocol
    try:
        while True:
            data=read_exact(0,C.sizeof(Request))
            if data is None:break
            r=Request.from_buffer_copy(data)
            try:response=service.handle(r)
            except Exception as exc:
                witness=dict(status='failed_provider_request',error=repr(exc),traceback=traceback.format_exc(),
                    epoch=int(r.epoch),request_id=int(r.request_id),command=int(r.command))
                print(json.dumps(witness),file=sys.stderr,flush=True)
                if service.root is not None:
                    write_json(service.root/f'failure_request_{r.request_id}.json',witness)
                response=Response(magic=MAGIC,version=ABI,epoch=r.epoch,run_token=r.run_token,
                    request_id=service.sequence,status=1,ncells=r.ncells,
                    ld_bytes=np.dtype(np.longdouble).itemsize,ld_mantissa=np.finfo(np.longdouble).nmant+1)
                response.message=str(exc).encode()[:255]
                write_exact(protocol,bytes(response));return 1
            write_exact(protocol,bytes(response))
            if service.stopped:break
    finally:service.pool.shutdown(wait=True,cancel_futures=True);os.close(protocol)
    return 0

if __name__=='__main__':sys.exit(main())

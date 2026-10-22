"""Private positive-support workers around the native full preprocessing.

Workers receive moment matrices, return nonnegative weights, and never evolve
positions, velocities, mass or heat. Native C++/CUDA computes all physical data.
"""
import json
import multiprocessing as mp
import os
import pathlib
import socket
import struct
import subprocess
import sys
import tempfile
import time

def receive(sock,count):
    data=bytearray(count); view=memoryview(data)
    while view:
        n=sock.recv_into(view)
        if not n: raise EOFError('native moment stream closed')
        view=view[n:]
    return data

def worker(path):
    os.environ['OPENBLAS_NUM_THREADS']='4'
    os.environ['MKL_NUM_THREADS']='4'
    os.environ['OMP_NUM_THREADS']='4'
    import numpy as np
    from scipy.optimize import nnls
    listener=socket.socket(socket.AF_UNIX,socket.SOCK_STREAM)
    listener.bind(path);os.chmod(path,0o600);listener.listen(1)
    connection,_=listener.accept()
    try:
        while True:
            header=receive(connection,40)
            magic,rows,n,source,_=struct.unpack('<5Q',header)
            if magic!=0x4d46504f53525131 or rows>400 or n>1048576 or not rows or not n:
                raise RuntimeError('invalid bounded native moment packet')
            B=np.frombuffer(receive(connection,8*rows*n),dtype='<f8').reshape(rows,n)
            b=np.frombuffer(receive(connection,8*rows),dtype='<f8')
            start=time.monotonic_ns();error=0
            try:
                proposal,_=nnls(B,b,maxiter=3*n)
                if not np.isfinite(proposal).all() or np.any(proposal<0):raise RuntimeError('nonfinite/nonpositive candidate')
                active=np.flatnonzero(proposal>0).astype('<u4')
                values=np.asarray(proposal[active],dtype='<f8')
            except Exception:
                active=np.empty(0,dtype='<u4');values=np.empty(0,dtype='<f8');error=1
            connection.sendall(struct.pack('<5Q',0x4d46504f53525031,source,len(active),error,time.monotonic_ns()-start))
            connection.sendall(active.tobytes());connection.sendall(values.tobytes())
    except EOFError:
        pass
    finally:
        connection.close();listener.close()

def main():
    dest=pathlib.Path(sys.argv[1]).resolve()
    root=dest.parent
    assert dest.name=='material_full_positive_rule_v1_ada_mpi4'
    (dest/'supervisor_pid.txt').write_text(str(os.getpid())+'\n')
    processes=[];context=mp.get_context('spawn')
    with tempfile.TemporaryDirectory(prefix='lf_positive_',dir='/tmp') as private:
        os.chmod(private,0o700)
        try:
            for i in range(16):
                p=context.Process(target=worker,args=(str(pathlib.Path(private)/f'w{i}.sock'),));p.start();processes.append(p)
            (dest/'worker_pids.json').write_text(json.dumps({'workers':16,'BLAS_threads_per_worker':4,'pids':[p.pid for p in processes],'socket_directory':private},indent=2)+'\n')
            deadline=time.monotonic()+120
            while not all((pathlib.Path(private)/f'w{i}.sock').exists() for i in range(16)):
                if not all(p.is_alive() for p in processes) or time.monotonic()>deadline:raise RuntimeError('worker launch failed')
                time.sleep(.1)
            command=['mpirun','-np','4',str(dest/'preprocess.exe'),
                     str(root/'material_full_initial_v2_h200_mpi4/run'),str(dest/'immutable_payload.sha256'),
                     str(root/'material_full_cached_continuation_v1_ada_mpi4/run/continuous_step4/accepted_r'),
                     str(root/'material_full_resident_rhs_v1_ada_mpi4/run/next_RK_stage1/prepared_k1_r'),
                     private,str(dest/'run'),str(root/'material_positive_support_v1'),'86400']
            (dest/'native_command.json').write_text(json.dumps(command,indent=2)+'\n')
            with (dest/'contract.json').open('x') as out,(dest/'stderr.log').open('x') as error:
                p=subprocess.Popen(command,stdout=out,stderr=error)
                (dest/'native_pid.txt').write_text(str(p.pid)+'\n')
                (dest/'status.txt').write_text('RUNNING_FULL_NATIVE_POSITIVE_PREPROCESSING_NO_EVOLUTION\n')
                code=p.wait()
            (dest/'native_exit_code.txt').write_text(str(code)+'\n')
            if code==0:
                result=json.loads((dest/'contract.json').read_text())
                status='FULL_CANDIDATE_READY_FOR_RK_ORACLE' if result['candidate_ready_for_actual_RK_oracle'] else 'PREPROCESS_PREFIX_OR_COST_GATE_STOP_NO_OPERATOR_COMMIT'
            else: status=f'FULL_NATIVE_PREPROCESS_FAILED_EXIT_{code}'
            (dest/'status.txt').write_text(status+'\n')
            return code
        finally:
            for p in processes:
                if p.is_alive():p.terminate()
            for p in processes:p.join(10)

if __name__=='__main__':sys.exit(main())

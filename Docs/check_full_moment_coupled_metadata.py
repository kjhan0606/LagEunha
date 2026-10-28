"""Actual saved-fit metadata/wire/physical-residual ABI check, zero RHS."""
import argparse
import ctypes as C
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from check_full_moment_coupled_residual import Potential,Residual
from full_moment_coupled_metadata import POTENTIAL_DTYPE


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    own=Path(__file__).resolve().parent;b=Path('/gpfs/kjhan/LagForce')
    sys.path.insert(0,str(b/'Codex_full_disk_native_evolution_build_v10_20261019'));sys.path.insert(0,str(own))
    from full_moment_provider_service import Service
    from full_moment_provider_wire import Request,MAGIC,ABI,DENSITY,COUPLED_METADATA
    cp=b/'Codex_full_disk_continuous_two_dt_to_t1_20261019/fine/checkpoints/fine_0010'
    files=[cp/'accepted.npz',cp/'geometry/mesh.npz',cp/'density/cells.jsonl']
    source={str(path):sha(path) for path in files}
    rows=[json.loads(line) for line in files[2].read_text().splitlines()]
    fit=args.output/'fit_rows.json';fit.write_text(json.dumps(rows,separators=(',',':')))
    with np.load(files[0]) as f:state=f['state'];mass=f['mass'];psi=f['Psi'];ids=f['ids']
    root=args.output/'configured_service';root.mkdir()
    service=object.__new__(Service)
    physics=b'0'*64
    service.run_token=123456789;service.sequence=0;service.epoch=7
    service.signature=(0,1,len(ids),physics);service.root=root;service.completed={DENSITY}
    service.ctx=dict(state=str(files[0]),geometry=str(cp/'geometry'),state_sha=source[str(files[0])],mesh_sha=source[str(files[1])])
    request=Request(magic=MAGIC,version=ABI,run_token=service.run_token,request_id=1,
        epoch=service.epoch,command=COUPLED_METADATA,rank=0,nranks=1,ncells=len(ids),
        ld_bytes=16,ld_mantissa=64,physics_sha=physics,fits_path=str(fit).encode())
    response=service.handle(request);output=Path(bytes(response.output_path).decode())
    report=json.loads((output/'report.json').read_text())
    assert report['status']=='exported_actual_coupled_residual_metadata'
    assert report['physical_providers']==report['accepted_time_steps']==0
    assert POTENTIAL_DTYPE.itemsize==C.sizeof(Potential)
    assert all(POTENTIAL_DTYPE.fields[name][1]==getattr(Potential,name).offset for name,_ in Potential._fields_)
    norms=np.fromfile(output/'potential_norms.bin',dtype=POTENTIAL_DTYPE)
    area=np.fromfile(output/'areas.bin',dtype=np.float64)
    assert np.array_equal(np.fromfile(output/'ids.bin',dtype=np.int64),ids)
    assert len(norms)==len(ids) and len(area)==len(ids)
    retained_precision=bool(np.any(norms['sigma']!=norms['sigma'].astype(np.float64).astype(np.longdouble)))
    assert retained_precision
    lib=C.CDLL(str(args.library));lib.full_moment_coupled_residual.argtypes=[C.c_int,C.c_double]+[C.c_void_p]*10
    lib.full_moment_coupled_residual.restype=C.c_int
    lib.full_moment_coupled_converged.argtypes=[C.POINTER(Residual)];lib.full_moment_coupled_converged.restype=C.c_int
    ptr=lambda a:a.ctypes.data_as(C.c_void_p)
    result=Residual();failed=C.c_int(-1)
    code=lib.full_moment_coupled_residual(len(ids),float(np.float32(1.666667)),ptr(state),ptr(mass),ptr(area),ptr(state),ptr(state),ptr(psi),ptr(psi),ptr(norms),C.byref(result),C.byref(failed))
    assert code==0 and max(result.maximum)==0 and lib.full_moment_coupled_converged(C.byref(result))
    request.request_id=2
    try:service.handle(request)
    except AssertionError:duplicate_rejected=True
    else:raise AssertionError('Duplicate metadata export was accepted')
    assert all(sha(Path(path))==digest for path,digest in source.items())
    summary=dict(status='passed_actual_coupled_metadata_wire_fixture',cells=len(ids),
        source_sha256=source,configured_context_fixture_only=True,physical_providers=0,accepted_time_steps=0,
        binary80_precision_beyond_binary64_retained=retained_precision,struct_bytes=POTENTIAL_DTYPE.itemsize,
        original_sources_unchanged=True,actual_fit_identity_residual_zero=True,duplicate_export_rejected=duplicate_rejected,
        exported_metadata=report,library_sha256=sha(args.library),
        limits='Unit-configured provider handler and native residual ABI only; no actual MPI/nonlinear time step.')
    (args.output/'report.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':main()

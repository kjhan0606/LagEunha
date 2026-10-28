"""Export actual fit scales/areas for the native nonlinear residual evaluator.

No reconstruction, physics evaluation or state update is performed here.
Binary80 values are parsed directly into NumPy longdouble, never via float.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

POTENTIAL_DTYPE=np.dtype([('origin',np.longdouble,(2,)),('scale',np.longdouble),
    ('reference',np.longdouble,(6,)),('projection',np.longdouble,(6,)),
    ('sigma',np.longdouble),('rank_error',np.longdouble),('dependent',np.int32)],align=True)


def potential_norms(rows):
    result=np.zeros(len(rows),dtype=POTENTIAL_DTYPE)
    for i,row in enumerate(rows):
        assert row['status']=='passed'
        p=row['parameter'];q=result[i]
        q['origin']=np.asarray(p['origin'],dtype=np.longdouble);q['scale']=np.longdouble(p['scale'])
        q['reference']=np.asarray(p['reference'],dtype=np.longdouble)
        q['projection']=q['reference']+np.asarray(p['projection_remainder'],dtype=np.longdouble)
        q['sigma']=np.longdouble(p['sigma']);q['rank_error']=np.longdouble(p['rank_error'])
        dependent=bool(p.get('polynomial_dependent',len(p['theta'])==5))
        assert len(p['theta'])==(5 if dependent else 6)
        q['dependent']=int(dependent)
        assert q['scale']>0 and all(np.isfinite(q[name]).all() for name in POTENTIAL_DTYPE.names if name!='dependent')
        if not dependent:assert q['sigma']>q['rank_error']>=0
    return result


def write_metadata(fits_path,geometry_path,ids,output,physics_path=None):
    fits_path=Path(fits_path);geometry_path=Path(geometry_path);output=Path(output)
    rows=json.loads(fits_path.read_text())
    assert [row['id'] for row in rows]==np.asarray(ids).tolist()
    mesh_path=geometry_path/'mesh.npz'
    with np.load(mesh_path) as f:
        assert np.array_equal(f['ids'],ids)
        area=np.asarray(f['moments'][:,0],dtype=np.float64)
    norms=potential_norms(rows)
    assert sys.byteorder=='little' and np.dtype(np.longdouble).itemsize==16 and np.finfo(np.longdouble).nmant==63
    assert POTENTIAL_DTYPE.itemsize==288
    assert np.all(np.isfinite(area)) and np.all(area>0)
    output.mkdir(exist_ok=False)
    fit_area=np.array([row['parameter']['area'] for row in rows],dtype=np.float64)
    payload=[('areas',area),('fit_areas',fit_area),('potential_norms',norms),('ids',np.asarray(ids,dtype=np.int64))]
    if physics_path is not None:
        physics=json.loads(Path(physics_path).read_text())['physics']
        # The production primitive closure uses this exact binary64 value.
        # Do not substitute a higher precision gamma from the JSON descriptor.
        gamma=float(np.float32(1.666667))
        assert np.isfinite(gamma) and gamma>1 and physics['gamma']>1
        payload.append(('gamma',np.array([gamma],dtype=np.float64)))
    for name,value in payload:
        with (output/(name+'.bin')).open('xb') as f:value.tofile(f)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    report=dict(status='exported_actual_coupled_residual_metadata',cells=len(ids),
        source_fit_sha256=sha(fits_path),source_mesh_sha256=sha(mesh_path),
        potential_struct_bytes=POTENTIAL_DTYPE.itemsize,
        potential_field_offsets={name:POTENTIAL_DTYPE.fields[name][1] for name in POTENTIAL_DTYPE.names},
        binary80_mantissa_bits=64,dependent_cells=int(norms['dependent'].sum()),
        maximum_fit_geometric_area_difference=float(abs(area-fit_area).max(initial=0.)),
        payload_sha256={name+'.bin':sha(output/(name+'.bin')) for name,_ in payload},
        physical_providers=0,accepted_time_steps=0)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return output

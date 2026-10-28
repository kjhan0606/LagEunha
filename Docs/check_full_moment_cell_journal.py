"""Observer conservation/sign/decomposition fixtures without a fluid RHS."""
import argparse
from decimal import Decimal,localcontext
import hashlib
import json
from pathlib import Path

import numpy as np
from full_moment_cell_journal import observe


def digest(a):return hashlib.sha256(np.ascontiguousarray(a).view(np.uint8)).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    rng=np.random.default_rng(1824);count=1024;n=3;ids=np.array([11,29,47],dtype=np.int64)
    i=np.arange(count)%n;j=(i+1)%n;wall=np.arange(count)%7==0;j[wall]=-1
    owners=np.column_stack([i,j]);normal=np.zeros((count,2));normal[:,0]=1.
    vl=rng.normal(size=(count,4,2));vr=rng.normal(size=(count,4,2));vr[wall]=vl[wall]*[-1.,1.]
    rhoL=.5+rng.random(count);rhoR=.5+rng.random(count);rhoR[wall]=rhoL[wall]
    wn=rng.normal(scale=.1,size=count);wn[wall]=0.;scale=.1+rng.random(count)
    unL=(vl@normal[:,:,None])[:,:,0];unR=(vr@normal[:,:,None])[:,:,0]
    wl=.25*rhoL[:,None]*np.maximum(unL-wn[:,None],0.)
    wr=.25*rhoR[:,None]*np.minimum(unR-wn[:,None],0.)
    nodes=np.zeros((count,33));nodes[:,:2]=rng.normal(size=(count,2));nodes[:,2:4]=normal;nodes[:,4]=scale
    nodes[:,5:13]=vl.reshape(count,8);nodes[:,13:21]=vr.reshape(count,8)
    nodes[:,21:25]=wl;nodes[:,25:29]=wr;nodes[:,29]=.1;nodes[:,30]=.2;nodes[wall,30]=.1
    nodes[:,31]=1.;nodes[:,32]=.05;nodes[wall,32]=0.
    geometry=dict(owners=owners,normals=normal)
    transport=dict(vl=vl,vr=vr,rl=rhoL,rr=rhoR)
    values=(wl,wr,wn);captured=dict(node_inputs=nodes,owners=owners)
    inputs=[nodes,owners,normal,vl,vr,rhoL,rhoR,wl,wr,wn,ids]
    before=[digest(x) for x in inputs]
    report=observe(geometry,transport,values,captured,ids,'synthetic-observer-fixture',.025,args.output)
    assert [digest(x) for x in inputs]==before
    with np.load(args.output/'cell_transport_journal.npz') as f:data={key:f[key] for key in f.files}
    with localcontext() as context:
        context.prec=80;expected=[Decimal(0) for _ in range(n)]
        for k,(left,right) in enumerate(owners):
            flux=Decimal.from_float(float(scale[k]))*Decimal.from_float(float(wn[k]))
            expected[left]+=flux
            if right>=0:expected[right]-=flux
        expected=np.array([np.longdouble(str(value)) for value in expected])
    error=abs(data['volume']-expected)
    assert np.all(error<=data['volume_stored_operand_arithmetic_bound'])
    for name in ('volume','mass_total','energy_internal','energy_kinetic','pressure_power'):
        bound=np.sum((128+64*data['incidence'])*np.finfo(np.longdouble).eps*data[name+'_primitive_absolute_sum'])
        assert abs(data[name].sum())<=bound,(name,str(data[name].sum()),str(bound))
    scale_sum=sum(data[name+'_primitive_absolute_sum'] for name in ('mass_total','mass_mean_advective','mass_stream_residual','mass_mesh_motion'))
    assert np.all(abs(data['mass_decomposition_difference'])<=512*np.finfo(float).eps*scale_sum)
    payload_sha=report['output_sha256']
    try:observe(geometry,transport,values,captured,ids,'synthetic-observer-fixture',.025,args.output)
    except FileExistsError:pass
    else:raise AssertionError('Observer overwrote an existing journal')
    assert hashlib.sha256((args.output/'cell_transport_journal.npz').read_bytes()).hexdigest()==payload_sha
    summary=dict(status='passed_cell_transport_journal_fixtures',synthetic_face_nodes=count,cells=n,
        input_arrays_unchanged=True,independent_decimal_volume_comparison_passed=True,
        opposite_owner_signs_and_wall_conservation_passed=True,mass_decomposition_passed=True,
        existing_output_not_overwritten=True,physical_providers=0,accepted_time_steps=0,
        limits='Observer arithmetic fixture only; real provider journal and finite-step geometry budgets remain to be captured during the approved new evolution.')
    (args.output/'report.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':main()

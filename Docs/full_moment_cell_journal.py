"""Observe actual converged face values into bounded per-cell signed budgets.

No endpoint inversion, new face evaluation, root solve or dynamics correction.
Longdouble sums are diagnostics of the stored binary64 operands.
"""
import hashlib
import json
from pathlib import Path

import numpy as np


def observe(geometry,transport,values,captured,ids,stage_key,time,output):
    output=Path(output);ids=np.asarray(ids);n=len(ids)
    assert n and np.all(np.diff(ids)>0)
    owners=geometry['owners'];nodes=captured['node_inputs'];wn=values[2]
    assert len(owners)==len(nodes)==len(wn)
    assert np.array_equal(owners,captured['owners'])
    names=('volume','mass_total','mass_mean_advective','mass_stream_residual','mass_mesh_motion',
           'energy_internal','energy_kinetic','pressure_power')
    left={key:np.zeros(n,np.longdouble) for key in names}
    right={key:np.zeros(n,np.longdouble) for key in names}
    left_absolute={key:np.zeros(n,np.longdouble) for key in names}
    right_absolute={key:np.zeros(n,np.longdouble) for key in names}
    incidence=np.zeros(n,np.int64)
    for start in range(0,len(nodes),65536):
        end=min(start+65536,len(nodes));slc=slice(start,end)
        i,j=owners[slc].T;inside=j>=0
        incidence+=np.bincount(i,minlength=n)+np.bincount(j[inside],minlength=n)
        normal=geometry['normals'][slc]
        vl64,vr64=transport['vl'][slc],transport['vr'][slc]
        unL64=(vl64@normal[:,:,None])[:,:,0];unR64=(vr64@normal[:,:,None])[:,:,0]
        chooseL=(unL64-wn[slc,None])>0;chooseR=(unR64-wn[slc,None])<0
        x=nodes[slc].astype(np.longdouble);speed=np.asarray(wn[slc],np.longdouble)
        ds=x[:,4];wl=x[:,21:25];wr=x[:,25:29]
        rl=np.asarray(transport['rl'][slc],np.longdouble);rr=np.asarray(transport['rr'][slc],np.longdouble)
        unL=unL64.astype(np.longdouble);unR=unR64.astype(np.longdouble)
        meanL=unL.sum(1)/4;meanR=unR.sum(1)/4
        countL=chooseL.sum(1);countR=chooseR.sum(1)
        residualL=(chooseL*(unL-meanL[:,None])).sum(1)
        residualR=(chooseR*(unR-meanR[:,None])).sum(1)
        mean_left=.25*rl*countL*meanL;mean_right=.25*rr*countR*meanR
        residual_left=.25*rl*residualL;residual_right=.25*rr*residualR
        effective_density=.25*(rl*countL+rr*countR)
        vl=x[:,5:13].reshape(-1,4,2);vr=x[:,13:21].reshape(-1,4,2)
        kL=.5*(vl*vl).sum(2);kR=.5*(vr*vr).sum(2)
        flux=dict(volume=ds*speed,mass_total=ds*(wl.sum(1)+wr.sum(1)),
            mass_mean_advective=ds*(mean_left+mean_right),
            mass_stream_residual=ds*(residual_left+residual_right),
            mass_mesh_motion=-ds*speed*effective_density,
            energy_internal=ds*(wl.sum(1)*x[:,29]+wr.sum(1)*x[:,30]),
            energy_kinetic=ds*((wl*kL).sum(1)+(wr*kR).sum(1)),pressure_power=ds*x[:,31]*x[:,32])
        absolute=dict(volume=abs(ds*speed),mass_total=ds*(abs(wl).sum(1)+abs(wr).sum(1)),
            mass_mean_advective=ds*(abs(mean_left)+abs(mean_right)),
            mass_stream_residual=ds*.25*(rl*(chooseL*abs(unL-meanL[:,None])).sum(1)+rr*(chooseR*abs(unR-meanR[:,None])).sum(1)),
            mass_mesh_motion=abs(ds*speed)*effective_density,
            energy_internal=ds*(abs(wl).sum(1)*abs(x[:,29])+abs(wr).sum(1)*abs(x[:,30])),
            energy_kinetic=ds*((abs(wl)*kL).sum(1)+(abs(wr)*kR).sum(1)),pressure_power=abs(flux['pressure_power']))
        for name in names:
            np.add.at(left[name],i,flux[name]);np.add.at(right[name],j[inside],flux[name][inside])
            np.add.at(left_absolute[name],i,absolute[name]);np.add.at(right_absolute[name],j[inside],absolute[name][inside])
    arrays=dict(ids=ids,incidence=incidence)
    for name in names:
        # Moving boundary gives +Vdot for outward i; conserved fluid rates
        # give minus outward flux for i. The same face has opposite j signs.
        arrays[name]=(left[name]-right[name])*(1 if name=='volume' else -1)
        arrays[name+'_primitive_absolute_sum']=left_absolute[name]+right_absolute[name]
        assert np.isfinite(arrays[name]).all()
    arrays['mass_decomposition_difference']=(arrays['mass_total']-
        (arrays['mass_mean_advective']+arrays['mass_stream_residual']+arrays['mass_mesh_motion']))
    eps=np.longdouble(np.finfo(np.longdouble).eps);operations=128+64*incidence
    gamma=operations*eps/(1-operations*eps)
    arrays['volume_stored_operand_arithmetic_bound']=gamma*arrays['volume_primitive_absolute_sum']
    file=output/'cell_transport_journal.npz'
    with file.open('xb') as f:np.savez(f,**arrays)
    report=dict(status='captured_actual_provider_cell_transport_journal',time=float(time),stage_key=stage_key,
        cells=n,face_nodes=len(nodes),source_node_inputs_sha256=hashlib.sha256(np.ascontiguousarray(nodes).view(np.uint8)).hexdigest(),
        output_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        global_volume_rate=str(arrays['volume'].sum()),
        maximum_mass_decomposition_difference=str(abs(arrays['mass_decomposition_difference']).max(initial=0.)),
        physical_providers_added=0,accepted_time_steps=0,
        limits='Direct actual wn/ds and transport samples only. Stored-operand arithmetic bound excludes geometry/speed evaluation errors. Longdouble diagnostic accumulation does not replace native rates. No finite-step attribution or time update.')
    (output/'cell_transport_journal.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

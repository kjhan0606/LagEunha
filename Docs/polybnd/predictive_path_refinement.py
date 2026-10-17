"""Nested257-node and focused65-node prescribed-path mass audit."""
import json,hashlib
from pathlib import Path
from concurrent.futures import as_completed
import numpy as np
from adaptive_coupled_mass import atomic,plain,simpson,budget_ratio


def combine_samples(whole,focus,mass,dt,horizon):
    flux=np.asarray([whole[str(i)]['mass_rate'] for i in range(257)])
    fflux=np.asarray([focus[str(i)]['mass_rate'] for i in range(65)])
    signed=simpson(flux,dt);absolute=simpson(abs(flux),dt)
    first_signed=simpson(flux[:5],dt/64);first_absolute=simpson(abs(flux[:5]),dt/64)
    focused_signed=simpson(fflux,dt/64);focused_absolute=simpson(abs(fflux),dt/64)
    corrected_signed=signed-first_signed+focused_signed;corrected_absolute=absolute-first_absolute+focused_absolute
    delta=abs(signed-simpson(flux[::2],dt))+abs(signed-simpson(flux[::4],dt))
    delta+=abs(focused_signed-first_signed)+abs(focused_signed-simpson(fflux[::2],dt/64))
    all_samples=list(whole.values())+list(focus.values())
    result=dict(rate=horizon*max(float(np.max(abs(np.asarray(s['mass_rate']))/mass)) for s in all_samples),
        generator=max(s['generator'] for s in all_samples),pressure=min(s['pressure'] for s in all_samples),
        factor_R=min(s['factor_R'] for s in all_samples),signed=corrected_signed,absolute=corrected_absolute,
        quad_delta=delta,refined_nodes=257,focus_nodes=65)
    result['ratio']=budget_ratio(result['rate'],float(np.max(corrected_absolute/mass)),float(np.max(delta/mass)),result['generator'],dt,horizon,5e-11)
    return result


def assess_refined(record,mass,horizon,provider,original_samples,cache_path,label,executor=None,record_path=None):
    from predictive_evolution_batch import path_sample,submit_path_samples
    digest=hashlib.sha256(json.dumps(plain(record),sort_keys=True,separators=(',',':')).encode()).hexdigest()
    path=Path(cache_path)
    data=json.loads(path.read_text()) if path.exists() else dict(record_sha256=digest,whole={},focus={})
    if data['record_sha256']!=digest:raise RuntimeError('refined audit cache belongs to another trial')
    whole=data['whole'];focus=data['focus']
    for i in range(65):whole[str(i*4)]=original_samples[str(i)]
    def fill(target,indices,divisor):
        missing=[i for i in indices if str(i) not in target]
        if executor is None:results=((i,path_sample(record,mass,i/divisor,provider)) for i in missing)
        else:
            futures=submit_path_samples(executor,record,mass,missing,divisor,record_path)
            results=((futures[f],f.result()) for f in as_completed(futures))
        for i,value in results:
            target[str(i)]=value
            if len(target)%16==0:atomic(path,data)
        atomic(path,data)
    fill(whole,range(257),4)
    for i in range(5):focus[str(i*16)]=whole[str(i)]
    fill(focus,range(65),64)
    result=combine_samples(whole,focus,mass,record['dt'],horizon)
    print(label,'refined257+focus65 ratio',result['ratio'],flush=True)
    return result


def combine_level(whole,focus,mass,dt,horizon):
    """Compare the last two quadrature levels after replacing the focus span."""
    flux=np.asarray([s['mass_rate'] for s in whole]);fflux=np.asarray([s['mass_rate'] for s in focus])
    count=(len(whole)-1)//64
    def integrals(values,focused):
        fine=simpson(values,dt)-simpson(values[:count+1],dt/64)+simpson(focused,dt/64)
        outside_difference=(simpson(values,dt)-simpson(values[:count+1],dt/64))-(simpson(values[::2],dt)-simpson(values[:count+1:2],dt/64))
        delta=abs(outside_difference)+abs(simpson(focused,dt/64)-simpson(focused[::2],dt/64))
        return fine,delta
    signed,signed_delta=integrals(flux,fflux);absolute,absolute_delta=integrals(abs(flux),abs(fflux))
    delta=np.maximum(signed_delta,absolute_delta);all_samples=whole+focus
    result=dict(rate=horizon*max(float(np.max(abs(np.asarray(s['mass_rate']))/mass)) for s in all_samples),
        generator=max(s['generator'] for s in all_samples),pressure=min(s['pressure'] for s in all_samples),
        factor_R=min(s['factor_R'] for s in all_samples),signed=signed,absolute=absolute+absolute_delta,
        estimated_absolute_integral=absolute,absolute_quadrature_delta=absolute_delta,quad_delta=delta,
        refined_nodes=len(whole),focus_nodes=len(focus))
    result['quadrature_ratio']=float(np.max(delta/mass))/(1e-12*dt/horizon)
    result['ratio']=budget_ratio(result['rate'],float(np.max(result['absolute']/mass)),float(np.max(delta/mass)),result['generator'],dt,horizon,5e-11)
    return result


def assess_adaptive(record,mass,horizon,provider,original_samples,cache_path,label,executor=None,min_nodes=257,max_nodes=1025,record_path=None):
    from predictive_evolution_batch import path_sample,submit_path_samples
    if max_nodes not in (1025,2049,4097) or min_nodes>max_nodes:raise ValueError('unsupported adaptive quadrature range')
    lattice=16*(max_nodes-1)
    digest=hashlib.sha256(json.dumps(plain(record),sort_keys=True,separators=(',',':')).encode()).hexdigest()
    path=Path(cache_path)
    data=json.loads(path.read_text()) if path.exists() else dict(version=2 if max_nodes==1025 else 3,record_sha256=digest,samples={})
    if data.get('version') not in (2,3) or data['record_sha256']!=digest:raise RuntimeError('adaptive quadrature cache mismatch')
    source_lattice=data.get('lattice',16384)
    if lattice%source_lattice:raise RuntimeError('quadrature cache lattice cannot be reduced')
    if source_lattice!=lattice:
        scale=lattice//source_lattice
        data=dict(version=3,record_sha256=digest,lattice=lattice,samples={str(int(k)*scale):v for k,v in data['samples'].items()})
    elif max_nodes!=1025:data['lattice']=lattice
    values=data['samples']
    for i in range(65):values[str(i*(lattice//64))]=original_samples[str(i)]
    levels=[]
    for nodes,focused in ((257,65),(513,129),(1025,257),(2049,513),(4097,1025)):
        if nodes>max_nodes:break
        whole_ids=[i*(lattice//(nodes-1)) for i in range(nodes)];focus_ids=[i*((lattice//64)//(focused-1)) for i in range(focused)]
        missing=sorted(set(whole_ids+focus_ids)-{int(k) for k in values})
        if executor is None:results=((j,path_sample(record,mass,j/(lattice/64),provider)) for j in missing)
        else:
            futures=submit_path_samples(executor,record,mass,missing,lattice/64,record_path)
            results=((futures[f],f.result()) for f in as_completed(futures))
        for j,value in results:
            values[str(j)]=value
            if len(values)%32==0:atomic(path,data)
        atomic(path,data)
        result=combine_level([values[str(j)] for j in whole_ids],[values[str(j)] for j in focus_ids],mass,record['dt'],horizon)
        levels.append(dict(nodes=nodes,focus_nodes=focused,ratio=result['ratio'],quadrature_ratio=result['quadrature_ratio']))
        print(label,'adaptive quadrature',nodes,focused,'ratio',result['ratio'],'quad ratio',result['quadrature_ratio'],flush=True)
        if nodes>=min_nodes and (result['ratio']<=1 or result['quadrature_ratio']<=1 or result['rate']>5e-9 or result['generator']>5e-9):break
    result['quadrature_levels']=levels;return result

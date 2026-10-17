"""Resumable canonical Kepler IC integration; no evolution is performed."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from check_kepler_moment_preflight import install
from kepler_disk_reference import ring_sites,physical_raw,PRESSURE
from factor_moment_evolve import rebuild


def atomic(path,value):
    p=Path(path);tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(value));tmp.replace(p)


def cached_cells(sites,cells,chunk_size,compute,stop_after=None):
    result=[]
    for first in range(0,len(cells),chunk_size):
        last=min(first+chunk_size,len(cells));p=Path(f'cells_{first:06d}.json')
        geometry=hashlib.sha256(json.dumps([np.asarray(sites[first:last]).tolist(),[np.asarray(x).tolist() for x in cells[first:last]]]).encode()).hexdigest()
        if p.exists():
            saved=json.loads(p.read_text())
            if saved['first']!=first or saved['geometry_sha256']!=geometry or len(saved['rows'])!=last-first:raise ValueError('IC chunk mismatch')
            rows=saved['rows']
        else:
            rows=[compute(i,sites[i],cells[i]) for i in range(first,last)]
            atomic(p,dict(first=first,geometry_sha256=geometry,rows=rows))
            print('Kepler initial cells',last,'/',len(cells),flush=True)
        result.extend(rows)
        if stop_after is not None and len(result)>=stop_after:return None
    return result


def compute_cell(i,site,poly):
    low,_,_=physical_raw(poly,16);raw,L,diag=physical_raw(poly,24)
    error=float(np.max(abs(raw-low))/max(np.max(abs(raw)),np.finfo(float).tiny))
    if error>1e-10:
        atomic('initial_failure.json',dict(cell=i,reference_error=error,polygon=np.asarray(poly).tolist()))
        raise ValueError('unresolved initial quadrature')
    return dict(mass=float(raw[0]),row=np.r_[site,0.,raw[1:]].tolist(),factor=L.tolist(),reference_error=error,reference_potential=diag['potential'])


def run(n,chunk_size=32,stop_after=None):
    if chunk_size<1:raise ValueError('invalid chunk size')
    install();sites=ring_sites(n,require_count=n==192)
    identity=dict(n=n,chunk_size=chunk_size,sites_sha256=hashlib.sha256(sites.tobytes()).hexdigest(),
        source_manifest_sha256=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest())
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('IC identity changed')
    else:atomic('identity.json',identity)
    if Path('geometry.json').exists():cells=[np.asarray(x) for x in json.loads(Path('geometry.json').read_text())]
    else:
        cells,_=base.geometry(sites,np.zeros(len(sites)),0.)
        atomic('geometry.json',[p.tolist() for p in cells])
    if len(cells)!=len(sites):raise ValueError('IC geometry count mismatch')
    rows=cached_cells(sites,cells,chunk_size,compute_cell,stop_after)
    if rows is None:return
    mass=np.asarray([r['mass'] for r in rows]);L=np.asarray([r['factor'] for r in rows]);original=np.r_[np.asarray([r['row'] for r in rows]).ravel(),np.zeros(4)]
    state=rebuild(original,mass,L)
    data=dict(t=0.,mass=mass.tolist(),state=state.tolist(),factors=L.tolist(),nominal_spacing_n=n,actual_cells=len(mass),official_resolution=n==192,pressure=PRESSURE,
        reference_quadrature_relative=max(r['reference_error'] for r in rows),Pi_initial_completion=float(np.max(abs(state-original))),
        reference_potential_moments=[r['reference_potential'] for r in rows])
    # Reference potential moments are IC diagnostics, not the optional evolved
    # potential_moments field used by the enriched closure.
    if Path('initial_state.json').exists():
        if json.loads(Path('initial_state.json').read_text())!=data:raise ValueError('existing initial state differs')
    else:atomic('initial_state.json',data)
    atomic('report.json',dict(initial_integration_complete=True,n=n,cells=len(mass),full_validation=False,
        scope='Initial physical moments only; density reconstruction, evolution and100orbits remain unverified.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--chunk-size',type=int,default=32);p.add_argument('--stop-after',type=int);a=p.parse_args();run(a.n,a.chunk_size,a.stop_after)

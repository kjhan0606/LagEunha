"""Compare selected native face supports/vectors/centroids and work independently."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from audit_kepler_native_radial_20261026 import ordered
from audit_kepler_wall_native_band_probe_20261029 import polygon
from audit_kepler_restart_bytes_20261029 import compare


def add(result,key,a,b):
    length=float(np.linalg.norm(b-a))
    value=np.r_[length,(b-a)[1],-(b-a)[0],length*.5*(a+b)]
    result[key]=result.get(key,np.zeros(5))+value


def analyze(root):
    results=[]
    for name in ['baseline','inside','on','outside','maximum_band']:
        r=root/name
        if int((r/'exit_status.txt').read_text())!=0:raise ValueError('Incomplete capture')
        p,_=ordered(r/'khdiag.step011580.wall_rhs_probe.dat')
        xy=np.column_stack((p['x'],p['y']));tree=cKDTree(xy)
        with (r/'face_dump_r0_call1.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        native={};independent={};work_error=0.;vector_error=0.
        owners=sorted({int(q['id_i']) for q in rows})
        if owners!=[55279,55283,55288,55291]:raise ValueError('Missing selected owner')
        for i in owners:
            poly,_,_=polygon(i,xy,tree)
            for k,a in enumerate(poly):
                b=poly[(k+1)%len(poly)]
                add(independent,(i,a[2]),np.asarray(a[:2])+xy[i],np.asarray(b[:2])+xy[i])
        for q in rows:
            if int(q['call'])!=1:raise ValueError('Unexpected force call')
            i=int(q['id_i']);j=int(q['id_j']);wall=int(q['wall'])
            a=np.array([float(q['x1']),float(q['y1'])]);b=np.array([float(q['x2']),float(q['y2'])])
            ds=np.array([float(q['dSx']),float(q['dSy'])]);vi=np.array([p['vx'][i],p['vy'][i]])
            add(native,(i,j),a,b)
            vector_error=max(vector_error,float(np.linalg.norm(ds-np.array([(b-a)[1],-(b-a)[0]]))))
            work=-float(p['pressure'][i])*float(q['udotds'])
            if wall==0 and j>=0:
                delta=xy[j]-xy[i];d=float(np.linalg.norm(delta))
                coefficient=.5*ds+float(np.linalg.norm(ds))*(.5*(a+b-xy[i]-xy[j]))/d
                work+=float(vi@((float(p['pressure'][i])-float(p['pressure'][j]))*coefficient))
            elif wall!=2 or j>=0:raise ValueError('Unexpected nonphysical ghost')
            work_error=max(work_error,abs(work-float(q['work'])))
        keys=set(native)|set(independent)
        maximum=max(float(np.max(abs(native.get(k,np.zeros(5))-independent.get(k,np.zeros(5))))) for k in keys)
        reciprocal=0.
        for (i,j),value in native.items():
            if j not in owners:continue
            other=native.get((j,i),np.zeros(5))
            # Length and first moment equal, area vector reverses.
            expected=other*np.array([1,-1,-1,1,1])
            reciprocal=max(reciprocal,float(np.max(abs(value-expected))))
        original=root.parent/'native_band_probes'/name/'khdiag.step011580.wall_rhs_probe.dat'
        parity=compare(original,r/'khdiag.step011580.wall_rhs_probe.dat')
        checks={'native_support_vectors_moments_match_independent':maximum<1.e-10,
                'native_selected_face_reciprocity':reciprocal<1.e-10,
                'native_area_vectors_match_endpoints':vector_error<1.e-12,
                'native_face_work_matches_variational_formula':work_error<1.e-12,
                'capture_does_not_change_full_RHS_state':parity['passed']}
        results.append({'case':name,'passed':all(checks.values()),'checks':checks,
            'native_rows':len(rows),'native_support_keys':len(native),'independent_support_keys':len(independent),
            'maximum_support_length_vector_first_moment_error':maximum,
            'selected_reciprocity_error':reciprocal,'area_vector_error':vector_error,
            'face_work_error':work_error,'capture_full_RHS_parity':parity})
    return {'cases':results,'passed':all(x['passed'] for x in results),
            'scope':'Four native owner IDs and their faces, including reciprocal crossing-site neighbors. Independent numerical geometry comparison, not bitwise CPU/GPU geometry equivalence. Repeated collinear segments aggregate by support; all supports retained and maximum absolute measure/moment error reported. No event or physical trajectory validation.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=analyze(args.root)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r,indent=2))
    if not r['passed']:raise SystemExit(1)

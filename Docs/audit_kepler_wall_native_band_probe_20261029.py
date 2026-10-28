"""Independent local bounded geometry, force/work and full-state probe checks."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from audit_kepler_native_radial_20261026 import ordered
from audit_kepler_bounded_varpress_adjoint_20261027 import clip
from audit_kepler_restart_bytes_20261029 import compare


def polygon(i,xy,tree):
    angles=np.arange(1024)*2*np.pi/1024
    seed=np.column_stack((3.+2.75*np.cos(angles),3.+2.75*np.sin(angles)))-xy[i]
    poly=[(float(v[0]),float(v[1]),-1000-k) for k,v in enumerate(seed)]
    distance,indices=tree.query(xy[i],k=128)
    for j in indices:
        if j==i:continue
        delta=xy[j]-xy[i]
        poly=clip(poly,delta,.5*float(delta@delta),int(j))
    vertices=np.asarray([v[:2] for v in poly])
    radius=float(np.max(np.linalg.norm(vertices,axis=1)))
    # Add every site which could cut the partial polygon. A bisector with
    # |qj-qi|>2R cannot cut any vertex inside this R ball about qi, even
    # when qi is outside the fixed polygon. No site-inside-cell assumption.
    candidates=tree.query_ball_point(xy[i],2*radius*(1+1.e-12))
    done=set(indices.tolist())
    for j in sorted((j for j in candidates if j not in done),
                    key=lambda j:float(np.sum((xy[j]-xy[i])**2))):
        delta=xy[j]-xy[i]
        poly=clip(poly,delta,.5*float(delta@delta),j)
    vertices=np.asarray([v[:2] for v in poly])
    area=.5*float(np.sum(vertices[:,0]*np.roll(vertices[:,1],-1)-
                        vertices[:,1]*np.roll(vertices[:,0],-1)))
    if not np.isfinite(area) or area<=0:raise ValueError('Nonpositive independent cell')
    delta=xy-xy[i]
    residual=delta@vertices.T-.5*np.sum(delta*delta,axis=1)[:,None]
    worst=float(residual.max())
    if worst>1.e-12:raise ValueError('Incomplete independent support search')
    return poly,area,worst


def reconstruct(i,xy,p,poly,edge_resolution=0.):
    vi=np.array([p['vx'][i],p['vy'][i]],dtype=float)
    force=np.zeros(2);vdot=0.;closure=np.zeros(2);rows=[]
    for k,a in enumerate(poly):
        b=poly[(k+1)%len(poly)];va=np.asarray(a[:2]);vb=np.asarray(b[:2])
        segment=vb-va;length=float(np.linalg.norm(segment))
        if length<=edge_resolution:continue
        outward=np.array([segment[1],-segment[0]])/length
        closure+=length*outward
        j=a[2]
        if j<0:continue
        delta=xy[j]-xy[i];d=float(np.linalg.norm(delta));normal=delta/d
        if float(outward@normal)<1-1.e-8:raise ValueError('Wrong support orientation')
        offset=.5*(va+vb-delta)
        coefficient=length*(.5*normal+offset/d)
        force+=(float(p['pressure'][i])-float(p['pressure'][j]))*coefficient
        vj=np.array([p['vx'][j],p['vy'][j]],dtype=float)
        vdot+=length*(.5*float((vi+vj)@normal)-float((vj-vi)@offset)/d)
        rows.append((j,va,vb))
    udot=-float(p['pressure'][i])*vdot
    edot=float(vi@force)+udot
    return force/float(p['mass'][i]),edot,udot,float(np.linalg.norm(closure)),rows


def analyze(root):
    stage=json.loads((root/'stage.json').read_text());cases=[]
    for cfg in stage['cases']:
        r=root/cfg['case'];rc=int((r/'exit_status.txt').read_text())
        data=r/'khdiag.step011580.wall_rhs_probe.dat'
        if cfg['expected_exit']:
            text=(r/'log').read_text(errors='replace')
            cases.append({'case':cfg['case'],'passed':rc!=0 and not data.exists() and '[GPU_TESS_INVALID]' in text,
                          'native_exit':rc,'scope':'Invalid diagnostic geometry rejected before K1 snapshot; rollback not tested here'})
            continue
        p,meta=ordered(data);source,m0=ordered(r/'khout.011580.dat')
        # The established snapshot reader exposes physical fields only;
        # use the native selected-stage TSV for compiled RHS fields.
        fields=['x','y','vx','vy','mass','volume','ie','den','pressure']
        with (r/'rk4_stage_r0.tsv').open() as stream:
            stage_rows=list(csv.DictReader(stream,delimiter='\t'))
        if len(stage_rows)!=1 or stage_rows[0]['id']!='55283' or stage_rows[0]['rhs_valid']!='1':
            raise ValueError('Missing unique native K1 RHS row')
        stage_row=stage_rows[0]
        checks={'exit_zero':rc==0,'owned_count_and_IDs_complete':len(p)==55351,
                'no_time_advance':meta['t']==m0['t'],
                'input_coordinates_and_mass_retained':all(np.array_equal(p[k],source[k]) for k in ['x','y','mass']),
                'finite_all_owned_fields':all(np.isfinite(p[k]).all() for k in fields),
                'positive_all_owned_volume_U':bool(np.all(p['volume']>0) and np.all(p['ie']>0)),
                'W0_exact':bool(np.all(p['w2']==0)),
                'no_invalid_native_diagnostic':'[RK4_INVALID_STATE]' not in (r/'log').read_text(errors='replace')}
        xy=np.column_stack((p['x'],p['y']));tree=cKDTree(xy)
        polys={};areas={};residuals={}
        poly,area,residual=polygon(55283,xy,tree)
        wanted={55283}|{a[2] for a in poly if a[2]>=0}
        for i in wanted:
            poly,area,residual=polygon(i,xy,tree)
            polys[i]=poly;areas[i]=area;residuals[i]=residual
        area_error=max(abs(areas[i]/float(p['volume'][i])-1) for i in wanted)
        acc,edot,udot,closure,faces=reconstruct(55283,xy,p,polys[55283])
        native_acc=np.array([float(stage_row['ax_hydro']),float(stage_row['ay_hydro'])])
        vi=np.array([p['vx'][55283],p['vy'][55283]],dtype=float)
        native_edot=float(stage_row['die_stored'])
        native_udot=native_edot-float(p['mass'][55283])*float(vi@native_acc)
        reciprocal=0.
        for j,a,b in faces:
            candidates=[(np.asarray(v[:2])+xy[j]-xy[55283],
                         np.asarray(polys[j][(k+1)%len(polys[j])][:2])+xy[j]-xy[55283])
                        for k,v in enumerate(polys[j]) if v[2]==55283]
            if len(candidates)!=1:raise ValueError('Missing/duplicated independent reciprocal support')
            c,d=candidates[0]
            reciprocal=max(reciprocal,float(np.linalg.norm(a-d)),float(np.linalg.norm(b-c)))
        acc_error=float(np.linalg.norm(acc-native_acc))
        edot_error=abs(edot-native_edot);udot_error=abs(udot-native_udot)
        checks.update({'independent_volumes_match':area_error<1.e-6,
                       'independent_force_matches':acc_error<1.e-9*max(1.,float(np.linalg.norm(native_acc))),
                       'independent_total_energy_work_matches':edot_error<1.e-11,
                       'independent_thermal_work_matches':udot_error<1.e-11,
                       'independent_face_closure_and_reciprocity':closure<1.e-10 and reciprocal<1.e-10})
        cases.append({'case':cfg['case'],'checks':checks,'passed':all(checks.values()),
                      'independently_checked_IDs':sorted(wanted),'area_max_relative':area_error,
                      'all_site_bisector_residual_max':max(residuals.values()),
                      'pressure_acceleration_independent':acc.tolist(),'pressure_acceleration_native':native_acc.tolist(),
                      'acceleration_error':acc_error,'total_energy_work_error':edot_error,
                      'thermal_work_error':udot_error,'face_closure':closure,'face_reciprocity_error':reciprocal})
    old=root.parent/'rhs_extension_probes/inside/khout.011581.dat'
    parity=compare(old,root/'flag_off_parity/khout.011581.dat')
    return {'cases':cases,'flag_off_full_restart_parity':parity,
            'passed':all(x['passed'] for x in cases) and parity['passed'],
            'scope':'All owned states checked for first RHS, independent complete geometry/force/work of crossing site and direct neighbors. Invalid cases reject before output. No event controller, rollback, reflection or long-time stability proof.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=analyze(args.root)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r,indent=2))
    if not r['passed']:raise SystemExit(1)

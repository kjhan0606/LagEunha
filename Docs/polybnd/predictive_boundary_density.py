"""Boundary density from numerical moments/neighbor data only; no oracle.
The full face-bin profile is fixed to inferred segment means while the bulk
LP preserves m,q,Q and checks the same near-optimal volume-L1 budget.
"""
import time
import numpy as np
from shapely.geometry import Polygon
from neighbor_density_audit import interface_prior
from moment_aligned_density_audit import moment_prior,cut_at_fraction
from boundary_layer_feasibility_audit import boundary_parts,boundary_face_matrix
from face_density_selection_audit import line_prior_mean
from face_moment_feasibility_audit import build,solve


class BoundaryDensity:
    def __init__(self,mode='neighbor',segments=64,radial=1e-7,reuse_witness=False,vector_geometry=False,core_level=2,align_core=False,lp_retry=False,construct_prior_witness=False,smooth_phase_floor=False,phase_fraction_floor=0.,refine_moment_direction=False):
        if mode not in ('neighbor','moment','global_moment'):raise ValueError('unknown density prior')
        self.mode=mode;self.segments=segments;self.radial=radial
        self.cache_key=None;self.cache=None;self.calls=0
        self.reuse_witness=reuse_witness;self.previous=[]
        self.vector_geometry=vector_geometry
        self.core_level=core_level;self.align_core=align_core
        self.lp_retry=lp_retry
        self.construct_prior_witness=construct_prior_witness
        self.smooth_phase_floor=smooth_phase_floor
        if not np.isfinite(phase_fraction_floor) or not 0<=phase_fraction_floor<.01:
            raise ValueError('invalid phase fraction regularization scale')
        if phase_fraction_floor and not smooth_phase_floor:
            raise ValueError('finite phase fraction floor requires smooth_phase_floor')
        self.phase_fraction_floor=phase_fraction_floor
        self.refine_moment_direction=refine_moment_direction

    def __call__(self,polys,edges,raw):
        key=(tuple(v.tobytes() for v in polys),raw[:,:7].tobytes())
        if key==self.cache_key:return self.cache
        started=time.monotonic();self.calls+=1;N=len(raw);S=self.segments
        means=np.array([r[0]/Polygon(v).area for r,v in zip(raw,polys)])
        low=float(means.min());high=float(means.max())
        centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in polys]);neighbors=[set() for _ in raw];faces=[[] for _ in raw]
        for edge_id,(i,j,a,b,n,length) in enumerate(edges):
            faces[i].append((edge_id,a,b))
            if j>=0:faces[j].append((edge_id,a,b));neighbors[i].add(j);neighbors[j].add(i)
        priors=[interface_prior(i,v,centers,means,neighbors,smooth_phase_floor=self.smooth_phase_floor,
            phase_fraction_floor=self.phase_fraction_floor)[1] for i,v in enumerate(polys)]
        if self.mode=='global_moment' and high-low>=1e-12:
            # Explicit single-interface-direction hypothesis for planar tests.
            # Infer from numerical global moments, never from a reference plane.
            from shapely.geometry import MultiPoint
            hull=MultiPoint(np.vstack(polys)).convex_hull
            if abs(hull.area-sum(Polygon(v).area for v in polys))>1e-10:
                raise RuntimeError('global direction requires a convex covered domain')
            vertices=np.asarray(hull.exterior.coords[:-1]);center=centers.mean(0)
            gradient=np.linalg.lstsq(centers-center,means-means.mean(),rcond=None)[0]
            norm=np.linalg.norm(gradient);seed=dict(flat=bool(norm<1e-12))
            if not seed['flat']:seed['normal']=(gradient/norm).tolist()
            global_prior=moment_prior(vertices,raw[:,:7].sum(0),means,seed,refine=self.refine_moment_direction)
            if global_prior['flat']:raise RuntimeError('global moment direction unavailable')
            normal=np.asarray(global_prior['normal']);shared=[]
            for v,p in zip(polys,priors):
                if p['flat']:shared.append(p);continue
                offset,_=cut_at_fraction(v,normal,p['fraction'])
                shared.append(dict(p,normal=normal.tolist(),offset=offset))
            priors=shared;self.last_global_prior=global_prior
        if self.mode=='moment':
            def infer(v,r,p):
                return moment_prior(v,r[:7],means,p,refine=True) if self.refine_moment_direction else moment_prior(v,r[:7],means,p)
            if self.phase_fraction_floor:
                # Fit the same normal on both sides of the upper blend endpoint.
                # Switching neighbor -> moment there would add another jump.
                priors=[p if p['flat'] else dict(infer(v,r,p),blend=p.get('blend',1.))
                    for v,r,p in zip(polys,raw,priors)]
            else:priors=[p if self.smooth_phase_floor and p.get('blend',1.)<1. else infer(v,r,p) for v,r,p in zip(polys,raw,priors)]
        self.last_attempted_priors=priors
        profiles={};fits=[];cell_reconstructions=[];next_previous=[];lp_calls=0;witness_reuses=0;prior_witnesses=0
        for i,v in enumerate(polys):
            part_options=dict(core_level=self.core_level,core_prior=priors[i] if self.align_core else None)
            if self.vector_geometry:
                from vector_boundary_geometry import triangle_moments_many,phase_averages
                parts=boundary_parts(v,self.radial,S,moment_function=triangle_moments_many,**part_options)
                geo=triangle_moments_many(parts)
                prior_values=None if priors[i]['flat'] else phase_averages(parts,priors[i],low,high)
                p=build(v,raw[i,:7],parts,priors[i],means[i],low,high,geometry_moments=geo,area_prior=prior_values)
            else:
                parts=boundary_parts(v,self.radial,S,**part_options);p=build(v,raw[i,:7],parts,priors[i],means[i],low,high)
            p['moment_row_scale']=1e4
            p['retry_numerical_failure']=self.lp_retry
            bulk=p;mappings=[];p=dict(p);p['bounds']=list(p['bounds'])
            for edge_id,a,b in faces[i]:
                T0,T1=boundary_face_matrix(v,len(parts),S,a,b,core_level=self.core_level,core_count=len(parts)-len(v)*2*S);ids=np.flatnonzero(T0)
                reverse=T1[ids[0]]>0
                target=np.array([line_prior_mean(a+(b-a)*j/S,a+(b-a)*(j+1)/S,priors[i],means[i],low,high) for j in range(S)])
                local=target[::-1] if reverse else target
                for part,value in zip(ids,local):p['bounds'][part]=(value/high,value/high)
                mappings.append((edge_id,ids,reverse,target))
            reused=None;from_prior=False
            if self.construct_prior_witness:
                from boundary_density_continuation import prior_witness
                reused=prior_witness(p)
                if reused is not None:prior_witnesses+=1;from_prior=True
            if reused is None and self.reuse_witness:
                from boundary_density_continuation import witness
                previous=self.previous[i] if len(self.previous)==N else None
                reused=witness(p,previous,v)
            if reused is None:
                bulk['keep_dual_certificate']=self.reuse_witness
                baseline=solve(bulk);lp_calls+=1
                if not baseline['feasible']:
                    self.last_failed_cell=i;self.last_failed_stage='bulk';self.last_failed_problem=bulk
                    from density_failure_diagnosis import diagnose_density_failure
                    self.last_failure_context=diagnose_density_failure(bulk,raw[i,:7],'bulk')
                    raise RuntimeError('bulk fit unavailable: '+str(dict(cell=i,core_level=self.core_level,
                        solver={k:x for k,x in baseline.items() if k!='mass_fractions'},diagnosis=self.last_failure_context)))
                fit=solve(p);lp_calls+=1
                if not fit['feasible']:
                    self.last_failed_cell=i;self.last_failed_stage='face';self.last_failed_problem=p
                    from density_failure_diagnosis import diagnose_density_failure
                    self.last_failure_context=diagnose_density_failure(p,raw[i,:7],'face')
                    raise RuntimeError('face/bulk fit unavailable: '+str(dict(cell=i,core_level=self.core_level,
                        solver={k:x for k,x in fit.items() if k!='mass_fractions'},diagnosis=self.last_failure_context)))
                u=np.asarray(fit['mass_fractions'])*raw[i,0]/p['raw'][:,0]/high
            else:
                fit,baseline,u=reused
                if not from_prior:witness_reuses+=1
            if self.reuse_witness:next_previous.append(dict(u=u,high=high,vertices=v.copy(),certificate=baseline['certificate']))
            cap_excess=max(0.,fit['objective']-baseline['objective']-1e-5)
            if cap_excess>1e-8:raise RuntimeError('numerical-only density volume cap failed: '+str(cap_excess))
            mass_fractions=np.asarray(fit['mass_fractions']);rho=mass_fractions*raw[i,0]/p['raw'][:,0]
            trace_error=0.
            for edge_id,ids,reverse,target in mappings:
                profile=rho[ids][::-1] if reverse else rho[ids]
                trace_error=max(trace_error,float(np.max(abs(profile-target)))/max(1.,high))
                profiles[(i,edge_id)]=profile
            if trace_error>1e-8:raise RuntimeError('density profile constraint failed')
            cell_reconstructions.append((parts,rho))
            fits.append(dict(moment=fit['moment_error'],minimum=fit['min_density'],maximum=fit['max_density'],
                cap_excess=cap_excess,trace_error=trace_error,objective=fit['objective'],baseline_bound=baseline['objective'],
                solver_attempts=baseline.get('solver_attempts',[])+fit.get('solver_attempts',[]),
                solver_retries=max(0,len(baseline.get('solver_attempts',[]))-1)+max(0,len(fit.get('solver_attempts',[]))-1)))
        def trace(cell,edge_id,parameter):
            j=min(S-1,max(0,int(parameter*S)))
            return float(profiles[(cell,edge_id)][j])
        diag=dict(closure_moment_error=max(f['moment'] for f in fits),
            closure_min_density=min(f['minimum'] for f in fits),closure_max_density=max(f['maximum'] for f in fits),
            closure_volume_cap_excess=max(f['cap_excess'] for f in fits),closure_trace_error=max(f['trace_error'] for f in fits),
            closure_build_seconds=time.monotonic()-started,closure_lp_calls=lp_calls,closure_witness_reuses=witness_reuses,
            closure_solver_retries=sum(f['solver_retries'] for f in fits),closure_prior_witnesses=prior_witnesses,
            closure_phase_blended_cells=sum(p.get('blend',1.)<1. for p in priors))
        self.last_profiles=profiles;self.last_priors=priors
        # Expose the already solved reconstruction for read-only diagnostics.
        self.last_cells=cell_reconstructions
        self.last_fits=fits
        self.previous=next_previous
        self.cache_key=key;self.cache=(trace,diag)
        return self.cache

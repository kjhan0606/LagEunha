"""Compact Psi density parameters and a stage-scoped bounded volume cache.

Restores the same exponential fit without Newton. Positive mass probabilities
are normalized by the original softmax arithmetic, not replaced by w*rho.
Only the requested few cells retain quadrature arrays; traces remain six-feature.
"""
from collections import OrderedDict
import numpy as np
from contact_potential_density import ContactPotentialDensity,features

class VolumeCache:
    def __init__(self,capacity=4):
        if capacity<1:raise ValueError('positive cache capacity required')
        self.capacity=int(capacity);self.stage_key=None;self.entries=OrderedDict();self.hits=0;self.misses=0
    def begin_stage(self,stage_key):
        if not isinstance(stage_key,str) or not stage_key:raise ValueError('explicit nonempty stage hash required')
        if stage_key!=self.stage_key:self.entries.clear();self.stage_key=stage_key
    def volume(self,density):
        if density.stage_key!=self.stage_key:raise RuntimeError('density volume cache stage mismatch')
        key=(density.stage_key,density.cell_id,density.order)
        if key in self.entries:
            self.hits+=1;value=self.entries.pop(key);self.entries[key]=value;return value
        self.misses+=1
        points,weights,prior=density.quadrature(density.order)
        G=density._evaluation_features(points);exponent=G@density.theta+np.log(prior);shift=float(exponent.max())
        wm=weights*np.exp(exponent-shift);Z=float(wm.sum());prob=wm/Z
        logZ=shift+np.log(Z)
        if logZ!=density.logZ:raise RuntimeError(dict(status='compact_density_normalization_not_reproduced',saved=density.logZ,restored=logZ))
        if not np.isfinite(prob).all() or np.any(prob<=0):raise RuntimeError('compact density probability unresolved')
        value=(points,weights,prior,density.mass*prob);self.entries[key]=value
        while len(self.entries)>self.capacity:self.entries.popitem(last=False)
        return value

class CompactPotentialDensity(ContactPotentialDensity):
    def __init__(self,domain,raw,prior_metadata,parameter,potential,cache,stage_key,cell_id,order=24):
        self.mass=float(raw[0]);self.area=float(parameter['area']);self.origin=np.array(parameter['origin']);self.scale=float(parameter['scale'])
        self.cut_metadata=dict(prior_metadata);self.theta=np.array(parameter['theta']);self.logZ=float(parameter['logZ']);self.potential=potential
        self.reference=np.array(parameter['reference'],dtype=np.longdouble);self.remainder_projection=np.array(parameter['projection_remainder'],dtype=np.longdouble)
        self.projection=self.reference+self.remainder_projection;self.sigma=np.longdouble(parameter['sigma']);self.rank_error=np.longdouble(parameter['rank_error']);self.Psi_target=np.longdouble(parameter['target'])
        self.residual=parameter['residual'];self.condition=parameter['condition'];self.iterations=parameter['iterations'];self.polynomial_dependent=False;self.radii=(.25,.5)
        if self.theta.shape!=(6,):raise ValueError('actual six-feature density required')
        if prior_metadata['flat']:pieces=[(domain,prior_metadata['value'])]
        else:
            n=np.array(prior_metadata['normal']);o=prior_metadata['offset']
            pieces=[(domain.clip_halfplane(n,o),prior_metadata['low']),(domain.clip_halfplane(-n,-o),prior_metadata['high'])]
        self.parts=[];self.part_priors=[]
        for part,rho in pieces:
            if part is None:continue
            self.parts.append(part);self.part_priors.append(float(rho))
        self.cache=cache;self.stage_key=stage_key;self.cell_id=int(cell_id);self.order=int(order)
        ld=np.longdouble;center=np.asarray(raw[1:3],dtype=ld)/ld(self.mass);o=self.origin.astype(ld);Q=np.asarray(raw[3:7],dtype=ld).reshape(2,2)
        second=Q/ld(self.mass)-np.outer(o,center)-np.outer(center,o)+np.outer(o,o)
        self.target_A=np.r_[ld(1),(center-o)/ld(self.scale),second[0,0]/ld(self.scale)**2,second[0,1]/ld(self.scale)**2,second[1,1]/ld(self.scale)**2]
    @property
    def points(self):return self.cache.volume(self)[0]
    @property
    def weights(self):return self.cache.volume(self)[1]
    @property
    def prior(self):return self.cache.volume(self)[2]
    @property
    def weighted_mass(self):return self.cache.volume(self)[3]

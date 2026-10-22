"""Same failed disks, exact rational algebra and scalar exp corroboration."""
from pathlib import Path
import json,struct
import mpmath as mp
from audit_shared_atlas_witness_v1 import wire
from audit_shared_pending29_projection_v1 import ball,add,mul,div,encloses,verify
from audit_shared_regular_q_transfer_v1 import rows,Q
from fractions import Fraction as F

ROOT=Path(__file__).resolve().parent

def m(x):return mp.mpf(x.numerator)/x.denominator

def main():
    run=ROOT/'shared_density_witness_v1_mpi1';c=verify(run)
    old=ROOT/'shared_baseline_fields_v2_mpi1'
    _,axes=wire(old/'regular_axis_proofs_rank0.bin')
    failed={struct.unpack_from('<5Q',a):a for a in axes if not struct.unpack_from('<Q',a,40)[0]}
    _,points=wire(run/'same_radius_density_witnesses_rank0.bin')
    param={int(r['parameter']):Q(r['value']) for r in rows(run/'original_density_parameters.csv')}
    seen=set();reflections=0
    with mp.workdps(160):
        for p in points:
            key=struct.unpack_from('<5Q',p);assert key in failed and key not in seen;seen.add(key)
            assert len(p)==672 and struct.unpack_from('<2Q',p,40)==(1,1)
            assert p[64:96]==failed[key][224:256]==p[96:128]
            radius=ball(p,64);radial=add(radius,(param[7],F(0)));encloses(ball(p,128),radial)
            assert radial[0]>radial[1]
            for at,index,increasing in ((160,4,1),(304,5,0)):
                z=div(mul((F(2),F(0)),add(radius,(-param[index],F(0)))),(param[8],F(0)))
                encloses(ball(p,at),z)
                reflected,inc=struct.unpack_from('<2Q',p,at+128)
                assert reflected==int(ball(p,at)[0]>0) and inc==increasing
                reflections+=reflected
                E=ball(p,at+32);den=add((F(1),F(0)),E)
                encloses(ball(p,at+64),den)
                assert den[0]>den[1]
                numerator=E if bool(increasing)!=bool(reflected) else (F(1),F(0))
                encloses(ball(p,at+96),div(numerator,ball(p,at+64)))
                center,spread=ball(p,at)
                if reflected:center=-center
                true_center=mp.exp(m(center))
                required=abs(m(E[0])-true_center)+true_center*mp.expm1(m(spread))
                assert required<=m(E[1])+mp.mpf('1e-130')*max(mp.mpf(1),required)
            rho=add((param[6],F(0)),mul(mul(div((param[9],F(0)),ball(p,128)),ball(p,256)),ball(p,400)))
            encloses(ball(p,448),rho)
            rr=m(radius[0]);zi=2*(rr-m(param[4]))/m(param[8]);zo=2*(rr-m(param[5]))/m(param[8])
            original=m(param[6])+m(param[9])/(rr+m(param[7]))/(1+mp.exp(-zi))/(1+mp.exp(zo))
            assert abs(original-m(ball(p,448)[0]))<=m(ball(p,448)[1])
    assert len(seen)==len(failed)==c['saved_old_failures']==c['same_disk_density_passes']==64
    print(json.dumps(dict(scope='same64 old failed radius disks; exact rational logistic identities/denominator separation/profile accounting,160digit scalar exponential and real-profile corroboration; native libm enclosure premise retained',
        independent_same_radius_density_accounting_pass=True,old_density_failure_disks=64,
        same_density_disk_passes=64,reflected_exponential_cases=reflections,
        density_profile_changed=False,graph_domains_or_precision_or_tau_changed=False,
        new_complex_field_queries=0,new_physical_quadrature_queries=0,
        unattempted_suffix_axes_still_require_full_fields=11264,
        all12_complete=False,full_quadrature_allowed=False,long_Kepler_validated=False),indent=2))

if __name__=='__main__':main()

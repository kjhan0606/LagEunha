"""Exact mask/preserved records and weighted-angular/root-gap inequalities."""
from fractions import Fraction as F
from pathlib import Path
import json
import struct
import sys
from audit_shared_atlas_witness_v1 import wire,quad
from audit_shared_regular_q_transfer_v1 import rows,sha
ROOT=Path(__file__).resolve().parent
def ball(b,at):return quad(b,at),quad(b,at+16)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[0],a[1]
def mul(a,b):return a[0]*b[0],abs(a[0])*b[1]+abs(b[0])*a[1]+a[1]*b[1]
def div(a,b):
    assert abs(b[0])>b[1];center=a[0]/b[0]
    return center,(a[1]+abs(center)*b[1])/(abs(b[0])-b[1])
def weighted_from_trace(r):
    s,eta=ball(r,32),ball(r,64);wl=(1-eta[0],F(0));wr=(eta[0],F(0))
    bl,al,cl=ball(r,160),ball(r,192),ball(r,224);br,ar,cr=ball(r,400),ball(r,432),ball(r,464)
    beta=add(mul(wl,bl),mul(wr,br));alpha=add(mul(wl,add(al,cl)),mul(wr,add(ar,cr)))
    db=add(br,neg(bl));da=add(add(ar,cr),neg(add(al,cl)))
    return add(add(beta,div(alpha,s)),mul((F(0),eta[1]),add(db,div(da,s))))
def main():
    run=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'shared_weighted_ray_v1_mpi1');assert run.parent==ROOT
    for name in ('source.sha256','input.sha256','binary.sha256','output.sha256'):
        for line in (run/name).read_text().splitlines():
            digest,path=line.split('  ',1);assert sha(Path(path))==digest
    for name in ('previous_source_integrity.txt','previous_output_integrity.txt','final_source_integrity.txt','final_input_integrity.txt'):
        lines=(run/name).read_text().splitlines();assert lines and all(s.endswith(': OK') for s in lines)
    assert (run/'native_exit_code.txt').read_text().strip()=='0';c=json.loads((run/'contract.json').read_text())
    _,old=wire(ROOT/'shared_ray_suffix_v1_mpi1/suffix_axis_proofs_rank0.bin')
    kept={(struct.unpack_from('<Q',b,24)[0],struct.unpack_from('<Q',b,32)[0]):b for b in old if struct.unpack_from('<Q',b)[0]==3203 and struct.unpack_from('<Q',b,40)[0]}
    assert len(kept)==381
    _,retained=wire(run/'unchanged_previous_passing_axes_rank0.bin');assert len(retained)==381 and all(b==kept[struct.unpack_from('<Q',b,24)[0],struct.unpack_from('<Q',b,32)[0]] for b in retained)
    _,oldplan=wire(ROOT/'shared_origin_atlas_v2_mpi1/graded_source133_atlas_rank0.bin');oldplan=[b for b in oldplan if struct.unpack_from('<Q',b)[0]==3203]
    _,plan=wire(run/'unchanged_original_plan_rank0.bin');assert plan==oldplan and len(plan)==304
    _,traces=wire(ROOT/'shared_ray_correlation_v1_mpi1/regular_q_correlations_rank0.bin');assert len(traces)==3
    trace_by_id={index:r for index,r in zip((191,192,193),traces)}
    h,new=wire(run/'new_weighted_axis_proofs_rank0.bin');assert h[3]==768
    expected=[(i,a) for i in range(2,304) for a in range(2) if (i,a) not in kept];assert len(expected)==223
    passcount=methods=checks=0;seen=[]
    for r in new:
        slot,source,target,index,axis,passed=struct.unpack_from('<6Q',r);assert (slot,source,target)==(3203,133,68000)
        key=index,axis;assert key not in kept and key not in seen;seen.append(key);passcount+=passed
        t=ball(r,144);chart=struct.unpack_from('<i',plan[index],48)[0];method=struct.unpack_from('<Q',r,592)[0]
        if axis==0 and index in trace_by_id:
            center,radius=weighted_from_trace(trace_by_id[index]);assert t[0]-t[1]<=center-radius<=center+radius<=t[0]+t[1];checks+=1
        if method==1:
            methods+=1;an,ae=ball(r,608);bn,be=ball(r,640);cn,ce=ball(r,672)
            exact=(F(9,1048576),F(0),F(1,1048576)) if chart==2 else (F(40,4194304),F(-32,4194304),F(10,4194304))
            assert an-ae<=exact[0]<=an+ae and bn-be<=exact[1]<=bn+be and cn-ce<=exact[2]<=cn+ce
            distance,gap,sqrtc,derivative=quad(r,704),quad(r,720),quad(r,736),quad(r,752)
            n,e=t
            if chart==2:distance2=n*n+9
            else:distance2=(n-F(8,5))**2+F(36,25)
            assert 0<distance and distance*distance<=distance2
            if passed:
                assert gap>0 and gap<=distance-e and 0<sqrtc and sqrtc*sqrtc<=exact[2]
                assert derivative>=abs(exact[1]+2*exact[2]*n)+2*exact[2]*e
        elif passed:
            argument=ball(r,240);assert argument[0]>argument[1]>=0
        assert method in (0,1)
    assert seen==expected[:len(seen)] and passcount==c['new_axis_passes'] and len(new)==c['new_axis_evaluations']
    assert methods==c['quadratic_root_distance_certificates'] and checks==3
    if c['actual_geometry_density_clearance_complete']:
        assert len(new)==passcount==223 and set(seen)|set(kept)|{(0,0),(0,1),(1,0),(1,1)}=={(i,a) for i in range(304) for a in range(2)}
    assert c['new_physical_quadrature_queries']==c['new_inverse_solves']==c['new_canonical_root_jobs']==c['normal_panel_replays']==0 and c['state_unchanged']
    assert not c['full_quadrature_allowed'] and not c['all12_analytic_moment_remainders_complete']
    print(json.dumps({'scope':'exact unchanged304panel/381record/mask accounting and same-weighted-function plus root-distance/derivative inequalities; retains native canonical field and <=8ulp libm premises','run':str(run),'contract':c,'independent_weighted_ray_accounting_pass':True,'exact_original_three_failed_angular_disk_checks':checks,'root_gap_certificate_count':methods,'all_axis_proofs_accounted_without_replay':c['actual_geometry_density_clearance_complete'],'mathematical_branch_review_pending':True,'full_quadrature_allowed':False},indent=2))
if __name__=='__main__':main()

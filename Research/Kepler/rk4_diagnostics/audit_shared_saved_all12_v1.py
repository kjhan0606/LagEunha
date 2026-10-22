"""Exact accounting of saved all12 analytic residuals and target subtotals.

The native same-function moment/libm/implicit domain premises remain premises.
This audit does not supply endpoint or native summation/serialization budgets.
"""
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import json,math,struct
from audit_shared_atlas_witness_v1 import quad
from audit_shared_pending29_projection_v1 import ball,add,mul,verify
from audit_shared_regular_q_transfer_v1 import rows,Q

ROOT=Path(__file__).resolve().parent
NAMES=('mass','px','py','angular','kinetic','physical_U','volume','qx','qy','polar','physical_Udot','potential')

def main():
    run=ROOT/'shared_saved_all12_v1_mpi1';c=verify(run)
    registry=ROOT/'shared_atlas_adoption_v7'
    registered=[r for r in rows(registry/'panel_evidence_registry.csv') if r['panel_kind']!='retained_group']
    assert len(registered)==90121
    data=(run/'saved_field_all12_panel_proofs_rank0.bin').read_bytes();h=struct.unpack_from('<8Q',data)
    assert h==(0x4753415749544e31,1,1,1664,90121,0,24,384) and len(data)==64+1664*90121
    coefficient=Q(c['cauchy_coefficient'])
    exact=F(2,33)*F(math.factorial(16)**2,math.factorial(32))**2*F(4,3)**32
    assert coefficient>=exact>0
    counts=Counter();sums=defaultdict(lambda:[F(0) for _ in range(12)]);methods=Counter()
    for i,r in enumerate(registered):
        p=memoryview(data)[64+1664*i:64+1664*(i+1)]
        slot,source,target,index,whole,passed,failed_axis=struct.unpack_from('<7Q',p)
        assert (slot,source,target,index)==tuple(int(r[k]) for k in ('slot_id','source','target','native_panel'))
        assert whole in (0,1) and passed==1 and failed_axis==(1<<64)-1
        if r['panel_kind']=='oblique':assert not whole
        counts[target]+=1;methods['whole_polynomial_measure' if whole else 'normalized_integrand_constant']+=1
        for k in range(12):
            center=quad(p,896+16*k);major=quad(p,1088+16*k);remainder=quad(p,1280+16*k)
            required=F(0)
            for axis in (0,1):
                z=ball(p,64+axis*384+32*k);measure=ball(p,832+axis*32)
                if whole:residual=mul(add(z,(-center,F(0))),measure)
                else:residual=add(mul(z,measure),(-center,F(0)))
                required=max(required,abs(residual[0])+residual[1])
            assert major>=required>=0 and remainder>=coefficient*major
            sums[target][k]+=remainder
    tau={(int(r['target']),r['moment']):Q(r['physical_tau']) for r in rows(ROOT/'shared_normalizers_v1_mpi4/all12_physical_tau.csv')}
    reports=[];seen=set()
    for r in rows(run/'target_analytic_subtotals.csv'):
        target=int(r['target']);k=NAMES.index(r['moment']);key=target,k
        assert key not in seen and target!=68526;seen.add(key)
        subtotal=Q(r['known_analytic_subtotal']);criterion=Q(r['ORIGINAL_tau'])
        assert subtotal>=sums[target][k] and criterion==tau[target,r['moment']]>0
        assert int(r['active_panels'])==int(r['completed_panels'])==counts[target]
        assert r['analytic_subtotal_complete']==r['analytic_subtotal_within_tau']=='1' and subtotal<=criterion
        assert r['endpoint_domain_error']==r['native_forward_sum_final_error']==r['total_error_status']=='budget_unproved'
        reports.append(dict(target=target,moment=r['moment'],analytic_subtotal_over_tau=float(subtotal/criterion)))
    assert len(seen)==276 and len(counts)==23 and sum(counts.values())==90121
    assert c['completed_all12_panel_majorants']==90121 and c['panels_with_missing_analytic_bounds']==0
    assert c['complete_target_moment_subtotals_over_original_tau']==c['target_moments_with_missing_panel_bounds']==0
    assert c['new_complex_field_queries']==c['new_real_center_queries']==c['new_physical_quadrature_queries']==0
    assert not c['full_quadrature_allowed'] and not c['endpoint_domain_error_complete'] and not c['native_forward_sum_final_complete']
    reports.sort(key=lambda r:r['analytic_subtotal_over_tau'],reverse=True)
    print(json.dumps(dict(scope='all90121 saved-field all12 residual/measure polynomial subtraction and exact Cauchy/target subtotal accounting; native same-function moment/libm/implicit premises retained; not endpoint/native/total transfer error or trajectory',
        independent_saved_all12_analytic_accounting_pass=True,active_panels=90121,
        original1047_retained_panels_replayed=0,completed_new_target_moment_subtotals=276,
        all_new_analytic_subtotals_within_original_tau=True,subtraction_methods=dict(methods),
        worst_analytic_subtotal_fractions=reports[:12],
        new_complex_or_center_or_physical_queries=0,P2_derivative_queries_needed_for_this_subtotal=0,
        endpoint_domain_error_complete=False,native_forward_sum_final_complete=False,
        full_quadrature_allowed=False,accepted_RK4_evolution=False,long_Kepler_validated=False),indent=2))

if __name__=='__main__':main()

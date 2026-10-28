"""Exact endpoint product identity for pressure and kinetic-energy changes.

This diagnoses a finite interval using saved states, without inventing its
missing stage fluxes or attributing the error to a particular update operator.
"""
import argparse,json
from pathlib import Path
import numpy as np
from compare_full_moment_checkpoint_errors import load,REGIONS
from diagnose_full_moment_subcell_energy import energies


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('a','b','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();assert not args.output.exists()
    a,ra,pa=load(args.a);b,rb,pb=load(args.b)
    assert rb['time']>ra['time'] and pa['physics_sha256']==pb['physics_sha256']
    assert np.array_equal(a['ids'],b['ids']) and np.array_equal(a['actual_mass'],b['actual_mass'])
    ld=np.longdouble;mass=a['actual_mass'].astype(ld);ea,eb=energies(a,'actual_'),energies(b,'actual_')
    ph=json.loads(Path('/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/report.json').read_text())['physics']
    const=(ld(ph['gamma'])-1)/ld(ph['P0']);ia=1/a['area'].astype(ld);ib=1/b['area'].astype(ld)
    scale=const*.5*(ia+ib)
    delta=lambda key:eb[key].astype(ld)-ea[key].astype(ld)
    terms={'volume':const*.5*(ea['thermal'].astype(ld)+eb['thermal'].astype(ld))*(ib-ia),
        'stored_total_energy':scale*delta('stored_total_energy'),
        'center_kinetic':-scale*delta('center_kinetic'),
        'affine_kinetic':-scale*delta('affine_kinetic'),
        'unresolved_kinetic':-scale*delta('unresolved_kinetic'),
        'decomposition_residual':-scale*delta('energy_decomposition_residual')}
    total=sum(terms.values());terms['total_plus_center_kinetic']=terms['stored_total_energy']+terms['center_kinetic']
    terms['all_kinetic']=terms['center_kinetic']+terms['affine_kinetic']+terms['unresolved_kinetic']
    era,erb=energies(a,'reference_'),energies(b,'reference_')
    reference_delta=lambda key:erb[key].astype(ld)-era[key].astype(ld)
    reference_terms={'volume':const*.5*(era['thermal'].astype(ld)+erb['thermal'].astype(ld))*(ib-ia),
        'stored_total_energy':scale*reference_delta('stored_total_energy'),
        'center_kinetic':-scale*reference_delta('center_kinetic'),
        'affine_kinetic':-scale*reference_delta('affine_kinetic'),
        'unresolved_kinetic':-scale*reference_delta('unresolved_kinetic'),
        'decomposition_residual':-scale*reference_delta('energy_decomposition_residual')}
    reference_total=sum(reference_terms.values())
    reference_terms['total_plus_center_kinetic']=reference_terms['stored_total_energy']+reference_terms['center_kinetic']
    reference_terms['all_kinetic']=reference_terms['center_kinetic']+reference_terms['affine_kinetic']+reference_terms['unresolved_kinetic']
    error_terms={key:terms[key]-reference_terms[key] for key in terms}
    observed=(b['actual_pressure'].astype(ld)-a['actual_pressure'].astype(ld))/ld(ph['P0'])
    reference_observed=(b['reference_pressure'].astype(ld)-a['reference_pressure'].astype(ld))/ld(ph['P0'])
    error_change=b['relative_pressure_error'].astype(ld)-a['relative_pressure_error'].astype(ld)
    regions={}
    for name in REGIONS:
        mask=a[name]&b[name]
        if not mask.any():continue
        w=mass[mask];obs=observed[mask]
        def stats(value):
            x=value[mask];k=int(np.argmax(abs(x)));norm=np.sqrt(np.sum(w*x*x)*np.sum(w*obs*obs))
            return dict(RMS=float(np.sqrt(np.sum(w*x*x)/w.sum())),signed_mean=float(np.sum(w*x)/w.sum()),
                maximum=float(abs(x[k])),maximum_signed=float(x[k]),maximum_id=int(a['ids'][mask][k]),
                alignment_with_actual_pressure_change=float(np.sum(w*x*obs)/norm) if norm>0 else None)
        regions[name]=dict(cells=int(mask.sum()),actual_relative_pressure_change=stats(observed),
            reference_error_change=stats(error_change),terms={k:stats(v) for k,v in terms.items()},
            actual_cell_reference_terms={k:stats(v) for k,v in reference_terms.items()},
            error_terms_after_actual_cell_reference_subtraction={k:stats(v) for k,v in error_terms.items()},
            reference_product_identity_residual=stats(reference_total-reference_observed),
            error_product_identity_residual=stats(total-reference_total-observed+reference_observed),
            product_identity_residual=stats(total-observed),
            relative_error_vs_actual_pressure_change_residual=stats(error_change-observed))
    args.output.mkdir();np.savez(args.output/'finite_signed_terms.npz',ids=a['ids'],observed=observed,total=total,
        **terms,**{'reference_'+key:value for key,value in reference_terms.items()},
        **{'error_'+key:value for key,value in error_terms.items()})
    report=dict(status='completed_saved_finite_pressure_energy_budget',start=ra['time'],end=rb['time'],
        source_A=pa,source_B=pb,regions=regions,physical_providers=0,accepted_time_steps=0,thresholds_added=False,
        identity='Delta(T/V)=mean(1/V)*Delta(T)+mean(T)*Delta(1/V); T=E-Kcm-Kaffine-Kunresolved-residual',
        limits='Endpoint algebra only. Energy terms can strongly cancel. This is not a finite-step geometry-volume defect, integrated flux budget, causal attribution, or time-integrator comparison. Minimum-rotation lower bound is never added as a separate energy.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],cavity=regions['cavity'],inner_contact=regions['inner_contact'])),flush=True)


if __name__=='__main__':main()

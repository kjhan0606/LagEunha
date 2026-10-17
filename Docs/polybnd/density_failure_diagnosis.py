"""Distinguish contradictory bounds, finite basis failure and solver failure."""
import numpy as np
from density_moment_bounds import density_moment_bounds
from bounded_moment_feasibility import minimum_moment_residual


def diagnose_density_failure(problem,target,stage='bulk'):
    context=dict(stage=stage,classification='undetermined_reconstruction_failure')
    low=.99*problem['low'];high=1.01*problem['high']
    try:
        certificate=density_moment_bounds(np.asarray(problem['raw']).sum(0),target)
        lo=certificate['density_eigen_min'];hi=certificate['density_eigen_max']
        context['density_bound_certificate']=dict(model_lower=low,model_upper=high,
            required_lower_at_most=lo,required_upper_at_least=hi)
        tolerance=1e-7*max(1.,abs(low),abs(high))
        if lo<low-tolerance or hi>high+tolerance:
            context['classification']='density_bounds_incompatible_with_moments'
            context['proof_scope']='Necessary squared-polynomial bound, independent of the finite subcell basis.'
            return context
    except Exception as exc:context['bound_diagnostic_error']=repr(exc)
    try:
        proof=minimum_moment_residual(problem['A'],problem['b'],problem['bounds'][:problem['k']])
        context['minimum_residual_diagnosis']={k:v for k,v in proof.items() if not isinstance(v,np.ndarray)}
        if proof['classification']=='finite_basis_infeasible_certificate':
            context['classification']='finite_basis_infeasible_certificate' if stage=='bulk' else 'finite_basis_and_trace_infeasible_certificate'
        elif proof['classification']=='moment_witness':
            context['classification']='moment_gate_witness_solver_unresolved'
        context['proof_scope']='Same finite basis and density box, moment gate1e-8; fixed face values included at the face stage. A gate witness does not prove exact equality, full objective feasibility or the volume budget.'
    except Exception as exc:context['residual_diagnostic_error']=repr(exc)
    return context

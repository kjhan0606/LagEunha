"""Audit helper: restore an already fitted density without Newton iterations."""
import ast
import inspect
import textwrap
import numpy as np
import contact_prior_density as density_module


def restore_density(poly,raw,metadata,saved_fit):
    """Original constructor assembly/evaluation, with stored converged theta.

    Only the Newton loop is replaced. Original quadrature, target, exponential
    normalization and probability calculations remain the constructor's AST.
    """
    tree=ast.parse(textwrap.dedent(inspect.getsource(density_module.ContactPriorDensity.__init__)))
    fn=tree.body[0];fn.name='_restore'
    fn.args.args.extend([ast.arg(arg='saved_theta'),ast.arg(arg='saved_iterations')])
    # New arguments must be required; the old order default is supplied by us.
    fn.args.defaults=[]
    index=next(i for i,node in enumerate(fn.body) if isinstance(node,ast.For)
               and isinstance(node.target,ast.Name) and node.target.id=='iteration')
    statements=ast.parse('theta=np.asarray(saved_theta,dtype=float).copy()\nresidual,H,logZ,prob=evaluate(theta)\nnorm=float(abs(residual).max())\niteration=int(saved_iterations)-1').body
    fn.body[index:index+1]=statements
    namespace=dict(density_module.__dict__)
    exec(compile(ast.fix_missing_locations(tree),'<stored-theta density restoration>','exec'),namespace)
    result=density_module.ContactPriorDensity.__new__(density_module.ContactPriorDensity)
    r=np.asarray(raw)
    namespace['_restore'](result,poly,r[0],r[1:3],r[3:7].reshape(2,2),metadata,12,
                          saved_fit['theta'],saved_fit['iterations'])
    assert result.residual==saved_fit['feature_residual']
    return result

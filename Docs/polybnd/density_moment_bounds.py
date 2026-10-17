"""Necessary pointwise density bounds from mass and position moments.

For phi=(1,(x-cx)/s,(y-cy)/s), L<=rho<=U implies
L*G <= R <= U*G in the positive-semidefinite order, where G=int phi phi^T
and R=int rho phi phi^T. Generalized eigenvectors supply nonnegative squared
polynomial witnesses. This diagnoses input bounds; it never changes moments.
"""
import numpy as np
from scipy.linalg import eigh


def moment_matrix(raw,center,scale):
    raw=np.asarray(raw);m=raw[0];q=raw[1:3];Q=raw[3:7].reshape(2,2)
    first=(q-m*center)/scale
    second=(Q-np.outer(center,q)-np.outer(q,center)+m*np.outer(center,center))/scale**2
    matrix=np.empty((3,3));matrix[0,0]=m;matrix[0,1:]=first;matrix[1:,0]=first;matrix[1:,1:]=second
    return (matrix+matrix.T)/2


def density_moment_bounds(geometry,target):
    geometry=np.asarray(geometry);target=np.asarray(target)
    if geometry[0]<=0:raise ValueError('positive cell area required')
    center=geometry[1:3]/geometry[0];scale=np.sqrt(geometry[0])
    gram=moment_matrix(geometry,center,scale);weighted=moment_matrix(target,center,scale)
    values,vectors=eigh(weighted,gram)
    return dict(density_eigen_min=float(values[0]),density_eigen_max=float(values[-1]),
        center=center,scale=scale,lower_polynomial=vectors[:,0],upper_polynomial=vectors[:,-1],
        upper_gram_norm=float(vectors[:,-1]@gram@vectors[:,-1]))

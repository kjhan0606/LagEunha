"""Differentiate SciPy DOP853 dense-output polynomial, independently of RHS."""
import numpy as np

def segment_derivative(segment,t):
    x=(t-segment.t_old)/segment.h
    value=np.zeros_like(segment.y_old);derivative=np.zeros_like(value)
    for i,f in enumerate(reversed(segment.F)):
        value=value+f
        if i%2==0:
            derivative=derivative*x+value;value=value*x
        else:
            derivative=derivative*(1-x)-value;value=value*(1-x)
    return derivative/segment.h

def path_derivative(solution,t):
    i=int(np.clip(np.searchsorted(solution.ts,t,side='left')-1,0,len(solution.interpolants)-1))
    return segment_derivative(solution.interpolants[i],t)

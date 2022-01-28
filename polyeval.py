import numpy as np

def polyeval(pcoefs,z):
    """
    Vectorized polynomial evalutation using horner's method.
    pcoefs has shape (P,K) where P is the number of coefficients in the
    polynomial and K is the number of polynomials. pcoefs[0] is the coefficient
    of the argument to the greatest power, pcoefs[1] the second greatest power
    etc.
    z has shape (N,) where N is the number of evaluations to compute
    the result has shape (K,N) and represents N evaluations of K polynomials.
    """
    p_x=pcoefs[0][:,None]
    for p in pcoefs[1:]:
        p_x=p_x*z+p[:,None]
    return p_x

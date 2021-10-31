# Some known discrete fourier transforms

import numpy as np
from common import cast_array

def dirichlet(k,W,N,normalized=True):
    C=1
    if normalized:
        C=1/W
    ret = C*np.sin((np.pi*W*k)/N)/np.sin(np.pi*k/N)
    ret[(k%N) == 0] = C*W
    return ret

def dirichlet_dk(k,W,N,normalized=True):
    C=1
    if normalized:
        C=1/W
    ret = C*W*np.pi/(N*np.sin(np.pi*k/N))*(np.cos(np.pi*W*k/N) - np.sin(np.pi*W*k/N)*np.cos(np.pi*k/N)/(W*np.sin(np.pi*k/N)))
    ret[k == 0] = 0
    return ret

def sum_of_cos_dft(
    k, # the bins at which to evaluate
    A, # a 1d vector of coefficients
    # the period of the cosine whose pth harmonics (integer multiples) 
    # correspond to the coefficients A[p] (see some_sig.sum_of_cos)    
    L, 
    # the length of the non-zero part of the window, odd if a purely real
    # transform is desired
    W,
    # the length of the DFT
    N
):
    A=cast_array(A)
    r = np.zeros_like(k).astype(A.dtype)
    B = N/L # bin width
    D = lambda k,W,N: dirichlet(k,W,N,normalized=False)
    for  p, a in enumerate(A):
        r += a*0.5*(D(k-p*B,W,N)+D(k+p*B,W,N))
    return r


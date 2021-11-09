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
    ret[(k%N) == 0] = 0
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
    N,
    # a dirichlet-like function
    D = lambda k,W,N: dirichlet(k,W,N,normalized=False)
):
    A=cast_array(A)
    r = np.zeros_like(k).astype(A.dtype)
    B = N/L # bin width
    for  p, a in enumerate(A):
        r += a*0.5*(D(k-p*B,W,N)+D(k+p*B,W,N))
    return r

def mod_sum_of_cos_dft(v,k,A,L,W,N):
    k_v=N*v
    return sum_of_cos_dft(k-k_v,A,L,W,N)

def mod_sum_of_cos_dft_k(k0,k,A,L,W,N):
    return sum_of_cos_dft(k-k0,A,L,W,N)

def sum_of_cos_dft_dk(k,A,L,W,N):
    return sum_of_cos_dft(
        k,A,L,W,N,D=lambda k,W,N: dirichlet_dk(k,W,N,normalized=False))

def mod_sum_of_cos_dft_dk(k0,k,A,L,W,N):
    # like sum_of_cos_dft_dk but of a signal modulated by a sinusoid of bin
    # frequency k0
    return sum_of_cos_dft_dk(k-k0,A,L,W,N)

def mod_sum_of_cos_dft_df(f0,Fs,f,A,L,W,N):
    # Like mod_sum_of_cos_dft_dk but in units of frequency (modulated by a
    # sinusoid of bin frequency f0 and with a sampling rate of Fs)
    v0=f0/Fs
    k0=v0*N
    k=f/Fs*N
    return mod_sum_of_cos_dft_dk(k0,k,A,L,W,N)

def _multi(B,V,k,A,L,W,N,fun):
    A=cast_array(A)
    w=np.zeros(len(k),dtype=A.dtype)
    for b,v in zip(B,V):
        w += b*fun(v,k,A,L,W,N)
    return w

def multi_mod_sum_of_cos_dft_k(B,K0,k,A,L,W,N):
    return _multi(B,K0,k,A,L,W,N,mod_sum_of_cos_dft_k)

def multi_mod_sum_of_cos_dft_dk(B,K0,k,A,L,W,N):
    return _multi(B,K0,k,A,L,W,N,mod_sum_of_cos_dft_dk)

def normalize_sum_of_cos_A(A,L,W,N,normalized=False):
    """
    Multiply all elements of A with a constant such that the fourier transform
    of a sum of cosine window using A as the coefficients will be 1 at k=0
    This is confusing: normalized refers to the "normalized" argument of
    dirichlet, which there causes it to be divided by W. So it depends what was
    passed to dirichlet when the sum-of-cosine window's FT was computed (using
    the definitions, not the DFT). For all the *sum_of_cos* functions above,
    normalized is set to False by default (see the D argument of sum_of_cos_dft).
    """
    D = lambda k: dirichlet(k,W,N,normalized=normalized)
    B = N/L # bin width
    A=cast_array(A)
    C=np.sum(A*D(B*np.arange(len(A))))
    A/=C
    return A



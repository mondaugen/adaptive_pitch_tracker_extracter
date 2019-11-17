import numpy as np
from scipy import signal

def a_to_r(a):
    """ Convert feedback filter coefficients to reflection coefficients. This assumes you
    include the coefficient for y[n], where it divides all the coefficients by
    it (so the first coefficient is 1), then discards the y[n] coefficient
    (because it knows it is 1). """
    a=a.copy()
    a=a/a[0]
    a=a[1:]
    p=len(a)
    r=np.zeros_like(a)
    r[-1]=a[-1]
    for j in np.arange(p-1)[::-1]:
        s=1/(1-r[j+1]*np.conj(r[j+1]))
        a[:j+1]=s*(a[:j+1]-r[j+1]*np.conj(a[:j+1][::-1]))
        r[j]=a[j]
    return r

def r_to_A(r):
    """
    Convert reflection coefficients to feedback coefficients. This returns all
    the intermediate filter orders along the rows, so the first row will be
    [1,0,0, .., 0], the second row [1,a_1[1], 0, 0, ..., 0], etc. So if you just
    want the highest filter order obtainable with these r coefficients just take
    A[-1,:] of the resulting matrix A.
    """
    P=len(r)
    A=np.zeros((P+1,P+1),dtype=r.dtype)
    A[0,0]=1
    for i in np.arange(1,P+1):
        A[i,:i+1]=A[i-1,:i+1]+r[i-1]*np.conj(A[i-1,:i+1][::-1])
    return A

def b_A_to_c(b,A):
    """
    Convert feedforward filter coefficients to 'C coefficients'. See p. 306 in
    Hayes book. These are used for implementing a lattice filter with poles and
    zeros.
    """
    assert(len(b)==A.shape[0])
    # This is to match the notation in Hayes
    q=len(b)-1
    c=np.zeros_like(b)
    c[q]=b[q]
    for k in np.arange(q)[::-1]:
        c[k]=b[k]-np.sum(c[k+1:q+1]*np.conj(np.diag(A,-k)[1:]))
    return c

def c_A_to_b(c,A):
    """ Convert 'C coefficients' to feedforward coefficients. """
    b=np.zeros_like(c)
    q=len(b)-1
    for k in np.arange(q+1):
        b[k]=np.sum(c[k:q+1]*np.conj(np.diag(A,-k)))
    return b

def c_r_to_b_a(c,r):
    A=r_to_A(r)
    b=c_A_to_b(c,A)
    return (b,A[-1,:])

def b_a_to_r_c(b,a):
    r=a_to_r(a)
    A=r_to_A(r)
    c=b_A_to_c(b,A)
    return (r,c)

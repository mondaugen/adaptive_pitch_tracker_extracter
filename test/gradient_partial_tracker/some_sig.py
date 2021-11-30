# Some signals

import numpy as np
from scipy import signal
from common import cast_array

j=complex('j')

def dk_ramp(N,p):
    return np.power(np.concatenate((np.arange(N//2),np.arange(N//2)-N//2)),p)

def dk_scale(N,p):
    return (-j*2*np.pi/N)**p

def multiply_ramp(x,N,p):
    return dk_ramp(N,p)*dk_scale(N,p)*x

def rectangular(n,W,N):
    ret = np.zeros(len(n))
    ret[(2*n < W)|(2*(N-n)<W)] = 1./W
    return ret

def sinusoid_about_0(k,N):
    n=np.concatenate((np.arange(N//2),np.arange(-N//2,0)))
    return np.exp(j*2*np.pi/N*np.multiply.outer(k,n))

def complex_chirp(*args):
    raise NotImplementedError
    # TODO: broken
    ret=signal.chirp(*args)
    return np.exp(j*np.arccos(ret))

def comb(t,v0,t1,v1,method='linear',phi=0,vertex_zero=True,v_max=0.5,partial_amp=lambda p: 1./p):
    V0=np.arange(v0,v_max,v0)
    V1=np.arange(v1,v_max,v1)
    min_len=min(len(V0),len(V1))
    V0=V0[:min_len]
    V1=V1[:min_len]
    X=np.zeros((len(V0),len(t)),dtype='complex128')
    for i, (v0_,v1_) in enumerate(zip(V0,V1)):
        X[i,:]=signal.chirp(t,v0_,t1,v1_,method,phi,
                            vertex_zero) * partial_amp(i+1)
    ret=X.sum(axis=0)
    return ret
        
def sum_of_cos(
    # array of cosine coefficients
    a,
    # the period of the window. This means that a[1] is the coefficient of a
    # sampled cosine of angular velocity v whose cycle takes L samples, i.e.,
    # cos(0) == cos(vL)
    L,
    # the length of boxcar whose center as at sample time 0. This should be odd
    # if a purely real fourier transform is desired.
    W,
    # the length of the signal, usually a power of 2 > W (so the DFT of the
    # resulting signal can be taken easily)
    N,
    # If centered_at_time_zero is False, then the window is centered at time N//2
    centered_at_time_zero=True
):
    a=cast_array(a)
    Q=len(a)
    x=np.zeros(N)
    l=np.arange(-(W//2),W//2+1)
    w=np.sum(a[:,None]*np.cos(
            2*np.pi/L*np.multiply.outer(np.arange(Q),l)),axis=0)
    x[:W//2+1]=w[-(W//2+1):]
    x[-(W//2):]=w[:W//2]
    if not centered_at_time_zero:
        return np.concatenate((x[N//2:],x[:N//2]))
    return x

def mod_sum_of_cos(v,a,L,W,N):
    l=v*N
    # the modulator must be conjugate-symmetric around time-0 (the centre of
    # the boxcar)
    nl=np.arange(W//2+1)
    nr=np.arange(-(W//2),0)
    w=sum_of_cos(a,L,W,N).astype('complex128')
    w[:W//2+1]*=np.exp(j*2*np.pi*l*nl/N)
    w[-(W//2):]*=np.exp(j*2*np.pi*l*nr/N)
    return w

def _multi(B,V,A,L,W,N,fun):
    w=np.zeros(N,dtype='complex128')
    for b,v in zip(B,V):
        w += b*fun(v,A,L,W,N)
    return w

def multi_mod_sum_of_cos(B,V,A,L,W,N):
    return _multi(B,V,A,L,W,N,mod_sum_of_cos)

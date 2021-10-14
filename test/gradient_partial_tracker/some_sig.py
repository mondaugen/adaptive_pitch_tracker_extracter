# Some signals

import numpy as np
from scipy import signal

j=complex('j')

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
        

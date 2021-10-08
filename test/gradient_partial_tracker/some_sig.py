# Some signals

import numpy as np

j=complex('j')

def rectangular(n,W,N):
    ret = np.zeros(len(n))
    ret[(2*n < W)|(2*(N-n)<W)] = 1./W
    return ret

def sinusoid_about_0(k,N):
    n=np.concatenate((np.arange(N//2),np.arange(-N//2,0)))
    return np.exp(j*2*np.pi/N*np.multiply.outer(k,n))

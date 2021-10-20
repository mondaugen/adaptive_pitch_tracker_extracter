# Some known discrete fourier transforms

import numpy as np

def dirichlet(k,W,N,normalized=True):
    C=1
    if normalized:
        C=1/W
    ret = C*np.sin((np.pi*W*k)/N)/np.sin(np.pi*k/N)
    ret[k == 0] = 1
    return ret

def dirichlet_dk(k,W,N,normalized=True):
    C=1
    if normalized:
        C=1/W
    ret = C*W*np.pi/(N*np.sin(np.pi*k/N))*(np.cos(np.pi*W*k/N) - np.sin(np.pi*W*k/N)*np.cos(np.pi*k/N)/(W*np.sin(np.pi*k/N)))
    ret[k == 0] = 0
    return ret



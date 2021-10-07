# discrete fourier transforms and their derivatives

import numpy as np

j=complex('j')

def dk_ramp(N,p):
    return np.power(np.concatenate((np.arange(N//2),np.arange(N//2)-N//2)),p)

def dk_scale(N,p):
    return (-j*2*np.pi/N)**p

class dft_dk:
    """ The pth derivative of a DFT w.r.t. bin k """
    def __init__(self,N,p=1):
        self.N=N
        self.p=p
        self.n=dk_ramp(N,p)
        self.scale=dk_scale(self.N,self.p)
    def __call__(self,x):
        nx=self.n*x
        return np.fft.fft(self.scale*nx)

# TODO:
# An optimized version of this might exist where the inner-product is done in
# the frequency domain.
# To do that, does it work to use the inner product between the fourier
# transform of the zero padded signal x where every Mth value (if zero padding
# by M*N) is from the the original fourier transform (length N), and a sinusoid
# convolved with the M*N length Fourier transform of a low side-lobed window
# (e.g., blackman)?
def dft_bin(x,k,p=0):
    """ Given the signal x, find the DFT's value at k. k must
    be a vector and x is multiplied by a matrix. """
    N=x.shape[0]
    Q=np.hstack((
    np.exp(-j*2*np.pi/N*np.multiply.outer(k,np.arange(0,N//2))),
    np.exp(-j*2*np.pi/N*np.multiply.outer(k,np.arange(-N//2,0))),
    ))
    return Q@(dk_scale(N,p)*dk_ramp(N,p)*x)

def ps_dk(X,dX,extract=np.real):
    """ Power spectrum derivative. """
    return extract(X*np.conj(dX)+dX*np.conj(X))

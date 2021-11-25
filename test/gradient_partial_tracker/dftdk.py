# discrete fourier transforms and their derivatives

import numpy as np
from functools import partial
from scipy import interpolate, signal
from some_ft import normalize_sum_of_cos_A
from some_sig import sum_of_cos, multi_mod_sum_of_cos, dk_ramp, dk_scale, multiply_ramp

j=complex('j')

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
    """ Given the signal x, find the value of the DFT's pth derivative at k. k
    must be a vector and x is multiplied by a matrix. """
    N=x.shape[0]
    Q=np.hstack((
    np.exp(-j*2*np.pi/N*np.multiply.outer(k,np.arange(0,N//2))),
    np.exp(-j*2*np.pi/N*np.multiply.outer(k,np.arange(-N//2,0))),
    ))
    return Q@(dk_scale(N,p)*dk_ramp(N,p)*x)

def dft_bin_grad(p):
    return partial(dft_bin,p=p)

def ps_dk(X,dX,extract=np.real):
    """ Power spectrum derivative. """
    return extract(X*np.conj(dX)+dX*np.conj(X))

def log_ps_dk(X,dX,extract=np.real):
    return extract(np.conj(dX)/np.conj(X) + dX/X)

def dft_bin_log_ps_dk(x,k):
    X=dft_bin(x,k)
    dX=dft_bin(x,k,p=1)
    return log_ps_dk(X,dX)

def gradient_ascent_step(x,k0,mu,grad=dft_bin_grad(1)):
    k1 = k0 + mu * np.real(grad(x,k0))
    return k1

def gradient_ascent_step_harm_lock(x,k0,mu,grad=dft_bin_grad(1),grad_weight='equal'):
    g=np.real(grad(x,k0))
    if grad_weight == 'equal':
        step=g.mean()
    elif grad_weight == '1/p':
        step=(g/(1+np.arange(len(g)))).mean()
    k_step = mu * step * (1+np.arange(len(k0)))
    k1 = k0 + k_step
    return k1

def newton_ascent_step(x,k0,grad=dft_bin_grad(1),grad2=dft_bin_grad(2)):
    k1 = k0 + np.real(grad(x,k0)/grad2(x,k0))
    return k1

def self_adjusting_ga_step(x,k0,step,default_step,mu,f=dft_bin_grad(0),grad=dft_bin_grad(1)):
    """
    Step by "step", does the gradient shrink?
    If it doesn't, but we moved up, then we're stepping in the correct direction
    (if the step is small and the function is not changing too fast). If we
    didn't move up, move the other way. In these two cases, just move by step
    because we can't make any better estimate of the step size.
    If it shrank, then we can use the amount it shrank and the distance we
    stepped to try to determine the step we must take to get the gradient to 0
    (assuming it is evolving linearly, i.e., the function is a parabola).
    default_step is positive by convention
    returns the next k and the revised step
    """
    if step == 0:
        step = default_step
    k1 = k0 + step
    if np.sign(step)*grad(x,k0) < np.sign(step)*grad(x,k1):
        # gradient grew
        if f(x,k1) > f(x,k0):
            # but the function also grew, so we're going in the correct direction
            # but we can't make any estimates on step
            step = np.sign(step)*default_step
        else:
            # function shrank, so we're going the wrong way
            step = -(default_step*np.sign(step))
            k1 = k0 # revert useless step
    else:
        # how far does it seem we have to step to make the gradient 0?
        d_grad = grad(x,k1) - grad(x,k0)
        d_grad_d_k = d_grad / step
        if d_grad_d_k != 0:
            step = -mu * np.real(grad(x,k0) / d_grad_d_k)
        k1 = k0 + step
    return (step,k1)

def odd_ceil(x):
    n=np.ceil(x)
    n[n % 2 == 0] += 1
    return n

def odd_round_positive(x):
    """
    Works only if x is positive.
    If x is rounded to an odd number we have nothing to do.
    If x is rounded to an even number then it has 1 added to it if x - round(x)
    is positive and -1 added if x - round(x) is negative.
    """
    if np.any(x < 0):
        raise ValueError
    n=np.round(x)
    d=x-n
    msk=n % 2 == 0
    n[msk] += np.sign(d[msk]) + (d[msk] == 0)
    return n

def window_scale(w):
    return w / w.sum()

class harm_grad_td:
    """
    Find a gradient in the frequency domain at bin k0 by using a linear
    combination of the gradient at harmonically spaced frequencies (k0, 2*k0,
    etc.).
    For example, this can be used to find the average gradient (the coefficients
    are 1/number_of_harmonics in this case).
    The gradients are found using inner-products with sinusoids windowed by
    sum-of-cosine windows.
    """
    @staticmethod
    def default_harm_sig(B=np.ones(10)/10.):
        def _synth(k0,L,W,N):
            p_max=len(B)
            p=np.arange(p_max)+1
            v=k0/N*p
            # TODO: implement a version that uses interpolation
            return multi_mod_sum_of_cos(B,v,[1],L,W,N)
        return _synth
        
    def __init__(self,
        A,L,W,N,
        # A function synthesizing a harmonic signal with fundamental k0
        # using the L, W, N parameters passed to __init__.
        harm_sig,
        normalize_A=True,
    ):
        """
        See some_ft.sum_of_cos_dft, some_ft.sum_of_cos_dft_dk for descriptions
        of the function arguments.
        """
        if normalize_A:
            self.A=normalize_sum_of_cos_A(A,L,W,N)
        else:
            self.A=A
        self.L=L
        self.W=W
        self.N=N
        self._w=sum_of_cos(self.A,self.L,self.W,self.N)
        self.harm_sig=harm_sig
    def _kern(self,k0):
        s=self.harm_sig(k0,self.L,self.W,self.N)
        return s*self._w
    def dX_dk_p(self,x,k0,p):
        return (np.conj(self._kern(k0))*multiply_ramp(x,self.N,p)).sum()
    def X(self,x,k0):
        return self.dX_dk_p(x,k0,0)
    def dX_dk(self,x,k0):
        return self.dX_dk_p(x,k0,1)

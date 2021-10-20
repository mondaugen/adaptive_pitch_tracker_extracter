# discrete fourier transforms and their derivatives

import numpy as np
from functools import partial
from scipy import interpolate, signal

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

class harm_grad_tracker:
    """
    A partial tracker that tracks partials related by the harmonic series.
    To use, initialize with the lowest and highest frequency you are interested
    in tracking, this will initialize the wavetables.
    Then to make analyses, simply call the dft_bin functions on signals.
    """
    def __init__(self,
        min_f=55, # A1
        max_f=110, # A2
        window=lambda N: window_scale(signal.get_window('blackman',N,fftbins=False)),
        q=lambda p: 1 + (p-1)/10,
        sample_rate=44100
    ):
        """
        Initialize the wavetables for doing gradient ascent tracking of harmonic
        signals.
        min_f is the frequency (Hz) of the lowest partial (or fundamental)
        we are interested in tracking.
        max_f is the frequency of the highest fundamental we are
        interested in tracking
        window is a function accepting a number a returning a window of that
        length. It is recommended that the window be symmetric and have its
        maximum at floor(N/2) for odd N.
        q is a function accepting a partial number and returning a "Q" factor.
        Here Q is simply interpreted as a divider of maximum window length. So
        if the window length determined for min_f is N, then the second partial
        would have length N*q(2) etc. In fact q(p) can be any positive number
        for p in [1, 2, ...] so N is actually determined using min_f, the
        sample rate and q(1).
        """
        # window sizes
        max_p=int(np.floor(sample_rate/2/max_f))
        N_w=odd_round_positive(sample_rate/min_f/q(np.arange(1,max_p+1))).astype('int')
        W=np.zeros((max_p,N_w.max()),dtype='complex128')
        center=N_w.max() // 2
        self.min_v=min_f/sample_rate
        for r,n_w in enumerate(N_w):
            nr=np.arange(-(n_w//2),n_w//2+1)
            c=center+nr
            W[r,c] = window(n_w) * np.exp(-j*2.*np.pi*nr*self.min_v*(r+1))
        self._w=W.mean(axis=0)
        # TODO: not sure how complex values are interpolated so splitting into
        # real and imaginary.
        self.w=np.vstack((self._w.real,self._w.imag))
        self.n=np.arange(N_w.max()) - center
        self.w_lu=interpolate.interp1d(self.n,self.w,bounds_error=False,fill_value=0)
        self.N_w_max=N_w.max()
        self.max_p=max_p
        
    def dft_v(self,x,v0,p=0,N=1):
        """
        Compute the inner product of x with self.w. self.w is over- or
        undersampled so as to adjust its fundamental (in normalized frequency)
        to v0.
        x should have length == len(self.n)
        p (>= 0) is the order of the derivative of the fourier transform
        computed.
        N is a scalar passed to dk_scale, so it is used only when p is > 0 (when
        derivatives are being computed). Basically it can be used to scale the
        resulting gradient, because if we were differentiating the fourier
        transform w.r.t. the bin number, the length of the fourier transform has
        to get divided out of the result. In this case N would be set to the
        length of the fourier transform. So len(x) is a reasonable value for N
        but it can be set as needed.
        """
        n_steps=v0/self.min_v
        n=np.multiply.outer(n_steps,self.n)
        w_=self.w_lu(n)
        # convert back to complex
        w=w_[0] + j*w_[1]
        ret=w@(dk_scale(N,p)*np.power(self.n,p)*x)
        return ret
        
    def dft_bin_log_ps_dk(self,x,v,N=1):
        X=self.dft_v(x,v,N=N)
        dX=self.dft_v(x,v,p=1,N=N)
        return log_ps_dk(X,dX)
        

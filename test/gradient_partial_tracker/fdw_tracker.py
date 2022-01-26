# A tracker using freq_dom_window

from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type
import numpy as np
import matplotlib.pyplot as plt
from dftdk import gradient_ascent_step_harm_lock
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal

def sq_mod(x):
    return np.real(x*np.conj(x))

def sq_mod_dv(X,dX):
    return np.real(X*np.conj(dX)+np.conj(X)*dX)

def log_sq_mod_dv(X,dX):
    return np.real(np.conj(dX)/np.conj(X)+dX/X)

class fdw_tracker:

    def __init__(self,
        # frame size
        N=2048,
        # window size
        N_w=2048,
        # hop size
        N_h=2048//4,
        # oversampling of DFT when computing the window tables in the frequency
        # domain
        dft_os=32,
        # coefficients of a sum of cosine window
        win_coeffs=[0.35875, 0.48829, 0.14128, 0.01168], # blackmanharris
        # the number of bins to one side of the peak that are included in inner
        # products in the time-domain for sparse matrix operations approximating
        # fourier transforms using the window described by win_coeffs
        win_lobe_radius=6,
    ):

        self.N=N
        self.N_w=N_w
        self.N_h=N_h
        self.win_type=sum_of_cos_dft_win_type(win_coeffs,win_lobe_radius)
        self.fdw=freq_dom_window(N,
            self.win_type,
            dft_os)

    def analyse(self,
        # signal to analyse
        x,
        # starting frequency
        vstart,
        # the number of hops where no gradient step is made (to ignore noisy
        # gradients caused by transients)
        warm_up_hops=4,
        # the number of harmonics to track
        n_harms=10,
        # gradient step coefficient
        mu=0.5,
        # number of gradient steps per frame
        n_steps=1,
    ):
        """
        Returns tuple of frequency estimates and sample times at which they
        occur.
        """
        
        N_x = len(x)
        
        h=np.arange(0,N_x-self.N,self.N_h,dtype='int')


        def grad(i):
            def _grad(x,k0):
                if i < warm_up_hops:
                    return np.zeros_like(k0)
                v0=k0/self.N
                X_dft_R=self.fdw.R(v0)
                X_dft_fdw=dft(x,X_dft_R)
                # derivative
                X_dft_dv_fdw=dft_dv(x,X_dft_R)
                return log_sq_mod_dv(X_dft_fdw,X_dft_dv_fdw)
            return _grad

        k0=vstart*self.N*(1+np.arange(n_harms))
        buf=np.zeros(self.N,dtype=x.dtype)
        k=np.zeros((len(k0),len(h)+1))
        k[:,0] = k0
        for i, n_h in enumerate(h):
            buf[:]=x[n_h:n_h+self.N]
            for m in range(n_steps):
                k0=gradient_ascent_step_harm_lock(buf,k0,mu,grad=grad(i))
            k[:,i+1]=k0

        return (k,h)

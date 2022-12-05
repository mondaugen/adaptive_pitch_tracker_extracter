# A tracker using freq_dom_window

from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type, dft_frame
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

class v_dev_stop_criterion:
    """ Allows stopping if v changes by more than a set amount """
    def __init__(self,initial_v,v_select=0,max_dev_semitones=1):
        """
        initial_v is the initial v value whose dimensionality should match v (in
        the __call__ method) after v_select has been used to select the subset.
        For example, if v is [1,2,3] and v_select is any integer, then initial_v
        should be [1] (say).
        v_select selects a subset of v:
            if v_select is an integer, this index of v is chosen to compare
            other values of v_select are not yet supported
        max_dev_semitones is the maximum deviation in semitones
        """
        self.v=initial_v
        self.v_select=v_select
        self.max_dev_semitones=max_dev_semitones
    def __call__(self,h,v,g,X):
        abs_log_dev_ratio=np.abs(np.log(v[self.v_select]/self.v)/np.log(2))*12.
        ret = False
        if abs_log_dev_ratio > self.max_dev_semitones:
            ret = True
        self.v=v[self.v_select]
        return ret

class soc_fdw_lookup:
    """
    Initialize sum-of-cosine "frequency-domain windows" so they can be looked up
    easily.
    """
    def __init__(self,
        # frame size
        N=2048,
        # window size
        N_w=2048,
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
        self.win_type=sum_of_cos_dft_win_type(win_coeffs,win_lobe_radius)
        self.fdw=freq_dom_window(N,
            self.win_type,
            dft_os)


class fdw_tracker(soc_fdw_lookup):
    # TODO: fdw_tracker with the analyse method should be a subclass of a simple
    # class that just comprises __init__ so other stuff can use this
    # initialization mechanism.

    def analyse(self,
        # signal to analyse
        x,
        # starting frequencies
        vstarts,
        # hops "iterator" (must implement __len__ method if no stopping
        # criterion provided)
        h=None,
        # Hop size (deprecated)
        N_h=None,
        # frequency groups
        v_groups=None,
        # the number of hops where no gradient step is made (to ignore noisy
        # gradients caused by transients)
        warm_up_hops=4,
        # gradient step coefficient
        mu=0.5,
        # number of gradient steps per frame
        n_steps=1,
        # gradient weighting argument
        grad_weight='equal',
        # a function passed a tuple of (h, v, gradients, X) where
        # h is the current hop
        # v are the current normalized frequencies of the partial tracks,
        # gradients are the gradients evaluated at those frequencies and X the
        # DFT evaluated at those frequencies)
        # The function should return True if the tracking should stop or True to
        # continue
        stop_criterion=lambda h,v,g,X: False
    ):
        """
        Returns tuple of frequency estimates and sample times at which they
        occur.
        """
        

        if h is None:
            N_x = len(x)
            if N_h is None:
                N_h = self.N//4
            
            if N_h < 0:
                h_start=N_x-self.N
                h_end=0
            elif N_h > 0:
                h_start=0
                h_end=N_x-self.N
            else:
                raise ValueError(f"N_h={N_h} but must be non-zero.")

            h=np.arange(h_start,h_end,N_h,dtype='int')

        grads=[]

        def store_grad(func):
            def fun(*args):
                ret=func(*args)
                grads.append(ret)
                return ret
            return fun

        def grad(i):
            def _grad(x,k0):
                if i < warm_up_hops:
                    return np.zeros_like(k0)
                v0=k0/self.N
                X_dft_R=self.fdw.R(v0)
                X_dft_fdw=x.dft(X_dft_R)
                # derivative
                X_dft_dv_fdw=x.dft_dv(X_dft_R)
                ret=log_sq_mod_dv(X_dft_fdw,X_dft_dv_fdw)
                grads.append(ret)
                return ret
            return store_grad(_grad)

        def k_to_v(k):
            return k/self.N
            
        def bin_vals(x,k0):
            v0=k_to_v(k0)
            X_dft_R=self.fdw.R(v0)
            X_dft_fdw=x.dft(X_dft_R)
            return X_dft_fdw
            
        k0=vstarts*self.N
        buf=np.zeros(self.N,dtype=x.dtype)
        k=np.zeros((len(k0),len(h)+1))
        X=np.zeros((len(k0),len(h)+1),dtype='complex128')
        k_list=[k0]
        X_list=[]
        for i, n_h in enumerate(h):
            buf[:]=x[n_h:n_h+self.N]
            buf_dft_frame=dft_frame(buf)
            X_list.append(bin_vals(buf_dft_frame,k0))
            for m in range(n_steps):
                k0=gradient_ascent_step_harm_lock(buf_dft_frame,
                k0,mu,grad=grad(i),
                groups=v_groups,
                grad_weight=grad_weight)
            k_list.append(k0)
            v0=k_to_v(k0)
            if stop_criterion(n_h,v0,grads[-1],X[-1]):
                break
        X_list.append(bin_vals(buf_dft_frame,k0))
        min_list_len=min([len(l) for l in [X_list,k_list,grads,h]])
        k=np.hstack([k_[:,None] for k_ in k_list[:min_list_len]])
        X=np.hstack([X_[:,None] for X_ in X_list[:min_list_len]])
        grads_ret=np.hstack([grad_[:,None] for grad_ in grads[:min_list_len]])
        h=h[:min_list_len]
        return (k,h,X,grads_ret)

def multi_analyse_combine(fdwt,**args):
    hs=args.pop('h')
    stop_criterion=args.pop('stop_criterion')
    ks=[]
    Xs=[]
    grads=[]
    final_hs=[]
    for h in hs:
        args['h']=h
        args['stop_criterion']=stop_criterion()
        k,h_,X,g=fdwt.analyse(**args)
        sort_i=np.argsort(h_)
        ks.append(k[:,sort_i])
        Xs.append(X[:,sort_i])
        grads.append(g[:,sort_i])
        final_hs.append(h_[sort_i])
        if h[1] > h[0]:
            print('final_hs',h_)
    hs_ret=np.concatenate(final_hs)
    ks_ret=np.hstack(ks)
    Xs_ret=np.hstack(Xs)
    grads_ret=np.hstack(grads)
    return ks_ret,hs_ret,Xs_ret,grads_ret

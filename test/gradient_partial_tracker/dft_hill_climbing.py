import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from exp_approx_line import exp_approx_line_coeffs
from recursive_dft import rec_dv_dft_sumcos

j=complex('j')

def dft_bin(x,v):
    N=len(x)
    v_shape=np.shape(v)
    if len(v_shape) > 0:
        # TODO Maybe not safe to use if len(v_shape) > 1 ?
        return np.exp(-j*2*np.pi*v[:,None]*np.arange(N))@x
    return np.sum(np.exp(-j*2*np.pi*v*np.arange(N))*x)

def dft_bin_dv(x,v):
    N=len(x)
    x_=x*-j*2*np.pi*np.arange(N)
    return dft_bin(x_,v)

class dft_bin_dv_approx:
    """
    Approximate the derivative using the product of the signal and an
    exponential that approximates the ramp. This is to see how well the
    derivative could work in a recursive implementation.
    """
    def __init__(self,N,err_max):
        """
        N is 1 more than the maximum value of the ramp.
        err_max is roughly the maximum error tolerated.
        """
        self.N=N
        alpha, beta, gamma, s0, s1 = exp_approx_line_coeffs(N,N,err_max)
        self.alpha=alpha
        self.beta=beta
        self.gamma=gamma
        self.s0=s0
        self.s1=s1
        self.approx_ramp=-j*2*np.pi*(np.exp(self.alpha*np.arange(self.N)+self.beta)
                            +self.gamma)
    def __call__(self,x,v):
        """ x must have length N """
        x_=x*self.approx_ramp
        return dft_bin(x_,v)
    
def dft_bin_d2v(x,v):
    N=len(x)
    x_=x*(-j*2*np.pi*np.arange(N))**2
    return dft_bin(x_,v)

def dft_bin_pow(x,v):
    X=dft_bin(x,v)
    return np.real(X*np.conj(X))

def dft_bin_pow_dv(x,v):
    X=dft_bin(x,v)
    dX=dft_bin_dv(x,v)
    return np.real(X*np.conj(dX)+np.conj(X)*dX)

def dft_bin_log_pow_dv(x,v,thresh_to_zero=1e-4,exp_dv=False):
    # thresh_to_zero forces the derivative to be zero if the magnitude of the
    # DFT at bin centred on v is too small. This is to avoid wild gradients
    # due to low-powered noise which should actually be considered silence.
    X=dft_bin(x,v)
    dX=dft_bin_dv(x,v)
    ret = np.real(X*np.conj(dX)+np.conj(X)*dX)/np.real(X*np.conj(X))
    if len(np.shape(X)) > 0:
        ret[np.abs(X) < thresh_to_zero] = 0
    else:
        if np.abs(X) < thresh_to_zero:
            return 0
    return ret

class dft_bin_log_pow_exp_dv:
    # thresh_to_zero forces the derivative to be zero if the magnitude of the
    # DFT at bin centred on v is too small. This is to avoid wild gradients
    # due to low-powered noise which should actually be considered silence.
    # Like dft_bin_log_pow_dv but the ramp multiplying the signal to estimate the
    # power spectrum derivative is generated with an exponential and so is
    # approximate (but faster to compute for recursive implementations)
    def __init__(self,N,err_max):
        self.dft_bin_dv = dft_bin_dv_approx(N,err_max)
    def __call__(self,x,v,thresh_to_zero=1e-4):
        X=dft_bin(x,v)
        dX=self.dft_bin_dv(x,v)
        ret = np.real(X*np.conj(dX)+np.conj(X)*dX)/np.real(X*np.conj(X))
        if len(np.shape(X)) > 0:
            ret[np.abs(X) < thresh_to_zero] = 0
        else:
            if np.abs(X) < thresh_to_zero:
                return 0
        return ret

def xs_dvxs_log_pow_dv(X,dX,thresh_to_zero=1e-4):
    """ Given the fourier transform and its derivative w.r.t. frequency, return the derivative of the log-power spectrum w.r.t frequency.
    """
    ret = np.real(X*np.conj(dX)+np.conj(X)*dX)/np.real(X*np.conj(X))
    ret[np.abs(X) < thresh_to_zero] = 0
    return ret

def dft_bin_pow_d2v(x,v):
    X=dft_bin(x,v)
    dX=dft_bin_dv(x,v)
    d2X=dft_bin_d2v(x,v)
    return np.real(X*np.conj(d2X)+dX*np.conj(dX)+np.conj(X)*d2X+np.conj(dX)*dX)

def get_padded_window(name,N,os=1):
    # Like get_window, but oversample by factor os and zero pad the window
    w_=signal.get_window(name,N)
    N_=N*fos
    zp=(N_-N)//2
    w=np.zeros(N*fos)
    w[zp:-zp]=w_
    return w

def dft_hill_climb(x,v_k,K=10,mu=0.01,newton=False):
    # works not great ?
    v_ks = np.zeros(K)
    pows = np.zeros(K)
    for k in range(K):
        if newton:
            v_k = v_k + mu * dft_bin_pow_dv(x,v_k) / dft_bin_pow_d2v(x,v_k)
        else:
            v_k = v_k + mu * dft_bin_pow_dv(x,v_k)
        v_ks[k] = v_k
        pows[k] = dft_bin_pow(x,v_k)
    return v_ks, pows

def adaptive_ghc_slow(x,v_k,w,mu=1e-6,max_step=1/16000):
    # Slow version of adaptive peak tracking using gradient
    # Max step is specified in normalized frequency (default 1 Hz at 16KHz Fs)
    buf=np.zeros_like(w)
    buf[1:]=x[:len(buf)-1]
    v_ks = np.zeros_like(x[len(buf)-1:])
    Xs = np.zeros_like(x[len(buf)-1:],dtype='complex128')
    grad = np.zeros_like(x[len(buf)-1:])
    for n, xn in enumerate(x[len(buf)-1:]):
        buf[:-1] = buf[1:]
        buf[-1] = xn
        grad[n]=dft_bin_pow_dv(buf*w,v_k)
        v_k = v_k + np.clip(mu * grad[n],-max_step,max_step)
        v_ks[n] = v_k
        Xs[n] = dft_bin(buf*w,v_k)
    return v_ks, Xs, grad

def adaptive_ghc_slow_log_pow(x,v_k,w,mu=1e-6,max_step=float('inf')):
    # Slow version of adaptive peak tracking using gradient of log power-spectrum
    # Max step is specified in normalized frequency (default 1 Hz at 16KHz Fs)
    # TODO: Sometimes the tracking goes really astray because the partial seems
    # to disappear at times (most likely due to resonance with a closely tuned
    # node). A solution might be to restrict the absolute deviation from the
    # initial guess.
    buf=np.zeros_like(w)
    buf[1:]=x[:len(buf)-1]
    v_ks = np.zeros_like(x[len(buf)-1:])
    Xs = np.zeros_like(x[len(buf)-1:],dtype='complex128')
    grad = np.zeros_like(x[len(buf)-1:])
    for n, xn in enumerate(x[len(buf)-1:]):
        buf[:-1] = buf[1:]
        buf[-1] = xn
        grad[n]=dft_bin_log_pow_dv(buf*w,v_k)
        v_k = v_k + np.clip(mu * grad[n],-max_step,max_step)
        v_ks[n] = v_k
        Xs[n] = dft_bin(buf*w,v_k)
    return v_ks, Xs, grad

def adaptive_ghc_slow_log_pow_v(x,v_k,w,mu=1e-6,max_step=float('inf'),
                                exp_dv=False,err_max=1e-4,verbose=False):
    # The same as adaptive_ghc_slow_log_pow but v_k can be a vector
    # if exp_dv is true, use ramp approximated by exponential to estimate derivative
    # err_max is a parameter for computing this ramp approximation
    buf=np.zeros_like(w)
    buf[1:]=x[:len(buf)-1]
    N_v=len(v_k)
    N_ret=len(x[len(buf)-1:])
    v_ks = np.zeros((N_ret,N_v))
    Xs = np.zeros((N_ret,N_v),dtype='complex128')
    grad = np.zeros((N_ret,N_v))
    grad_compute=dft_bin_log_pow_dv
    if verbose:
        print("exp_dv:",exp_dv)
    if exp_dv:
        if verbose:
            print("Using ramp-approximating exponential, err_max=%e."%(err_max,))
        grad_compute=dft_bin_log_pow_exp_dv(len(buf),err_max)
    for n, xn in enumerate(x[len(buf)-1:]):
        buf[:-1] = buf[1:]
        buf[-1] = xn
        cur_grad=grad_compute(buf*w,v_k)
        grad[n]=cur_grad
        v_k = v_k + np.clip(mu * grad[n],-max_step,max_step)
        v_ks[n] = v_k
        Xs[n] = dft_bin(buf*w,v_k)
    return v_ks, Xs, grad

def adaptive_ghc_recsumcos_log_pow_v(x,v_k,wp,Nw,mu=1e-6,max_step=float('inf'),
                                     ramp_err_max=1e-4,grad_step_warmup=True):
    """
    Tracks partials using a recursively updated adaptive gradient algorithm.
    x is the signal whose partials it tracks.
    v_k is a vector containing the initial (normalized) frequency estimates (0.5
    is the half the sample rate).
    wp is a vector containing the sum-of-cosines coefficients. For example, for
    the Hann window, you would provide [0.5,-0.5].
    Nw is the length of the sum-of-cosines window.
    mu is the gradient scalar for the gradient ascent algorithm.
    max_step is the maximum absolute step size.
    ramp_err_max is qualitatively the maximum error tolerated for the
    ramp-approximating exponential. In practice it actually depends on Nw.
    grad_step_warmup, if true, means that the v_k will not be updated using the
    gradient step until Nw samples have been processed. This is to prevent using
    gradients that occur due to the transition from silence to signal for early
    steps. In this case, the length of Xv and dvXv will be len(x) - Nw + 1.
    """
    N_v=len(v_k)
    xv_dvxv_computer=rec_dv_dft_sumcos(Nw,wp,N_v,max_err=ramp_err_max)
    N_ret=len(x)
    if grad_step_warmup:
        for n, xn in enumerate(x[:Nw-1]):
            _,_=xv_dvxv_computer.update(xn,v_k)
            N_ret -= 1
    v_ks = np.zeros((N_ret,N_v))
    Xs = np.zeros((N_ret,N_v),dtype='complex128')
    grad = np.zeros((N_ret,N_v))
    for n, xn in enumerate(x[-N_ret:]):
        Xv,dvXv=xv_dvxv_computer.update(xn,v_k)
        cur_grad=xs_dvxs_log_pow_dv(Xv,dvXv)
        grad[n]=cur_grad
        v_k = v_k + np.clip(mu * grad[n],-max_step,max_step)
        v_ks[n] = v_k
        Xs[n] = Xv
    return v_ks, Xs, grad

if __name__ == '__main__':
    N=256
    Fs=16e3
    f0=1000
    f_k=1025
    v0=f0/Fs
    v1=1100/Fs
    v_k=f_k/Fs
    # frequency oversampling factor
    fos=32
    n=np.arange(0,N,1/fos)
    v=n/N
    w=get_padded_window('hann',N,fos)
    w/=np.sum(w)
    x_=np.exp(j*2*np.pi*v0*n*fos)+np.random.standard_normal(N*fos)*3
    x_+=np.exp(j*2*np.pi*v1*n*fos)
    x=x_*w
    X=np.abs(np.fft.fft(x))**2
    X_k=dft_bin_pow(x,v_k)
    print('v_k',v_k)
    print('X_k',X_k)
    # Gradient ascent works not bad. This update would occur every sample, remember.
    # However, why does mu have to be so small? We need to scale a variable somewhere.
    v_ks,X_ks=dft_hill_climb(x,v_k,K=30,mu=1e-6)
    # TODO: newton doesn't work that well...
    #v_ks,X_ks=dft_hill_climb(x,v_k,K=30,mu=0.03,newton=True)

    plt.plot(v,X)
    plt.plot(v_k,X_k,'.',label='start')
    plt.plot(v_ks,X_ks)
    plt.plot(v_ks[-1],X_ks[-1],'.',label='end')
    plt.legend()
    plt.show()

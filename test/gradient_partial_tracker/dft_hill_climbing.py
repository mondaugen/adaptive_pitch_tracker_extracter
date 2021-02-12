import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

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

def dft_bin_log_pow_dv(x,v,thresh_to_zero=1e-4):
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

def adaptive_ghc_slow_log_pow_v(x,v_k,w,mu=1e-6,max_step=float('inf')):
    # The same as adaptive_ghc_slow_log_pow but v_k can be a vector
    buf=np.zeros_like(w)
    buf[1:]=x[:len(buf)-1]
    N_v=len(v_k)
    N_ret=len(x[len(buf)-1:])
    v_ks = np.zeros((N_ret,N_v))
    Xs = np.zeros((N_ret,N_v),dtype='complex128')
    grad = np.zeros((N_ret,N_v))
    for n, xn in enumerate(x[len(buf)-1:]):
        buf[:-1] = buf[1:]
        buf[-1] = xn
        cur_grad=dft_bin_log_pow_dv(buf*w,v_k)
        grad[n]=cur_grad
        v_k = v_k + np.clip(mu * grad[n],-max_step,max_step)
        v_ks[n] = v_k
        Xs[n] = dft_bin(buf*w,v_k)
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

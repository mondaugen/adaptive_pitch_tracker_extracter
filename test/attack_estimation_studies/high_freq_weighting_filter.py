# try to design a filter that approximates a triangular frequency response

import numpy as np
from scipy import signal
import common
import matplotlib.pyplot as plt

def est_autocorr(x):
    """ Estimate autocorrelation by inverse transforming the powerspectrum """
    X=np.fft.fft(x)
    S=X*np.conj(X)
    r=np.fft.ifft(X)
    return r

def periodic_sinc(N,a):
    t_eval=np.concatenate((np.arange(N/2),np.arange(1,N/2+1)[::-1]*-1))
    return np.sinc(t_eval*a)

def a_to_r(a):
    """ Convert filter coefficients to reflection coefficients. This assumes you
    include the coefficient for y[n], where it divides all the coefficients by
    it (so the first coefficient is 1), then discards the y[n] coefficient
    (because it knows it is 1). """
    a=a.copy()
    a=a/a[0]
    a=a[1:]
    p=len(a)
    r=np.zeros_like(a)
    r[-1]=a[-1]
    for j in np.arange(p-1)[::-1]:
        s=1/(1-r[j+1]*np.conj(r[j+1]))
        a[:j+1]=s*(a[:j+1]-r[j+1]*np.conj(a[:j+1][::-1]))
        r[j]=a[j]
    return r

# From test/attack_estimation_studies/tri_freq_response.py it seems sinc^2 is
# its own autocorrelation, so we use this for the filter design number of poles
P=common.get_env('P',default=2,conv=int)
acc_mult=common.get_env('acc_mult',default=1,conv=int)

def fit_allpole_triangular_lowpass(
    # number of poles
    P,
    # length multiplier (because a longer autocorrelation sequence might give a
    # better fit, but it seems like actually not)
    acc_mult=1):

    # It seems the feedforward coefficients other than the first one often get
    # forced to 0, so we just use an allpole model
    Q=0

    # length of autocorrelation sequence we need
    N=common.next_pow_2(acc_mult*2*(max(P,Q)+1))

    # autocorrelation, assumed indices [0,N/2] are the points [0,N/2] and
    # [N/2+1,N-1] are [-N/2,-1]
    r=0.5*np.power(periodic_sinc(N,0.5),2)

    # build autocorrelation matrix
    # get indices
    k=np.arange(1,P+1)
    k_l=np.add.outer(k,-1*np.arange(1,P+1))
    # wrap so they are in [0,N)
    while np.any(k_l < 0):
        k_l[k_l < 0] += N
    R=r[k_l]
    r_x=r[k][:,None]

    # compute feedback coefficients
    a=np.linalg.solve(R,-r_x)

    # compute feedforward coefficients
    # in this case r is equivalent to x
    n=np.arange(0,Q+1)
    n_k=np.add.outer(n,-k)
    while np.any(n_k < 0):
        n_k[n_k < 0] += N
    X=r[n_k]
    x=r[n]
    b=x[:,None]+X@a

    a0=a[0]
    a=a/a0
    b=b/a0
    return b,a

def fit_allpole_triangular_highpass(
    # number of poles
    P,
    # length multiplier (because a longer autocorrelation sequence might give a
    # better fit)
    acc_mult=1):
    b,a=fit_allpole_triangular_lowpass(P,acc_mult=acc_mult)
    a*=np.power(-1,np.arange(len(a)))[:,None]
    return b,a

if __name__ == '__main__':

    b,a=fit_allpole_triangular_lowpass(P,acc_mult=acc_mult)

    print('a')
    print(a)
    if np.any(np.abs(a)>1):
        print('unstable')
    else:
        print('stable')
    print('b')
    print(b)
    r=a_to_r(a)
    print('r')
    print(r)
    # to check if a_to_r is correct
    # the reflection coefficients should be: [.5, .2, -.5]
    # see p. 237 of Statistical Signal Processing and Modeling by Monson
    print(a_to_r(np.array([1,.5,-0.1,-0.5])))

    # show frequency response
    w,h=signal.freqz(b,a)
    plt.plot(w,np.abs(h),label='low-pass')

    # show high-pass configuration
    b,a=fit_allpole_triangular_highpass(P,acc_mult=acc_mult)
    w,h=signal.freqz(b,a)
    plt.plot(w,np.abs(h),label='high-pass')
    plt.legend()

    plt.show()

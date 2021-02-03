import numpy as np
import matplotlib.pyplot as plt

def R(x):
    """
    Returns a sequence r where r[k] is the estimated autocorrelation at lag
    k for assumed wide-sense stationary process x.
    """
    N=len(x)
    d=1/N
    r=np.zeros_like(x)
    for k in np.arange(N):
        r[k]=np.sum(x[k:]*np.conj(x[:N-k]))*d
    return r

def R_dft(x):
    """
    The same as R but uses the DFT in the computation.
    """
    N=len(x)
    d=1/N
    x_=np.zeros(2*N)
    x_[:N]=x
    X=np.fft.fft(x_)
    X_=X*np.conj(X)
    r=np.fft.ifft(X_)*d
    return r[:N]

N=32
x=np.random.standard_normal(N)
r=R(x)
r_dft=R_dft(x)
n=np.arange(N)
plt.plot(n,r,label='R')
plt.plot(n,r_dft,label='R_dft')
plt.legend()
plt.show()

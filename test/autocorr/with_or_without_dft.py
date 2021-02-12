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

def R_dft(x,axis=-1):
    """
    The same as R but uses the DFT in the computation.
    """
    N=x.shape[axis]
    d=1/N
    X=np.fft.fft(x,2*N,axis=axis)
    X_=X*np.conj(X)
    r=np.fft.ifft(X_,axis=axis)*d
    ret=np.moveaxis(r,axis,0)[:N]
    ret=np.moveaxis(ret,0,axis)
    return ret

N=32
x=np.random.standard_normal(N)
r=R(x)
r_dft=R_dft(x)
n=np.arange(N)
plt.plot(n,r,label='R')
plt.plot(n,r_dft,label='R_dft')
plt.legend()

plt.figure()
D=10
x2=np.random.standard_normal((D,N))
r2=R_dft(x2)
print(r2.shape)
for row in r2:
    plt.plot(n,row)

plt.figure()
x3=x2.T
r3=R_dft(x3,axis=0)
print(r3.shape)
for col in r3.T:
    plt.plot(n,col)

plt.show()

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import common

# This should give a Fourier transform with a frequency spectrum that is
# entirely real, starts at 1 and goes to 0 at 0.5 the sample rate

def est_autocorr(x):
    """ Estimate autocorrelation by inverse transforming the powerspectrum """
    X=np.fft.fft(x)
    S=X*np.conj(X)
    r=np.fft.ifft(X)
    return r

N=common.get_env('N',default=512,conv=int)
n=np.arange(N)
t=np.arange(-N/2+0.5,N/2)
t_eval=np.concatenate((np.arange(N/2),np.arange(1,N/2+1)[::-1]*-1))
a=common.get_env('a',default=0.5,conv=float)
x=a*np.power(np.sinc(t_eval*a),2)
X=np.fft.fft(x)
r=est_autocorr(x)
fig,axs=plt.subplots(3,1)
axs[0].plot(n,x)
axs[1].plot(n,np.real(X))
axs[2].plot(n,np.abs(r))
# It really seems that the square of the sinc is its own autocorrelation
# function
print(np.sum(np.abs(x-r)))

plt.show()

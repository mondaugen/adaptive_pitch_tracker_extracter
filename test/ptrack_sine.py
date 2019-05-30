import ptracker
import numpy as np
import matplotlib.pyplot as plt

w=0.4
w0=1.5
t=np.arange(0,100000)
x=np.exp(t*w*complex("j"))+np.random.standard_normal(t.shape)*.3
apt=ptracker.apt(
0.9,
0.9,
0.9,
2)
a,phi,corr,err=apt.proc(x,w0)

y=a*np.cos(phi)

frq=np.diff(np.unwrap(phi))

fig,ax=plt.subplots(6,1)

ax[0].specgram(x, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')
ax[1].specgram(y, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')

ax[2].plot(t,corr)

ax[3].specgram(err, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')

ax[4].plot(t[:-1],frq)
ax[5].plot(t,a)

plt.show()



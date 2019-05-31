from ptracking import ptrackers
import numpy as np
import matplotlib.pyplot as plt

w=0.4
w0=.5
t=np.arange(0,10000)
x=np.cos(t*w)+np.random.standard_normal(t.shape)*.1

print(dir(ptrackers))

a,phi=ptrackers.ptracker(
x,
w0,
0.99,
0.9,
.1,
1e-3)

y=a*np.cos(phi)

frq=np.diff(np.unwrap(phi))

fig,ax=plt.subplots(4,1)

ax[0].specgram(x, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')
ax[1].specgram(y, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')

ax[2].plot(t[:-1],frq)
ax[3].plot(t,a)

plt.show()




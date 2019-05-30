import ptracker
import numpy as np
import matplotlib.pyplot as plt
import librosa
import attack_finder

sr=16000
x,sr=librosa.load("sounds/midi.wav",sr)
t=np.arange(len(x))/sr
w0=130.81/sr*2*np.pi # frequency of c3
apt=ptracker.apt(
0.95,
0.9,
0.9,
1e-4,
g_thresh=1e-3)
a,phi=apt.proc2(x,w0)

y=a*np.cos(phi)

frq=np.diff(np.unwrap(phi))

fig,ax=plt.subplots(6,1)

ax[0].specgram(x, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')
ax[1].specgram(y, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')

#ax[2].plot(t,corr)
#
#ax[3].specgram(err, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
#               cmap=None, xextent=None, pad_to=None, sides='default',
#               scale_by_freq=None, mode='default', scale='default')
ax[3].plot(t,attack_finder.fast_max(a,alph=0.999))

ax[4].plot(t[:-1],frq)
ax[5].plot(t,a)

plt.show()



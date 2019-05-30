import matplotlib.pyplot as plt
import ptracker
import librosa
import numpy as np

x,sr=librosa.load('sounds/guitar.wav',sr=16000)

ptrackers=ptracker.pitch_check_comb_array()

y=ptrackers.proc(x)
y=y*y
y/=np.max(x)

#plt.imshow(np.log(y),aspect='auto',origin='lower')
plt.plot(np.arange(len(x))/sr,np.log(y[48,:]),label="c3")
plt.plot(np.arange(len(x))/sr,np.log(y[47,:]),label="b2")
plt.legend()
plt.figure()
plt.specgram(x, NFFT=2048, Fs=sr, Fc=0, noverlap=512,
             cmap=None, xextent=None, pad_to=None, sides='default',
             scale_by_freq=None, scale='dB')
plt.figure()
maxs=np.argmax(y[:100,:],axis=0)
plt.plot(np.arange(len(x))/sr,maxs,'o')

plt.figure()
plt.imshow(np.log(y),aspect='auto',origin='lower')

plt.show()

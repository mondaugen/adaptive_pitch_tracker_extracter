from ptracking import ptrackers
import numpy as np
import matplotlib.pyplot as plt
import librosa
import attack_finder

sr=16000
#x,sr=librosa.load("sounds/one_guitar_note.wav",sr)
#x,sr=librosa.load("/tmp/sine.wav",sr)
x=np.fromfile('/tmp/one_guitar_note.f64',dtype='float64')
sr=16000
#x=np.fromfile('/tmp/sine.f64',dtype='float64')
#librosa.output.write_wav("/tmp/one_guitar_note_librosa.wav",x,sr)
t=np.arange(len(x))/sr
w0=65.406/sr*2*np.pi # frequency of c2
#w0=56*130.81/sr*2*np.pi # frequency of c3
#w0=123.47/sr*2*np.pi # frequency of b2
#w0=130.81/sr*2*np.pi * np.power(2,-.5/12.) # frequency of c3 detuned by 50 cents

a,phi,w=ptrackers.ptracker(
x,
w0,
0.999,
0.9,
1e-2,
1e-3)
#1e-3)
#0.999,
#0.999,
#1e-3,
#1e-4)

y=a*np.cos(phi)

frq=np.diff(np.unwrap(phi))
#frq=np.diff(phi)

fig,ax=plt.subplots(5,1)

ax[0].specgram(x, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')
ax[1].specgram(y, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=4,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')

ax[2].plot(t,w)
ax[3].plot(t,a)


ax[4].plot(t,x)

plt.show()



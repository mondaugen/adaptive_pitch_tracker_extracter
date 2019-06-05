import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import librosa

sr=16000
x=np.fromfile('/tmp/snd.f64')

n_samps=min(len(x),10000)
x=x[:n_samps]

N_FFT=256
N_W=256
w=signal.get_window('hann',N_W)
H=1

x_frames=librosa.util.frame(x, frame_length=N_W, hop_length=H)
X=np.fft.rfft(x_frames,axis=0)
#ph=np.angle(X[:,1:]/X[:,:-1])
ph=np.diff(np.unwrap(np.angle(X)))

fig,ax=plt.subplots(2,1)

t=np.arange(X.shape[1])*H/sr
n_parts=32
for ph_ in ph[:n_parts,:]:
    ax[0].plot(t[:-1],ph_,'b')

ax[1].plot(t[:-1],x[:len(t[:-1])])

plt.show()


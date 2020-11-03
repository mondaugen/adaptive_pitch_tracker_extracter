# Plot the waveform as a timeseries on the left and the spectrogram on the right
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import os
from spectral_difference import local_max, local_max_mat

def local_max_2d(x,thresh=1e-5):
    x_c=local_max_mat(x,indices=False)
    x_r=local_max_mat(x.T,indices=False).T
    return np.where(x_c & x_r & (x > thresh))

INFILE=os.environ.get('INFILE','/tmp/in.f32')
FS=float(os.environ.get('FS','16000'))
NFFT=2048
x=np.fromfile(INFILE,'float32')

N=len(x)
n=np.arange(N)
t=n/FS

fig,ax=plt.subplots(1,2)

ax[0].plot(t,x)
S,f,t,_=ax[1].specgram(x,window=signal.get_window('blackmanharris',NFFT), NFFT=NFFT,noverlap=int(NFFT*.75),Fs=FS)
lmax_r,lmax_c=local_max_2d(S,thresh=1e-6)
ax[1].scatter(t[lmax_c],f[lmax_r],c='red',s=1)

fig.suptitle(INFILE)

plt.show()

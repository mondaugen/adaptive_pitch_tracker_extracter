# Plot the waveform as a timeseries on the left and the spectrogram on the right
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import os
from spectral_difference import local_max, local_max_mat, spectral_diff
import common

envget=os.environ.get

def local_max_2d(x,thresh=1e-5):
    x_c=local_max_mat(x,indices=False)
    x_r=local_max_mat(x.T,indices=False).T
    return np.where(x_c & x_r & (x > thresh))

INFILE=envget('INFILE','/tmp/in.f32')
FS=float(envget('FS','16000'))
NFFT=int(envget('NFFT','2048'))
SAMPTYPE=envget('SAMPTYPE','float32')
LMAX=int(envget('LMAX','0'))
YSCALE=envget('YSCALE','linear')
SSCALE=envget('SSCALE','dB')
TMIN=float(envget('TMIN','0'))
TMAX=float(envget('TMAX','3'))
FMAX=float(envget('FMAX','5000'))
FMIN=float(envget('FMIN','10'))
SD_H=int(envget('SD_H','256'))
SD_W=int(envget('SD_W','1024'))
SD_WIN=envget('SD_WIN','hann')
# The bin at which the spectral difference is scaled by half, falling exponentially
SD_BIN_WEIGHT=float(envget('SD_BIN_WEIGHT','inf'))

x=np.fromfile(INFILE,SAMPTYPE)

N=len(x)
n=np.arange(N)
t=n/FS

fig,ax=plt.subplots(3,1,sharex=True)

ax[0].plot(t,x)
S,f,t,_=ax[1].specgram(x,window=signal.get_window('blackmanharris',NFFT), NFFT=NFFT,noverlap=int(NFFT*.75),Fs=FS,scale=SSCALE)
if LMAX:
    lmax_r,lmax_c=local_max_2d(S,thresh=1e-6)
    ax[1].scatter(t[lmax_c],f[lmax_r],c='red',s=1)

xlim=(TMIN,TMAX)
ax[1].set_xlim(xlim)

ylim=(FMIN,FMAX)
ax[1].set_ylim(ylim)
ax[1].set_yscale(YSCALE)

bin_weights=None
if SD_BIN_WEIGHT != float('inf'):
    bin_weights=np.power(2,-np.arange(SD_W)/SD_BIN_WEIGHT)

spec_diff=spectral_diff(x,SD_H,SD_W,SD_WIN,bin_weights=bin_weights)
spec_diff_times=common.frame_times(len(x),SD_H,SD_W).astype('float64')
# Add half a window and half a hop to the frame times, this will place the frame times at the middle of the total time spanned by the 2 frames that are differenced
spec_diff_times+=(SD_H+SD_W)*0.5
# Convert to seconds
spec_diff_times/=FS
ax[2].plot(spec_diff_times[:-1],spec_diff)

fig.suptitle(INFILE)

plt.show()

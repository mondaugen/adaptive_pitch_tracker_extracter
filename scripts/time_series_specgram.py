# Plot the waveform as a timeseries on the left and the spectrogram on the right
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
from scipy import signal
import os
from spectral_difference import local_max, local_max_mat, spectral_diff
import common
from peak_finder import find_peaks

envget=os.environ.get

def local_max_2d(x,thresh=1e-5):
    x_c=local_max_mat(x,indices=False)
    x_r=local_max_mat(x.T,indices=False).T
    return np.where(x_c & x_r & (x > thresh))

def dB(x):
    return 20*np.log10(x)

INFILE=envget('INFILE','/tmp/in.f32')
FS=float(envget('FS','16000'))
NFFT=int(envget('NFFT','2048'))
HFFT=int(envget('HFFT',int(NFFT*0.25)))
SAMPTYPE=envget('SAMPTYPE','float32')
LMAX=int(envget('LMAX','0'))
YSCALE=envget('YSCALE','linear')
SSCALE=envget('SSCALE','dB')
TMIN=float(envget('TMIN','0'))
TMAX=float(envget('TMAX','3'))
FMAX=float(envget('FMAX','5000'))
FMIN=float(envget('FMIN','30'))
PLOT_SD=int(envget('PLOT_SD','0'))
SD_H=int(envget('SD_H','256'))
SD_W=int(envget('SD_W','1024'))
SD_WIN=envget('SD_WIN','hann')
# The bin at which the spectral difference is scaled by half, falling exponentially
SD_BIN_WEIGHT=float(envget('SD_BIN_WEIGHT','inf'))
# The amount a peak must exceed the nearest minimum to be considered a peak as a multiple of the minimum
PEAK_K=float(envget('PEAK_K','2'))
# An absolute value a peak must exceed to be considered a peak
PEAK_T=float(envget('PEAK_T','1e-5'))

x=np.fromfile(INFILE,SAMPTYPE)

N=len(x)
n=np.arange(N)
t=n/FS

fig,ax=plt.subplots(4,1)
sd_ax,sg_ax,slider_ax,psd_ax=ax

S,f,t_spec,_=sg_ax.specgram(x,window=signal.get_window('blackmanharris',NFFT), NFFT=NFFT,noverlap=NFFT-HFFT,Fs=FS,scale=SSCALE,mode='magnitude')
if LMAX:
    lmax_r,lmax_c=local_max_2d(S,thresh=1e-6)
    sg_ax.scatter(t_spec[lmax_c],f[lmax_r],c='red',s=1)


ylim=(FMIN,FMAX)
sg_ax.set_ylim(ylim)
sg_ax.set_yscale(YSCALE)

if PLOT_SD:
    bin_weights=None
    if SD_BIN_WEIGHT != float('inf'):
        bin_weights=np.power(2,-np.arange(SD_W)/SD_BIN_WEIGHT)

    spec_diff=spectral_diff(x,SD_H,SD_W,SD_WIN,bin_weights=bin_weights)
    spec_diff_times=common.frame_times(len(x),SD_H,SD_W).astype('float64')
    # Add half a window and half a hop to the frame times, this will place the frame times at the middle of the total time spanned by the 2 frames that are differenced
    spec_diff_times+=(SD_H+SD_W)*0.5
    # Convert to seconds
    spec_diff_times/=FS
    sd_ax.plot(spec_diff_times[:-1],spec_diff)
else:
    sd_ax.plot(t,x)

# Vertical line showing current frame shown in power-spectrum window
t_HFFT=HFFT/FS
tmp=sg_ax.get_ylim()
ylims=[tmp[0],tmp[0],tmp[1],tmp[1]]
# Store the line so we can move it with the slider
frame_line=sg_ax.axvline(0)

# Figure 2 shows the power spectral density of the selected frame
psd_ax.plot(f,dB(S[:,0]))
psd_ax.set_ylim([-100,dB(S.max())])
# If specgram has log-scaled y axis, then PSD has log-scaled x-axis
psd_ax.set_xlim(ylim)
psd_ax.set_xscale(YSCALE)

# Slider to choose a frame
sframe=Slider(slider_ax,'Frame (s)',0,t_spec.max(),valinit=0,valstep=t_HFFT)
def update_frame_line(val):
    next_val=val+t_HFFT
    frame_line.set_xdata([val,val])
    frame=int(round(val/t_HFFT))
    psd_ax.clear()
    psd_ax.plot(f,dB(S[:,frame]))
    peaks,max_over_min,max_i,min_i,ext_val_fhold,ext_val_bhold=find_peaks(S[:,frame],K=PEAK_K,
        T=PEAK_T)
    psd_ax.plot(f[peaks],dB(S[peaks,frame]),'.')
    psd_ax.set_ylim([-100,dB(S.max())])
    # If specgram has log-scaled y axis, then PSD has log-scaled x-axis
    psd_ax.set_xlim(ylim)
    psd_ax.set_xscale(YSCALE)
sframe.on_changed(update_frame_line)

xlim=(TMIN,TMAX)
for ax_ in [sd_ax,sg_ax,slider_ax]:
    ax_.set_xlim(xlim)

fig.suptitle(INFILE)

plt.show()

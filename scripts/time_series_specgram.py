# Plot the waveform as a timeseries on the left and the spectrogram on the right
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
from scipy import signal
import os
from spectral_difference import local_max, local_max_mat, spectral_diff
import common
from peak_finder import find_peaks
import dft_hill_climbing as dhc

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
# Enable/disable (1/0) partial tracking
PTRACK=int(envget('PTRACK','0'))
# Partial tracking starting time in seconds (will be rounded to nearest sample)
PTRACK_T0=float(envget('PTRACK_T0','0'))
# Partial tracking duration
PTRACK_DUR=float(envget('PTRACK_DUR','0.5'))
# Partial tracking starting frequency (will be divided by Fs to get normalized frequency)
PTRACK_F0=float(envget('PTRACK_F0','1000'))
# Partial tracking gradient ascent step size
PTRACK_MU=float(envget('PTRACK_MU','1e-8'))
# Partial tracking maximum step size in Hz
PTRACK_MAX_STEP=float(envget('PTRACK_MAX_STEP','1'))
# Partial tracking window type
PTRACK_WINTYPE=envget('PTRACK_WINTYPE','hann')
# Partial tracking window length
PTRACK_WINLEN=int(envget('PTRACK_WINLEN','4096'))
# If true, tries to remove partial from original sound and plots spectrogram
PTRACK_REMOVE=int(envget('PTRACK_REMOVE','0'))
# Where to write the partial to
PTRACK_OUT=envget('PTRACK_OUT','/tmp/ptrack_out.f64')
# Where to write the audio with the removed partial to
PTRACK_REMOVED_OUT=envget('PTRACK_REMOVED_OUT','/tmp/ptrack_removed_out.f64')


x=np.fromfile(INFILE,SAMPTYPE)

N=len(x)
n=np.arange(N)
t=n/FS

fig,ax=plt.subplots(7,1)
sd_ax,sg_ax,slider_ax,psd_ax,ptrack_grad_ax,ptrack_extracted_ax,ptrack_mag_ax=ax

ylim=(FMIN,FMAX)
def do_specgram(ax_,x_):
    ret = ax_.specgram(x_,
                window=signal.get_window('blackmanharris',NFFT),
                NFFT=NFFT,noverlap=NFFT-HFFT,Fs=FS,
                scale=SSCALE,mode='magnitude')
    ax_.set_ylim(ylim)
    ax_.set_yscale(YSCALE)
    return ret

S,f,t_spec,_=do_specgram(sg_ax,x)

if LMAX:
    lmax_r,lmax_c=local_max_2d(S,thresh=1e-6)
    sg_ax.scatter(t_spec[lmax_c],f[lmax_r],c='red',s=1)



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

if PTRACK:
    # Partial tracking
    # Track partial using supplied start and end times and starting frequency
    ptrack_n0=int(round(PTRACK_T0*FS))
    ptrack_n1=ptrack_n0+int(round(PTRACK_DUR*FS))
    ptrack_v0=PTRACK_F0/FS
    ptrack_w=signal.get_window(PTRACK_WINTYPE,PTRACK_WINLEN)
    # Normalize window
    ptrack_w/=np.sum(ptrack_w)
    # TODO: How to incorporate this time offset with the phase estimate in Xs?
    ptrack_t=(np.arange(ptrack_n0,ptrack_n1)+PTRACK_WINLEN*0.5)/FS
    # also add window's length of samples so we have an analysis for every t in
    # the specified range
    v_ks,Xs,grad=dhc.adaptive_ghc_slow_log_pow(x[ptrack_n0:ptrack_n1+PTRACK_WINLEN-1],ptrack_v0,
                                       ptrack_w,mu=PTRACK_MU,
                                       max_step=PTRACK_MAX_STEP/FS)
    sg_ax.plot([PTRACK_T0],[PTRACK_F0],'r.')
    sg_ax.plot(ptrack_t,v_ks*FS,lw=1)
    ptrack_grad_ax.plot(ptrack_t,np.log(np.abs(grad)))
    ptrack_grad_ax.set_title('log gradient')
    ptrack_extracted_ax.plot(ptrack_t,np.real(Xs))
    ptrack_extracted_ax.set_title('extracted partial')
    ptrack_mag_ax.plot(ptrack_t,np.abs(Xs))
    ptrack_mag_ax.set_title('extracted partial amplitude')
    Xs_out=np.real(Xs).copy()
    common.normalize(Xs_out).tofile(PTRACK_OUT)
    if PTRACK_REMOVE:
        fig_ptrack,ax_ptrack=plt.subplots(1,1)
        x_no_partial=np.zeros_like(x)
        x_no_partial[:] = x
        x_no_partial[ptrack_n0:ptrack_n1] -= 2*np.real(Xs)
        do_specgram(ax_ptrack,x_no_partial)
        ax_ptrack.set_ylim([PTRACK_F0-100,PTRACK_F0+100])
        ax_ptrack.set_xlim([PTRACK_T0,PTRACK_T0+PTRACK_DUR])
        common.normalize(np.real(x_no_partial)).tofile(PTRACK_REMOVED_OUT)


# Vertical line showing current frame shown in power-spectrum window
t_HFFT=HFFT/FS
tmp=sg_ax.get_ylim()
ylims=[tmp[0],tmp[0],tmp[1],tmp[1]]
# Store the line so we can move it with the slider
frame_line=sg_ax.axvline(0)

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
update_frame_line(0)
sframe.on_changed(update_frame_line)

xlim=(TMIN,TMAX)
for ax_ in [sd_ax,sg_ax,slider_ax]:
    ax_.set_xlim(xlim)

fig.suptitle(INFILE)

plt.show()

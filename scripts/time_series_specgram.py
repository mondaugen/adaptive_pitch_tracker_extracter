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

def get_window_coefficients(name):
    if name == 'hann':
        return np.array([0.5,-0.5])
    if name == 'blackmanharris':
        return np.array([0.35875,0.48829,0.14128,0.01168])
    else:
        raise ValueError("Unsupported window %s" % (name,))

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
# Choose the partial tracking method
PTRACK_METHOD=envget('PTRACK_METHOD','recursive')
# Partial tracking starting time in seconds (will be rounded to nearest sample)
PTRACK_T0=float(envget('PTRACK_T0','0'))
# Partial tracking duration
PTRACK_DUR=float(envget('PTRACK_DUR','0.5'))
# Partial tracking starting frequency (will be divided by Fs to get normalized frequency)
# Can also be the string 'auto' in which case the frequencies will be chosen
# using peak picking at time T0
PTRACK_F0=envget('PTRACK_F0','1000')
# Partial tracking gradient ascent step size
PTRACK_MU=float(envget('PTRACK_MU','1e-8'))
# Partial tracking maximum step size in Hz
PTRACK_MAX_STEP=float(envget('PTRACK_MAX_STEP','1'))
# Partial tracking window type
PTRACK_WINTYPE=envget('PTRACK_WINTYPE','hann')
# Partial tracking window length
PTRACK_WINLEN=int(envget('PTRACK_WINLEN','4096'))
# Partial tracking hop size (applies only to certain algorithms)
PTRACK_H=int(envget('PTRACK_H','1024'))
# If true, tries to remove partial from original sound and plots spectrogram
PTRACK_REMOVE=int(envget('PTRACK_REMOVE','0'))
# Where to write the partial to
PTRACK_OUT=envget('PTRACK_OUT','/tmp/ptrack_out.f64')
# Where to write the audio with the removed partial to
PTRACK_REMOVED_OUT=envget('PTRACK_REMOVED_OUT','/tmp/ptrack_removed_out.f64')
# Where to write the partial frequency to
PTRACK_F0_OUT=envget('PTRACK_F0_OUT','/tmp/ptrack_f0_out.f64')
# Where to write the partial amplitude to
PTRACK_A_OUT=envget('PTRACK_A_OUT','/tmp/ptrack_a_out.f64')
# If PTRACK_EXP_DV is truthy, this specifies the maximum error of the
# ramp-approximating exponential
PTRACK_EXP_DV_ERR_MAX=float(envget('PTRACK_EXP_DV_ERR_MAX','1e-4'))
PTRACK_A_LOG=int(envget('PTRACK_A_LOG','0'))

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
    if PTRACK_F0 == 'auto':
        # find frequencies by peak picking
        x_local=x[ptrack_n0:ptrack_n0+NFFT]
        dum_fig,dum_ax=plt.subplots()
        S_local,f_local,_,_=do_specgram(dum_ax,x_local)
        peaks,_,_,_,_,_=find_peaks(S_local[:,0],K=PEAK_K,T=PEAK_T)
        print("tracking %d peaks" % (len(peaks,)))
        ptrack_v0=f_local[peaks]/FS
    else:
        ptrack_v0=np.array([float(PTRACK_F0)/FS])
    x_ptrack=x[ptrack_n0:ptrack_n1+PTRACK_WINLEN-1]
    if PTRACK_METHOD == 'recursive':
        ptrack_wp=get_window_coefficients(PTRACK_WINTYPE)
        # Normalize window
        ptrack_wp/=(ptrack_wp[0]*PTRACK_WINLEN)
        # also add window's length of samples so we have an analysis for every t in
        # the specified range
        # TODO use vector of v0 to track multiple simultaneously
        v_ks,Xs,grad=dhc.adaptive_ghc_recsumcos_log_pow_v(x_ptrack,
                                           ptrack_v0,
                                           ptrack_wp,
                                           PTRACK_WINLEN,
                                           mu=PTRACK_MU,
                                           max_step=PTRACK_MAX_STEP/FS,
                                           ramp_err_max=PTRACK_EXP_DV_ERR_MAX)
        ptrack_t=(np.arange(ptrack_n0,ptrack_n1)+PTRACK_WINLEN*0.5)/FS
    elif PTRACK_METHOD == 'hop':
        print('Using "hop" method for partial tracking.')
        ptrack_w=signal.get_window(PTRACK_WINTYPE,PTRACK_WINLEN)
        # Normalize window
        ptrack_w/=np.sum(ptrack_w)
        v_ks,Xs,grad=dhc.adaptive_ghc_hop_log_pow_v(x_ptrack,ptrack_v0,ptrack_w,PTRACK_H,mu=PTRACK_MU,max_step=PTRACK_MAX_STEP/FS,verbose=False,harmonic_lock=False)
    # TODO: How to incorporate this time offset with the phase estimate in Xs?
    # Seems to not be necessary
        ptrack_t=(np.arange(ptrack_n0,ptrack_n1,PTRACK_H)+PTRACK_WINLEN*0.5)/FS
    sg_ax.plot(np.ones_like(ptrack_v0)*PTRACK_T0,
               ptrack_v0*FS,'r.')
    sg_ax.plot(ptrack_t,v_ks*FS,lw=1)
    ptrack_grad_ax.plot(ptrack_t,np.log(np.abs(grad)))
    ptrack_grad_ax.set_title('log gradient')
    ptrack_extracted_ax.plot(ptrack_t,np.real(Xs))
    ptrack_extracted_ax.set_title('extracted partial')
    if PTRACK_A_LOG:
        ptrack_mag_ax.plot(ptrack_t,20*np.log10(np.abs(Xs)))
    else:
        ptrack_mag_ax.plot(ptrack_t,np.abs(Xs))
    ptrack_mag_ax.set_title('extracted partial amplitude')
    if Xs.shape[1] == 1:
        # Right now you can't save multiple partial tracks
        # TODO make work for multiple tracks
        np.abs(Xs).tofile(PTRACK_A_OUT)
        v_ks.tofile(PTRACK_F0_OUT)
    Xs_out=np.sum(np.real(Xs).copy(),axis=1)
    common.normalize(Xs_out).tofile(PTRACK_OUT)
    if PTRACK_REMOVE:
        if PTRACK_METHOD == 'hop':
            raise NotImplementedError
        fig_ptrack,ax_ptrack=plt.subplots(1,1)
        x_no_partial=np.zeros_like(x)
        x_no_partial[:] = x
        x_no_partial[ptrack_n0:ptrack_n1] -= 2*np.sum(np.real(Xs),axis=1)
        do_specgram(ax_ptrack,x_no_partial)
        min_f0=ptrack_v0.min()
        max_f0=ptrack_v0.max()
        ax_ptrack.set_ylim([min_f0-100,max_f0+100])
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
fig.tight_layout()
plt.show()

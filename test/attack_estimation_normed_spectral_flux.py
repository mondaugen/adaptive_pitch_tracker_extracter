# Estimate attacks by computing the spectral flux
# This is simply the sum of the absolute values of the differences of magnitudes
# within 1 FFT frame
import numpy as np
from scipy import signal
import common
import matplotlib.pyplot as plt
import math

def spectral_flux(x,H,W,window_type='hann'):
    if W >= 4:
        w=signal.get_window(window_type,W)
    else:
        w=np.ones(W)
    w/=np.sum(W)
    spec_flux_i=0
    X=np.fft.fft(common.frame(x,H,W)*w[:,None],axis=0)
    spec_flux=np.sum(np.abs(np.diff(np.abs(X),axis=0)),axis=0)
    return spec_flux

INPUT=common.get_env('INPUT',check_if_none=True)
WINDOW_TYPE=common.get_env('WINDOW_TYPE',default='hann')
SAMPLE_RATE=common.get_env('SAMPLE_RATE',conv=float,default=16e3)
ATTACK_THRESH=common.get_env('ATTACK_THRESH',conv=float,default=0.05)
H=common.get_env('H',conv=int,default=2)
W=common.get_env('W',conv=int,default=8)
ALPHA=common.get_env('ALPHA',conv=float,default=.01)
# whether or not to normalize the spectral flux
NORMALIZE=common.get_env('NORMALIZE',conv=int,default=0)
# filter coefficient for smoothing the instantaneous amplitude
ALPHA_NORMALIZE=common.get_env('ALPHA_NORMALIZE',conv=float,default=0.01)
# threshold for the instantaneous amplitude below which no normalization is
# carried out (because very small amplitudes will cause huge values)
# value in dB
THRESH_NORMALIZE=common.get_env('THRESH_NORMALIZE',conv=float,default=-60)
X_LIM=common.get_env('X_LIM',conv=eval,default=(0,2))

x=np.fromfile(INPUT)

fig_signals,ax_signals=plt.subplots(5,1)
t=np.arange(len(x))/SAMPLE_RATE
ax_signals[0].plot(t,x)
ax_signals[0].set_title('original signal')
spec_flux=spectral_flux(x,H,W,window_type=WINDOW_TYPE)
n_spec_flux=(np.arange(len(spec_flux))*H+W/2).astype('int')
t_spec_flux=n_spec_flux/SAMPLE_RATE
if NORMALIZE:
    # normalize the spectral flux by a local average
    # the normalization is only applied if the instantaneous amplitude is
    # greater than THRESH_NORMALIZE
    th_norm_a=np.power(10,THRESH_NORMALIZE/20)
    inst_amp=np.abs(x[n_spec_flux])
    _a=ALPHA_NORMALIZE
    smooth_spec_flux_gain_scalar=1./(_a*(1+(1-_a)/(1+_a)))
    smooth_spec_flux_fc=([smooth_spec_flux_gain_scalar*_a],[1,-(1-_a)])
    smoothed_spec_flux=signal.lfilter(*smooth_spec_flux_fc,spec_flux)
    norm_scalars=np.zeros(smoothed_spec_flux.shape)
    inst_amp_th_mask=inst_amp>th_norm_a
    norm_scalars[inst_amp_th_mask]=1/smoothed_spec_flux[inst_amp_th_mask]
    spec_flux*=norm_scalars
    ax_signals[4].plot(t_spec_flux,np.abs(spec_flux))
    ax_signals[4].set_title('abs locally normalized spectral flux')
_a=ALPHA
smooth_flux_gain_scalar=1./(_a*(1+(1-_a)/(1+_a)))
smooth_filter_coeffs=([smooth_flux_gain_scalar*ALPHA],[1,-(1-ALPHA)])
spec_flux=signal.lfilter(*smooth_filter_coeffs,spec_flux)
ax_signals[3].plot(t_spec_flux,spec_flux)
ax_signals[3].set_title('smoothed spectral flux')
ax_signals[3].set_ylim(0,1)
spec_flux_diff=spec_flux[1:]-spec_flux[:-1]
attacks=spec_flux_diff.copy()
decay=spec_flux_diff.copy()
attacks[attacks<0]=0
decay[decay>0]=0
ax_signals[1].plot(t_spec_flux[1:],np.array(attacks>ATTACK_THRESH,dtype='float'))
ax_signals[1].set_title('estimated attacks')
#ax_signals[1].plot(t_spec_flux[1:],attacks)
ax_signals[2].plot(t_spec_flux[1:],decay)
ax_signals[2].set_title('estimated decays')
for ax in ax_signals:
    ax.set_xlim(*X_LIM)

plt.tight_layout()
plt.show()

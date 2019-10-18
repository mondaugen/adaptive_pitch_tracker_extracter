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
x=np.fromfile(INPUT)

filter_co_hz=2000
filter_order=16
filter_ripple=3
b,a=signal.cheby1(filter_order,filter_ripple,filter_co_hz,btype='highpass',fs=SAMPLE_RATE)
filter_coeffs=(b,a)

fig_mag_resp,ax_mag_resp=plt.subplots(1,1)
w,h=signal.freqz(*filter_coeffs)
ax_mag_resp.plot(w/(2*np.pi)*SAMPLE_RATE,20*np.log10(np.abs(h)))

x_filtered=signal.lfilter(*filter_coeffs,x)

fig_signals,ax_signals=plt.subplots(5,1)
t=np.arange(len(x))/SAMPLE_RATE
ax_signals[0].plot(t,x)
spec_flux=spectral_flux(x,H,W,window_type=WINDOW_TYPE)
t_spec_flux=(np.arange(len(spec_flux))*H+W/2)/SAMPLE_RATE
smooth_filter_coeffs=([ALPHA],[1,-(1-ALPHA)])
spec_flux=signal.lfilter(*smooth_filter_coeffs,spec_flux)
ax_signals[3].plot(t_spec_flux,spec_flux)
ax_signals[3].set_title('spec flux')
spec_flux_diff=spec_flux[1:]-spec_flux[:-1]
ax_signals[4].plot(t_spec_flux[1:],spec_flux_diff)
ax_signals[4].set_title('spec flux diff')
attacks=spec_flux_diff.copy()
decay=spec_flux_diff.copy()
attacks[attacks<0]=0
decay[decay>0]=0
ax_signals[1].plot(t_spec_flux[1:],np.array(attacks>ATTACK_THRESH,dtype='float'))
#ax_signals[1].plot(t_spec_flux[1:],attacks)
ax_signals[2].plot(t_spec_flux[1:],decay)
for ax in ax_signals:
    ax.set_xlim(0,2)

plt.tight_layout()
plt.show()

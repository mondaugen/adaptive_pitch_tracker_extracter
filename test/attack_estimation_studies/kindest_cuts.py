# Take a sound file
# Do some processing to get the attacks, be it spectral flux or high-frequency weighting
# Look for the local peaks. In order to filter out spurious peaks, we do a sort
# of noise-gating of the spectral flux signal by looking at the RMS amplitude

import numpy as np
from scipy import signal
import common
import spectral_difference
import subprocess
import high_freq_weighting_filter
import matplotlib.pyplot as plt

SAMPLE_RATE=common.get_env('SAMPLE_RATE',conv=float,default=16e3)
INPUT=common.get_env('INPUT',check_if_none=True)
OUTPUT=common.get_env('OUTPUT',default='/tmp/cut-%d.f64')


# window type for the spectral flux
WINDOW_TYPE_SF=common.get_env('WINDOW_TYPE_SF',default='hann')
# hop size of spectral flux
H_SF=common.get_env('H_SF',conv=int,default=256)
# window size of DFT implementation
W_SF=common.get_env('W_SF',conv=int,default=1024)
# smoothing factor (0-1) for SF
SMOOTH_SF=common.get_env('SMOOTH_SF',conv=float,default=1)
# max filter discount rate, the time in seconds until it reaches 1% of the
# maximum
LMAX_FILT_RATE=common.get_env('LMAX_FILT_RATE',conv=float,default=1)
# calculate number of hops until it hits 1%
# time per hop
s_per_H=H_SF/SAMPLE_RATE
# number of hops in the LMAX_FILT_RATE
lmfr_n_H=LMAX_FILT_RATE/s_per_H
if lmfr_n_H <= 0:
    lmfr_n_H=1
# what number to this power is .01 ?
lmfr=np.power(.01,1/lmfr_n_H)
print("hops per LMAX_FILT_RATE %f" % (lmfr_n_H,))
print("discount rate: %f" % (lmfr,))

# the threshold in dB for the noise gate
NG_TH=common.get_env('NG_TH',conv=float,default=-60)

x=np.fromfile(INPUT,dtype='float64')
x_t=np.arange(len(x))/SAMPLE_RATE

# the boundaries for plotting
X_LIM=common.get_env('X_LIM',conv=eval,default=(0,x_t[-1]))

sd=spectral_difference.spectral_diff(x,H_SF,W_SF,WINDOW_TYPE_SF)
sd=spectral_difference.iir_avg(sd,SMOOTH_SF)
sd_maxs,sd_max_thresh=spectral_difference.discount_local_max(sd,lmfr)
# set one_sided_max to 'left' so that the first point to become non-zero is
# included as a local minimum
sd_mins=spectral_difference.local_max(-sd,one_sided_max='left')
sd_mins_filtered=spectral_difference.closest_index_after(sd_maxs,sd_mins,reverse=True)
sd_t=((np.arange(len(sd))+1)*H_SF+W_SF*0.5)/SAMPLE_RATE

x_rms=spectral_difference.local_rms(x,H_SF,W_SF)
sd_gate=(x_rms>np.power(10,NG_TH/20)).astype('float64')
sd_gate_t=np.arange(len(sd_gate))*H_SF/SAMPLE_RATE
sd_maxs=spectral_difference.index_mask(sd_maxs,sd_gate)

sd*=sd_gate[1:]

fig,axs=plt.subplots(4,1)
axs[0].plot(x_t,x)
axs[0].set_title('Original')
axs[1].plot(sd_t,sd)
axs[1].set_title('Masked spectral flux')
axs[1].plot(sd_t[sd_maxs],sd[sd_maxs],'.r')
axs[1].plot(sd_t[sd_mins_filtered],sd[sd_mins_filtered],'.g')
axs[1].plot(sd_t,sd_max_thresh)
#axs[2].set_title('Local max threshold')
axs[3].plot(sd_gate_t,sd_gate)
axs[3].set_title('RMS gate')

for ax in axs:
    ax.set_xlim(*X_LIM)

# now cut out the sections
# sd_maxs are the indices of the peaks in the spectral flux function
# sd_mins_filtered are the indices of the closest local minimum just before the peaks
# if the first sd_maxs is greater than the first sd_mins_filtered, we put 0 as
# the first minimum. This is not great but usually we can catch the first
# minimum using one_sided_max='left'
if sd_maxs[0] > sd_mins_filtered[0]:
    sd_mins_filtered=np.concatenate(([0],sd_mins_filtered))
# if the last sd_maxs is greater than the last sd_mins_filtered, we add the end
# of the file as the last sd_mins_filtered index
if sd_maxs[-1] > sd_mins_filtered[-1]:
    sd_mins_filtered=np.concatenate((sd_mins_filtered,[len(x)-1]))
# now we combine all the indices
attack_points=np.sort(np.concatenate((sd_maxs,sd_mins_filtered)))
# convert to input time
attack_points=attack_points*H_SF
# If there is an even number of indices, there is an error
if (len(attack_points) % 2) == 0:
    raise Exception("There can't be an even number of attack_points")

# now cut out the samples by taking every 3 points
for k,n in enumerate(range(0,len(attack_points)-2,2)):
    s,m,e=attack_points[n:n+3]
    y=x[s:e+1]
    y.tofile(OUTPUT % (k,))

plt.tight_layout()
plt.show()

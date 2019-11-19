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
import window_tools

SAMPLE_RATE=common.get_env('SAMPLE_RATE',conv=float,default=16e3)
INPUT=common.get_env('INPUT',check_if_none=True)
OUTPUT=common.get_env('OUTPUT',default='/tmp/cut-%d.f64')
OUTPUT_SUMMARY=common.get_env('OUTPUT_SUMMARY',default='/tmp/cut_summary.txt')

# window type for the spectral flux
WINDOW_TYPE_SF=common.get_env('WINDOW_TYPE_SF',default='hann')
# hop size of spectral flux
H_SF=common.get_env('H_SF',conv=int,default=256)
# window size of DFT implementation
W_SF=common.get_env('W_SF',conv=int,default=1024)
# smoothing factor (0-1) for SF
SMOOTH_SF=common.get_env('SMOOTH_SF',conv=float,default=1)
SF_THRESH=common.get_env('SF_THRESH',conv=float,default=0.01)
MAX_SEG_LEN=common.get_env('MAX_SEG_LEN',conv=float,default=float('inf'))
# If not None, Instead of using the most recent local minimum in the spectral
# flux before an attack point as the start time, simply subtract
# FORCE_START_TIME from the attack time to get the start time
FORCE_START_TIME=common.get_env('FORCE_START_TIME',conv=int,default=None)
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
sd_maxs,sd_max_thresh=spectral_difference.discount_local_max(sd,lmfr,min_thresh=SF_THRESH)
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
if np.min(sd_maxs) < np.min(sd_mins_filtered):
    print("Adding extra sd_mins_filtered to beginning.")
    sd_mins_filtered=np.concatenate(([0],sd_mins_filtered))
# if the last sd_maxs is greater than the last sd_mins_filtered, we add the end
# of the file as the last sd_mins_filtered index
if np.max(sd_maxs) > np.max(sd_mins_filtered):
    print("Adding extra sd_mins_filtered to end.")
    sd_mins_filtered=np.concatenate((sd_mins_filtered,[len(x)-1]))
# now we combine all the indices
attack_points=np.sort(np.concatenate((sd_maxs,sd_mins_filtered)))
# check to see if there are any adjacent minima or maxima
sd_mins_mask=np.zeros(np.max(attack_points)+1)
sd_mins_mask[sd_mins_filtered]=1
sd_maxs_mask=np.zeros(np.max(attack_points)+1)
sd_maxs_mask[sd_maxs]=1
mask_accum=np.cumsum(sd_maxs_mask)-np.cumsum(sd_mins_mask)
print(np.where((mask_accum>0)|(mask_accum<-1))[0])
# convert to input time
attack_points_sf=attack_points*H_SF
# If there is an even number of indices, there is an error
if (len(attack_points_sf) % 2) == 0:
    print("WARNING: There are an even number of attack_points.")

# a thing to help taper the samples
twa=window_tools.taper_window_applier()
# also make a summary file
with open(OUTPUT_SUMMARY,'w') as f:
    # now cut out the samples by taking every 3 points
    for k,n in enumerate(range(0,len(attack_points_sf)-2,2)):
        s,m,e=attack_points_sf[n:n+3]
        if FORCE_START_TIME is not None:
            s=m-FORCE_START_TIME
        if ((e-s) > MAX_SEG_LEN) or ((e-s) < W_SF):
            continue
        # this is a hack to get an inverse ramp
        start_taper=1-twa.get_taper(m-s,where='end')
        end_taper_len=min(e-m,100)
        end_taper=twa.get_taper(end_taper_len,where='end')
        y=x[s:e]
        y[:m-s]*=start_taper
        y[-end_taper_len:]*=end_taper
        y.tofile(OUTPUT % (k,))
        f.write((OUTPUT % (k,)) + (" length=%d attack=%d\n" % (e-s,m-s)))

# plot the cut ranges on the graph
sd_mins_filtered_s=np.sort(((sd_mins_filtered+1)*H_SF+W_SF*0.5)/SAMPLE_RATE)
heights=[0.1,-0.1]
h_n=0
# compensate for spectral flux windowing when converting to seconds
attack_points_s=((attack_points+1)*H_SF+W_SF*0.5)/SAMPLE_RATE
for n in range(0,len(attack_points_s)-2,2):
    l,c,r=attack_points_s[n:n+3]
    axs[0].plot([l,l,r,r],[0,heights[h_n],heights[h_n],0],'k')
    axs[0].plot([c],[heights[h_n]],'.r')
    h_n = 1 - h_n

plot_local_x_max=False

if plot_local_x_max:

    # create an envelope using the discounted local max of the signal
    x_max_rate=spectral_difference.discount_local_max_rate_calc(SAMPLE_RATE*0.1)
    x_maxs,x_max_thresh=spectral_difference.discount_local_max(x,x_max_rate)
    x_maxs_sd_attacks=(sd_maxs+2)*H_SF+W_SF//2
    x_maxs_interp=np.interp(x_t,x_t[x_maxs],x[x_maxs])
    axs[2].plot(x_t,x_maxs_interp)
    axs[2].set_title('Envelope using local max')
    # where are the local maxima as given by the spectral flux on this function?
    axs[2].plot(x_t[x_maxs_sd_attacks],x_maxs_interp[x_maxs_sd_attacks],'.r')

else:
    axs[2].plot(np.arange(len(mask_accum))/SAMPLE_RATE,mask_accum)

plt.tight_layout()
plt.show()

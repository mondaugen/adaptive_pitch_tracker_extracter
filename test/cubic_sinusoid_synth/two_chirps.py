# Tests the cubic sinusoid synth by trying to synthesize chirps
# This is a qualitative test because an exact resynthesis of the original chirps
# is difficult using the simple analysis done here: taking the local maxima of
# the power spectrum and their corresponding phases in the first and last frame
# to get parameters for the chirp.

from cubic_sinusoid_synth import cubic_sinusoid_synth
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from common import normalize

## Synthesis
# Sample rate
FS=16e3
# Length of signal, total length will be N+W
N=1e4
# Window size for analysis
W=1024
# Hop size for plotting
H_plot=256
# Sample times
n=np.arange(N+W)
# Chirp start frequencies
f0=[100,200]
# Chirp end frequencies
f1=[200,100]
# Chirps, separate (not summed)
chirps=[signal.chirp(n,f0_/FS,N,f1_/FS) for f0_,f1_ in zip(f0,f1)]
# Number of chirps
n_chirps=len(chirps)
# Synthesize chirps and save to disk
chirp_sum=np.sum(np.vstack(chirps),axis=0)
chirp_sum_export=normalize(chirp_sum)
chirp_sum_export.tofile('/tmp/chirp_sum.f64')

## Analysis
# Analyses of first frame
frame0=[np.fft.rfft(chirp[:W]) for chirp in chirps]
# Analyses of last frame
frame1=[np.fft.rfft(chirp[-W:]) for chirp in chirps]
# Frequency of maximum bin, first frame (angular velocity, radians per sample)
frame0_max_f=[np.argmax(np.abs(frame))/W*2*np.pi for frame in frame0]
# Frequency of maximum bin, last frame (angular velocity, radians per sample)
frame1_max_f=[np.argmax(np.abs(frame))/W*2*np.pi for frame in frame1]
# Phase of maximum bin, first frame (radians)
frame0_max_f_phase=[np.angle(frame[np.argmax(np.abs(frame))]) for frame in frame0]
frame1_max_f_phase=[np.angle(frame[np.argmax(np.abs(frame))]) for frame in frame1]

## Resynthesis
# A cubic sinusoid synth where the hop size is N
part_synth=cubic_sinusoid_synth(n_chirps,N)
part_synth.set_theta_k0(np.array(frame0_max_f_phase))
part_synth.set_omega_k0(np.array(frame0_max_f))
part_synth.set_a_k0(np.ones(n_chirps))
Th,A=part_synth.process(np.array(frame1_max_f_phase),np.array(frame1_max_f),np.ones(n_chirps))
resynth_chirp_sum=np.sum(np.cos(Th),axis=0)
resynth_chirp_sum_export=normalize(resynth_chirp_sum)
resynth_chirp_sum.tofile('/tmp/resynth_chirp_sum.f64')


fig,axs=plt.subplots(1,2,squeeze=False)
axs[0,0].specgram(chirp_sum,NFFT=W,Fs=FS,noverlap=W-H_plot)
axs[0,0].set_title('chirp sum')
axs[0,1].specgram(resynth_chirp_sum,NFFT=W,Fs=FS,noverlap=W-H_plot)
axs[0,1].set_title('resynth chirp sum')
plt.show()



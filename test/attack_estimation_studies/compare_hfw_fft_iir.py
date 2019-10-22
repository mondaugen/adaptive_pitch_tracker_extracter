# Compare the performance of the DFT and IIR implementations of the high
# frequency weighting filter

import numpy as np
from scipy import signal
import common
import spectral_difference
import subprocess
import high_freq_weighting_filter
import matplotlib.pyplot as plt

INPUT=common.get_env('INPUT',check_if_none=True)
WINDOW_TYPE=common.get_env('WINDOW_TYPE',default='hann')
SAMPLE_RATE=common.get_env('SAMPLE_RATE',conv=float,default=16e3)
# hop size of DFT implementation
H=common.get_env('H',conv=int,default=2)
# window size of DFT implementation
W=common.get_env('W',conv=int,default=8)
# number of poles in IIR implementation
P=common.get_env('P',conv=int,default=2)
X_LIM=common.get_env('X_LIM',conv=eval,default=(0,2))
# Max search window sizes
H_LMAX=common.get_env('H_LMAX',conv=int,default=5)
W_LMAX=common.get_env('W_LMAX',conv=int,default=10)
# How many times greater the max has to be to get counted
A_LMAX=common.get_env('A_LMAX',conv=float,default=1.5)
# Coefficient of smoothing filter
A_SMOOTH=common.get_env('A_SMOOTH',conv=float,default=1e-2)
# Max search window sizes for averaging filter
H_ASLMAX=common.get_env('H_ASLMAX',conv=int,default=5)
W_ASLMAX=common.get_env('W_ASLMAX',conv=int,default=10)
A_ASLMAX=common.get_env('A_ASLMAX',conv=float,default=A_LMAX)

# Design high frequency weighting filter
b,a=high_freq_weighting_filter.fit_allpole_triangular_highpass(P)
r=high_freq_weighting_filter.a_to_r(a)
r.astype('float32').tofile('/tmp/r.f32')

# run IIR hfwf
iir_env={
    'B0': '%.18f' % (b[0],),
    'IN_PATH': INPUT
}

subprocess.run('src/test/bin/iir_lattice_filter_proc',env=iir_env)
x=np.fromfile(INPUT,dtype='float32')
t_x=np.arange(len(x))/SAMPLE_RATE
y_iir=np.fromfile('/tmp/out.f32',dtype='float32')[:len(x)]
# smoothed using IIR averaging filter
y_iir_smooth_avg=spectral_difference.iir_avg(np.abs(y_iir),A_SMOOTH)
y_iir=spectral_difference.local_rms(y_iir,H,W)
t_y_iir=np.arange(len(y_iir))*H/SAMPLE_RATE
y_iir_maxs=spectral_difference.filtered_local_max(y_iir,H_LMAX,W_LMAX,A_LMAX)
y_iir_maxs_inv=spectral_difference.filtered_local_max(-y_iir,H_LMAX,W_LMAX,1)
# When looking for max in smooth average, we use same total window size as when looking in RMS
y_iir_smooth_avg_maxs=spectral_difference.filtered_local_max(
    y_iir_smooth_avg,H_ASLMAX,W_ASLMAX,A_ASLMAX)
y_dft=spectral_difference.high_frequency_weight(x,H,W,WINDOW_TYPE)
t_y_dft=np.arange(len(y_dft))*H/SAMPLE_RATE
y_dft_maxs=spectral_difference.filtered_local_max(y_dft,H_LMAX,W_LMAX,A_LMAX)
y_dft_maxs_inv=spectral_difference.filtered_local_max(-y_dft,H_LMAX,W_LMAX,A_LMAX)

fig,axs=plt.subplots(4,1)
axs[0].plot(t_x,x)
axs[0].set_title('original')
axs[1].plot(t_y_iir,y_iir)
axs[1].plot(t_y_iir[y_iir_maxs],y_iir[y_iir_maxs],'.')
axs[1].plot(t_y_iir[y_iir_maxs_inv],y_iir[y_iir_maxs_inv],'.g')
axs[1].set_title('IIR HFW')
axs[2].plot(t_y_dft,y_dft)
axs[2].plot(t_y_dft[y_dft_maxs],y_dft[y_dft_maxs],'.')
axs[2].plot(t_y_dft[y_dft_maxs_inv],y_dft[y_dft_maxs_inv],'.g')
axs[2].set_title('DFT HFW')
axs[3].plot(t_x,y_iir_smooth_avg)
axs[3].plot(t_x[y_iir_smooth_avg_maxs],y_iir_smooth_avg[y_iir_smooth_avg_maxs],'.')
axs[3].set_title('IIR smooth average')

for ax in axs:
    ax.set_xlim(*X_LIM)

plt.tight_layout()
plt.show()

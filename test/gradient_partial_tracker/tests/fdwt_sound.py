from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type
import numpy as np
import matplotlib.pyplot as plt
from dftdk import gradient_ascent_step_harm_lock
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal, interpolate
from fdw_tracker import fdw_tracker, v_dev_stop_criterion, multi_analyse_combine
from cubic_sinusoid_synth import quadratic_phase_poly_interp, linear_amplitude_interp
from common import normalize, get_env
import partial_processing as pp
from buf import inf_buf

IN_FILE=get_env('IN_FILE',default='sounds/tanguillo-1ch-48k.f64')
OUT_FILE=get_env('OUT_FILE',default='/tmp/a.f64')
SR=get_env('SR',default=44100,conv=float)
T_START=get_env('T_START',default=0.6,conv=float)
F_START=get_env('F_START',default=150.,conv=float)
R=get_env('R',default=0.02,conv=float) # output rate

x=np.fromfile(IN_FILE)
sr=SR
t_start=T_START
n_start=int(np.round(t_start*sr))

# starting frequency
vstart=np.array([F_START/sr])
n_harms=30
vstarts=np.multiply.outer(vstart,(1+np.arange(n_harms))).flatten()
v_groups=np.multiply.outer(np.arange(len(vstart)),np.ones(n_harms)).flatten().astype('int')
print(v_groups)
N=2048
N_h=N//8
fdwt=fdw_tracker(N=N)

def sc():
    return v_dev_stop_criterion(vstarts[0],max_dev_semitones=0.5)

wrapped_x=inf_buf(x)
h_fw=np.arange(n_start,len(x),N_h)
h_bw=np.arange(n_start,0,-N_h)

k,h,X,_ = multi_analyse_combine(fdwt,x=wrapped_x,vstarts=vstarts,h=[h_bw,h_fw],v_groups=v_groups,warm_up_hops=0,mu=0.2,n_steps=3,grad_weight='equal',stop_criterion=sc)

h_plot = h

plot_tracks=True
if plot_tracks:
    for k_row in k:
        plt.plot(h_plot/sr,k_row[:len(h)]/N*sr)

plt.specgram(x,NFFT=N,Fs=sr,
window=fdwt.win_type.window(N,centered_at_time_zero=False),
noverlap=N-np.abs(N_h),sides='onesided')

X_synth=X
X_abs=np.abs(X_synth)
X_phs=np.angle(X_synth)
omega=k/N*2*np.pi
if N_h < 0:
    omega=omega[:,::-1]

# output rate for time stretching (<1) or compression (>1)
output_rate=R
orig_rate_len=min(100,omega.shape[1])
input_indices=np.arange(0,omega.shape[1])
output_indices=np.concatenate((
    np.arange(0,orig_rate_len),
    np.arange(orig_rate_len,omega.shape[1]-1+output_rate,output_rate))
)
X_abs=pp.squish_anomalies(X_abs)

phase_track_interp=interpolate.interp1d(
input_indices,
omega,
axis=1)
omega_interp=phase_track_interp(output_indices)

amp_track_interp=interpolate.interp1d(
input_indices,
X_abs,
axis=1)
X_abs_interp=amp_track_interp(output_indices)

phase_tracks=quadratic_phase_poly_interp(np.abs(N_h),X_phs[:,0],omega_interp.T)
amp_tracks=linear_amplitude_interp(np.abs(N_h),X_abs_interp.T)

x_synth=np.real((np.exp(complex('j')*phase_tracks)*amp_tracks).sum(axis=0))
normalize(x_synth).tofile(OUT_FILE)

plt.show()

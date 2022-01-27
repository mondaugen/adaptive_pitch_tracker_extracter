from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type
import numpy as np
import matplotlib.pyplot as plt
from dftdk import gradient_ascent_step_harm_lock
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal
from fdw_tracker import fdw_tracker

x=np.fromfile('sounds/bovel-1ch-48k.f64')
sr=48e3
t_start=1.5
n_start=int(np.round(t_start*sr))

# starting frequency
vstart=np.array([400./sr])
n_harms=10
vstarts=np.multiply.outer(vstart,(1+np.arange(n_harms))).flatten()
v_groups=np.multiply.outer(np.arange(len(vstart)),np.ones(n_harms)).flatten().astype('int')
print(v_groups)
N=2048
N_h=N//8
fdwt=fdw_tracker(N=N,N_h=N_h)

k,h = fdwt.analyse(x[n_start:],vstarts,v_groups=v_groups,warm_up_hops=0,mu=0.5,n_steps=3)

plot_tracks=False
if plot_tracks:
    for k_row in k:
        plt.plot((h+n_start)/sr,k_row[:len(h)]/N*sr)

plt.specgram(x,NFFT=N,Fs=sr,
window=fdwt.win_type.window(N,centered_at_time_zero=False),
noverlap=N-N_h,sides='onesided')

plt.show()

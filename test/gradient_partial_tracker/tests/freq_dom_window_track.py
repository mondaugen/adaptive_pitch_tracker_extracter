from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type
import numpy as np
import matplotlib.pyplot as plt
from dftdk import gradient_ascent_step_harm_lock
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal
from fdw_tracker import fdw_tracker

# synthesize signal to analyse
N_x=int(1e5)
n_x=np.arange(N_x)
v0=0.03
phi=0
v_max=0.5
partial_amp = lambda p: 1/(p+1)
x=comb_no_mod(n_x,v0,phi=phi,v_max=v_max,partial_amp=partial_amp)
x+=comb_no_mod(n_x,v0*np.power(2,4/12.),phi=phi,v_max=v_max,partial_amp=partial_amp)
x+=np.random.standard_normal((N_x,))*0.5

# analysis parameters
starting_error=1 # semitones
# starting frequency
vstart=np.array([v0*np.power(2,starting_error/12.)])
N=2048
N_h=N//4
fdwt=fdw_tracker(N=N,N_h=N_h)
k,h = fdwt.analyse(x,vstart)

for k_row in k:
    plt.plot(h,k_row[:len(h)]/N)

plt.specgram(x,NFFT=N,Fs=1,
window=fdwt.win_type.window(N,centered_at_time_zero=False),
noverlap=N-N_h,sides='onesided')

plt.show()

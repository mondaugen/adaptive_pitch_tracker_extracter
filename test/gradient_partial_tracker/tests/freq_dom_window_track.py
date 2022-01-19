from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type
import numpy as np
import matplotlib.pyplot as plt
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal

def sq_mod(x):
    return np.real(x*np.conj(x))

def sq_mod_dv(X,dX):
    return np.real(X*np.conj(dX)+np.conj(X)*dX)

def log_sq_mod_dv(X,dX):
    return np.real(np.conj(dX)/np.conj(X)+dX/X)

N=2048
n=np.arange(N)
v0=0.03
phi=0
v_max=0.5
partial_amp = lambda p: np.ones_like(p)
x=comb_no_mod(n,v0,phi=phi,v_max=v_max,partial_amp=partial_amp)

starting_error=0.2
# starting frequency
vstart=np.array([v0*np.power(2,starting_error/12.)])
print("starting frequency:",vstart*Fs)
# hop size
N_w=L
N_h=N_w//4
h=np.arange(0,N_x-N,N_h)
# gradient step coefficient
mu=1e-4 # if 0, investigating accuracy of gradient calculation
# number of gradient steps per frame
n_steps=10

fdw=freq_dom_window(N,
    sum_of_cos_dft_win_type(win_coeffs,win_lobe_radius),
    dft_os)

X_dft_R=fdw.R(dft_v)
X_dft_fdw=dft(x,X_dft_R)
# derivative
X_dft_dv_fdw=dft_dv(x,X_dft_R)
X_dft_dv_fdw_plot=log_sq_mod_dv(X_dft_fdw,X_dft_dv_fdw)

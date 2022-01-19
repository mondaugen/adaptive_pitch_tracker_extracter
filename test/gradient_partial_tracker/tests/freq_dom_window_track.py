from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type
import numpy as np
import matplotlib.pyplot as plt
from dftdk import gradient_ascent_step_harm_lock
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal

def sq_mod(x):
    return np.real(x*np.conj(x))

def sq_mod_dv(X,dX):
    return np.real(X*np.conj(dX)+np.conj(X)*dX)

def log_sq_mod_dv(X,dX):
    return np.real(np.conj(dX)/np.conj(X)+dX/X)

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
# frame size
N=2048
# hop size
N_w=N
N_h=N_w//4
h=np.arange(0,N_x-N,N_h,dtype='int')
# gradient step coefficient
mu=0.5
# number of gradient steps per frame
n_steps=1

dft_os=32
win_coeffs=[0.35875, 0.48829, 0.14128, 0.01168] # blackmanharris
win_lobe_radius=6
win_type=sum_of_cos_dft_win_type(win_coeffs,win_lobe_radius)
fdw=freq_dom_window(N,
    win_type,
    dft_os)

warm_up_hops=4

def grad(i):
    def _grad(x,k0):
        if i < warm_up_hops:
            return np.zeros_like(k0)
        v0=k0/N
        X_dft_R=fdw.R(v0)
        X_dft_fdw=dft(x,X_dft_R)
        # derivative
        X_dft_dv_fdw=dft_dv(x,X_dft_R)
        return log_sq_mod_dv(X_dft_fdw,X_dft_dv_fdw)
    return _grad

n_harms=10
k0=vstart*N*(1+np.arange(n_harms))
buf=np.zeros(N,dtype=x.dtype)
k=np.zeros((len(k0),len(h)+1))
k[:,0] = k0
for i, n_h in enumerate(h):
    buf[:]=x[n_h:n_h+N]
    for m in range(n_steps):
        k0=gradient_ascent_step_harm_lock(buf,k0,mu,grad=grad(i))
    k[:,i+1]=k0


for k_row in k:
    plt.plot(h,k_row[:len(h)]/N)

plt.specgram(x,NFFT=N,Fs=1,
window=win_type.window(N,centered_at_time_zero=False),
noverlap=N-N_h,sides='onesided')

plt.show()

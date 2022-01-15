from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft
import numpy as np
import matplotlib.pyplot as plt
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal

def sq_mod(x):
    return np.real(x*np.conj(x))

def sq_mod_dv(X,dX):
    return np.real(X*np.conj(dX)+np.conj(X)*dX)

N=2048
n=np.arange(N)
v0=0.03
phi=0
v_max=0.5
partial_amp = lambda p: np.ones_like(p)
x=comb_no_mod(n,v0,phi=phi,v_max=v_max,partial_amp=partial_amp)

# amount of DFT oversampling in frequency
dft_os=32
x_dft=np.zeros(N*dft_os,dtype=x.dtype)
x_dft[:N//2]=x[-(N//2):]
x_dft[-(N//2):]=x[:N//2]
win_type='blackmanharris'
win_coeffs=[0.35875, 0.48829, 0.14128, 0.01168] # blackmanharris
win_lobe_radius=4

L=N
W=N-1
w=sum_of_cos(normalize_sum_of_cos_A(win_coeffs,N,N-1,N*dft_os),N,N-1,N*dft_os)
x_dft[:N//2]*=w[:N//2]
x_dft[-(N//2):]*=w[-(N//2):]
X_dft=np.fft.fft(x_dft)
dft_v=np.arange(0,N,1./dft_os)/N
X_dft_plot=sq_mod(X_dft)

plt.plot(dft_v,X_dft_plot)
plt.title("DFT of signal")

plt.figure()
plt.plot(dft_v[1:]-0.5/dft_os,np.diff(X_dft_plot))
plt.title("Finite differences of DFT of signal")

# now find it with freq_dom_window

# callable win_type parameter
class sum_of_cos_dft_win_type:
    def __init__(self,A,lobe_radius):
        """ A are the sum-of-cos coefficients """
        self.A=A
        self.lobe_radius=lobe_radius
    def __call__(self,N_win,oversample):
        evalradius=self.lobe_radius*oversample
        evalbins=np.arange(-evalradius,evalradius+1)
        bins=evalbins/oversample
        L=N_win
        W=L-1
        N=N_win*oversample
        A_norm=normalize_sum_of_cos_A(self.A,L,W,N)
        vals=sum_of_cos_dft(evalbins,A_norm,L,W,N)
        return bins,vals

fdw=freq_dom_window(N,
    sum_of_cos_dft_win_type(win_coeffs,win_lobe_radius),
    dft_os)

X_dft_R=fdw.R(dft_v)
X_dft_fdw=dft(x,X_dft_R)
X_dft_fdw_plot=sq_mod(X_dft_fdw)
# derivative
X_dft_dv_fdw=dft_dv(x,X_dft_R)
X_dft_dv_fdw_plot=sq_mod_dv(X_dft_fdw,X_dft_dv_fdw)

plt.figure()
plt.plot(dft_v,X_dft_fdw_plot)
plt.title("DFT of signal via FDW")

plt.figure()
plt.plot(dft_v,X_dft_dv_fdw_plot)
plt.title("Derivative of DFT of signal")

plt.show()

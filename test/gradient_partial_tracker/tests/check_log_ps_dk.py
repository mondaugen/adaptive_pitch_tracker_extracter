# Check that we can correctly compute the derivative w.r.t. k of the log power
# spectrum of a sinusoid, regardless of its phase.

import numpy as np
import dftdk
from some_sig import mod_sum_of_cos, sum_of_cos
from some_ft import mod_sum_of_cos_dft, mod_sum_of_cos_dft_dk, mod_sum_of_cos_dft_k, normalize_sum_of_cos_A
from dftdk import dft_bin
import matplotlib.pyplot as plt

j=complex('j')

def sq_mag(X):
    return np.real(np.conj(X)*X)

def log_ps(X):
    return np.log(sq_mag(X))

def shift_x(x,N):
    return np.concatenate((x[N//2:],x[:N//2]))

# length of signal
N=1024
# phase of sinusoid to be analysed
ph=0.1*2*np.pi
# frequency of sinusoid to be analysed
v=0.03
# sample times
n=np.arange(N)
# sinusoid to be analysed
x=np.exp(j*(2*np.pi*v*n + ph))

plt.plot(n,np.real(x),label='analysed signal')
plt.title("Time-domain signal")

# see some_ft for descriptions of these variables
L=512
W=L-1
raw_A=[0.42,0.5,0.08] # blackman-harris
A=normalize_sum_of_cos_A(raw_A,L,W,N)
# bins
k=n
# its theoretical log power spectrum
# NOTE: the natural log is used (not log10 like for dB) because this function
# has a slightly simpler derivative
X_th=mod_sum_of_cos_dft(v,k,A,L,W,N)
X_lps_th=log_ps(X_th)

hgtd=dftdk.harm_grad_td(raw_A,L,W,N,
    harm_sig=dftdk.harm_grad_td.default_harm_sig(B=np.array([1.]))
)
X_lps=log_ps(np.array([hgtd.X(x,k_) for k_ in k]))
print("X_lps:",X_lps)

# computed with fourier transform
w=sum_of_cos(A,L,W,N)
X_lps_ft=log_ps(np.fft.fft(shift_x(x,N)*w))

plt.figure()
plt.plot(k,X_lps_th,label="X_lps_th")
plt.plot(k,X_lps,label="X_lps")
plt.plot(k,X_lps_ft,label="X_lps_ft")
plt.title("Log power spectrum")
plt.legend()

X_lps_dk_fd=np.diff(X_lps_th)

plt.figure()
plt.plot(k[1:]-0.5,X_lps_dk_fd,label="X_lps_dk_fd")
plt.title("d/dk log power spectrum")


dX_log_ps_dk=np.array([hgtd.d_log_ps_dk_fd(x,k_) for k_ in k])
plt.plot(k,dX_log_ps_dk,label="X_lps_dk")
plt.legend()

plt.figure()
plt.plot(k,w)
plt.title("window")

w_ft=np.fft.fft(w)
plt.figure()
plt.plot(k,np.real(w_ft))
plt.plot(k,np.imag(w_ft))
plt.title("window's fourier transform")

# I dunno why, but the log power spectrum derivative is wrong...
plt.show()

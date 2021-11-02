# find the derivative of the power spectrum of a signal build from a harmonic
# series of sinusoids using the fft and using the theoretical functions
# check this for different units of frequency

from some_sig import mod_sum_of_cos
from some_ft import mod_sum_of_cos_dft, mod_sum_of_cos_dft_dk, mod_sum_of_cos_dft_k
from dftdk import dft_bin
import numpy as np
import matplotlib.pyplot as plt

# synthesize time-domain signal and find DFT, DFT/df using fft
# compute theoretical DFT with derivative computed using finite differences
# compute theoretical DFT/df
# these should all be very close

def check(
    v=0.25,
    N=128,
    L=64,
    W=63,
    A=[0.5,0.5],
    Fs=16e3
):
    k=np.arange(N)
    k0=v*N
    w=mod_sum_of_cos(v,A,L,W,N)
    w_ft=np.fft.fft(w)
    w_ft_re=np.real(w_ft)
    w_ft_dk_fd=np.diff(w_ft_re)
    w_ft_th=mod_sum_of_cos_dft_k(k0,k,A,L,W,N)
    w_ft_ft=dft_bin(w,k,p=1)
    w_ft_ft_re=np.real(w_ft_ft)
    w_ft_dk_th=mod_sum_of_cos_dft_dk(k0,k,A,L,W,N)


    plt.figure()
    plt.title('DFT')
    plt.plot(k,w_ft,label='FFT')
    plt.plot(k,w_ft_th,label='DFT th')
    plt.legend()
    plt.figure()
    plt.title('dv DFT')
    plt.plot(k[:-1]+0.5,w_ft_dk_fd,label='dk FD')
    plt.plot(k,w_ft_ft_re,label='dk FFT')
    plt.plot(k,w_ft_dk_th,label='dk th')
    plt.legend()

    # different units
    w_ft_df_fd=np.diff(w_ft_re)/(Fs/N)
    w_ft_df_ft=w_ft_ft_re*(N/Fs)
    w_ft_df_th=w_ft_dk_th*(N/Fs)
    plt.figure()
    plt.title('dv DFT w.r.t. frequency in hz')
    plt.plot(k[:-1]+0.5,w_ft_df_fd,label='df fd') 
    plt.plot(k,w_ft_df_ft,label='df FFT') 
    plt.plot(k,w_ft_df_th,label='df th') 
    plt.legend()

check()
plt.show()

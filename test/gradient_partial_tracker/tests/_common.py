import numpy as np
from some_sig import mod_sum_of_cos
from some_ft import mod_sum_of_cos_dft, mod_sum_of_cos_dft_dk, mod_sum_of_cos_dft_k, normalize_sum_of_cos_A
from dftdk import dft_bin
import matplotlib.pyplot as plt

def check_dv_dft(
    v=0.25,
    N=128,
    L=64,
    W=63,
    A=[0.5,0.5],
    Fs=16e3,
    x_td=mod_sum_of_cos,
    x_dft=mod_sum_of_cos_dft_k,
    x_ddft=mod_sum_of_cos_dft_dk,
    normalize=True
):
    if normalize:
        A=normalize_sum_of_cos_A(A,L,W,N)
    k=np.arange(N)
    k0=v*N
    w=x_td(v,A,L,W,N)
    w_ft=np.fft.fft(w)
    w_ft_re=np.real(w_ft)
    w_ft_dk_fd=np.diff(w_ft_re)
    w_ft_th=x_dft(k0,k,A,L,W,N)
    w_ft_ft=dft_bin(w,k,p=1)
    w_ft_ft_re=np.real(w_ft_ft)
    w_ft_dk_th=x_ddft(k0,k,A,L,W,N)


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


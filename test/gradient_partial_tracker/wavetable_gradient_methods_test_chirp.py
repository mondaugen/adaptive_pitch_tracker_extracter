import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from dftdk import dft_dk, dft_bin, ps_dk, dft_bin_log_ps_dk, harm_grad_td, gradient_ascent_step, window_scale, normalize_sum_of_cos_A
from some_sig import complex_chirp, comb, sum_of_cos

j=complex('j')

# signal to analyse
# sampling rate
Fs=44100 # for something plausible
# length of signal in seconds
T=5
# length in samples
N_x=T*Fs
# starting normalized frequency
v0=220/Fs
# ending normalized frequency
v1=220/Fs

# sum of cosine analysis window coefficients
raw_A=[0.42,0.5,0.08]
# analysis parameters (in samples)
# transform length
N=2048
# window length
W=1023
# length of period of analysis kernel
L=W+1
# normalized A
A=normalize_sum_of_cos_A(raw_A,L,W,N)
# partial multipliers (to give harmonic analysis kernel)
P=10
p=(1+np.arange(P))
# amplitude scalars of each partial
B=1./p

n=np.arange(N_x)
# TODO: Why is this "flat" (i.e., its frequency is too low)
x_sig=comb(n,v0,N_x,v1,partial_amp=lambda p: np.ones_like(p))#/(p*p))
x_sig+=0*comb(n,v0*(2**4/12),N_x,v1*(2**4/12),partial_amp=lambda p: 1./(p*p))
x_noise=0*np.random.standard_normal(N_x)*0.1
SNR=10*np.log10(np.real(x_noise*x_noise.conj()).mean()/np.real(x_sig*x_sig.conj()).mean())
print("SNR:",SNR)
x=x_sig+x_noise
# analysis window
hgt=harm_grad_td(raw_A,L,W,N,
    harm_sig=harm_grad_td.default_harm_sig(B=np.ones_like(B)/P)
)
starting_error=0
# starting frequency
vstart=np.array([v0*np.power(2,starting_error/12.)])
print("starting frequency:",vstart*Fs)
# hop size
N_w=L
N_h=N_w//4
h=np.arange(0,N_x-N,N_h)
# gradient step coefficient
mu=1e-4

def do_gradient_method(
params=(vstart,),
do_step=lambda buf,p: gradient_ascent_step(buf,p[0],mu,
grad=lambda x,v: hgt.dX_dk(x,v*N)), # TODO: use log power spectrum gradient instead
extract_k=lambda p: p[0]):
    buf=np.zeros(N,dtype=x.dtype)
    k=np.zeros((len(vstart),len(h)+1))
    k[:,0] = extract_k(params)
    for i, n_h in enumerate(h):
        buf[:]=x[n_h:n_h+N]
        params=do_step(buf,params)
        if not isinstance(params,tuple):
            params=(params,)
        k[:,i+1]=extract_k(params)
    return k
v_ga=do_gradient_method()

# TODO: why does the spectrogram look like marde?
plt.specgram(x,NFFT=N,Fs=1,window=sum_of_cos(A,L,W,N),noverlap=N-N_h,sides='onesided')
print(v_ga)
for v_ga_ in v_ga:
    for v in np.multiply.outer(p,v_ga_):
        plt.plot(h,v[:-1])

plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from dftdk import dft_dk, dft_bin, ps_dk, dft_bin_log_ps_dk, harm_grad_tracker, gradient_ascent_step, window_scale
from some_sig import complex_chirp, comb

j=complex('j')

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
n=np.arange(N_x)
x_sig=comb(n,v0,N_x,v1,partial_amp=lambda p: 1./(p*p))
x_sig+=0*comb(n,v0*(2**4/12),N_x,v1*(2**4/12),partial_amp=lambda p: 1./(p*p))
x_noise=0*np.random.standard_normal(N_x)*0.1
SNR=10*np.log10(np.real(x_noise*x_noise.conj()).mean()/np.real(x_sig*x_sig.conj()).mean())
print("SNR:",SNR)
x=x_sig+x_noise
# analysis window
hgt=harm_grad_tracker(min_f=v0*Fs,max_f=2*v0*Fs,
window=lambda N: window_scale(signal.get_window('blackman',N)),
q=lambda p: np.ones(len(p)))
N_w=hgt.N_w_max
# starting error in semitones
starting_error=0
# starting frequency
vstart=np.array([v0*np.power(2,starting_error/12.)])
print("starting frequency:",vstart*Fs)
# hop size
N_h=512
h=np.arange(0,N_x-N_w,N_h)
# gradient step coefficient
mu=1

def do_gradient_method(
params=(vstart,),
do_step=lambda buf,p: gradient_ascent_step(buf,p[0],mu,
grad=lambda x,v: hgt.dft_bin_log_ps_dk(x,v,N=N_w)),
extract_k=lambda p: p[0]):
    buf=np.zeros(N_w,dtype=x.dtype)
    k=np.zeros((len(vstart),len(h)+1))
    k[:,0] = extract_k(params)
    for i, n_h in enumerate(np.arange(0,N_x-N_w,N_h)):
        buf[:]=x[n_h:n_h+N_w]
        params=do_step(buf,params)
        if not isinstance(params,tuple):
            params=(params,)
        k[:,i+1]=extract_k(params)
    return k
v_ga=do_gradient_method()

plt.specgram(x,NFFT=N_w,Fs=1,window=signal.get_window('blackman',N_w),noverlap=N_w-N_h,sides='onesided')
print(v_ga)
for v_ga_ in v_ga:
    for v in np.multiply.outer((1+np.arange(1)),v_ga_):
        plt.plot(h,v[:-1])

plt.show()





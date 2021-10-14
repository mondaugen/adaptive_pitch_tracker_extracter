import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from dftdk import dft_dk, dft_bin, ps_dk, gradient_ascent_step, gradient_ascent_step_harm_lock, self_adjusting_ga_step, dft_bin_log_ps_dk
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
x_sig+=comb(n,v0*(2**4/12),N_x,v1*(2**4/12),partial_amp=lambda p: 1./(p*p))
x_noise=np.random.standard_normal(N_x)*0.1
SNR=10*np.log10(np.real(x_noise*x_noise.conj()).mean()/np.real(x_sig*x_sig.conj()).mean())
print("SNR:",SNR)
x=x_sig+x_noise
# analysis window
N_w=4096
N_nz=2048
w=np.zeros(N_w)
w_tmp=signal.get_window('blackman',N_nz)
w_tmp/=w_tmp.sum()
w[N_w//2-N_nz//2:N_w//2+N_nz//2]=w_tmp
#w[:N_w//2]=w_tmp[N_w//2:]
#w[N_w//2:]=w_tmp[:N_w//2]
# starting error in semitones
starting_error=0.5
# starting bin
k0_=v0*np.power(2,starting_error/12.)*N_w
print("starting frequency:",k0_/N_w*Fs)
k0=np.arange(k0_,N_w*0.1,k0_)
# hop size
N_h=512
h=np.arange(0,N_x-N_w,N_h)
# gradient step coefficient
mu=1
# self adjusting ascent parameters
mu_saga=0.3 # ?
step_saga=1

def do_gradient_method(
params=(k0,),
do_step=lambda buf,p: gradient_ascent_step_harm_lock(buf,p[0],mu,grad=dft_bin_log_ps_dk),
extract_k=lambda p: p[0]):
    buf=np.zeros(N_w,dtype=x.dtype)
    k=np.zeros((len(k0),len(h)+1))
    k[:,0] = extract_k(params)
    for i, n_h in enumerate(np.arange(0,N_x-N_w,N_h)):
        buf[:]=x[n_h:n_h+N_w]*w
        params=do_step(np.concatenate((buf[N_w//2:],buf[:N_w//2])),params)
        if not isinstance(params,tuple):
            params=(params,)
        k[:,i+1]=extract_k(params)
    return k
k_ga=do_gradient_method()
#k_saga=do_gradient_method(
#params=(step_saga,k0),
#do_step=lambda buf,p: self_adjusting_ga_step(buf,p[1],p[0],1,mu_saga),
#extract_k=lambda p: p[1])

plt.specgram(x,NFFT=N_w,Fs=1,window=w,noverlap=N_w-N_h,sides='onesided')
for k_ga_ in k_ga:
    plt.plot(h,k_ga_[:-1]/N_w)
plt.figure()
plt.plot(np.arange(N_w),w)

plt.show()




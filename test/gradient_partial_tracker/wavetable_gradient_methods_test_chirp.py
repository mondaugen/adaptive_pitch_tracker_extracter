import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from dftdk import dft_dk, dft_bin, ps_dk, dft_bin_log_ps_dk, harm_grad_td, gradient_ascent_step, window_scale, normalize_sum_of_cos_A
from some_sig import complex_chirp, comb, sum_of_cos
from some_ft import multi_mod_sum_of_cos_dft_k

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
x_sig=comb(n,v0,N_x,v1,v_max=0.0625,partial_amp=lambda p: np.ones_like(p))#/(p*p))
x_sig+=0*comb(n,v0*(2**4/12),N_x,v1*(2**4/12),partial_amp=lambda p: 1./(p*p))
x_noise=0*np.random.standard_normal(N_x)*0.1
SNR=10*np.log10(np.real(x_noise*x_noise.conj()).mean()/np.real(x_sig*x_sig.conj()).mean())
print("SNR:",SNR)
x=x_sig+x_noise
# analysis window
hgt=harm_grad_td(raw_A,L,W,N,
    harm_sig=harm_grad_td.default_harm_sig(B=np.ones_like(B)/P)
)
starting_error=0.4
# starting frequency
vstart=np.array([v0*np.power(2,starting_error/12.)])
print("starting frequency:",vstart*Fs)
# hop size
N_w=L
N_h=N_w//4
h=np.arange(0,N_x-N,N_h)
# gradient step coefficient
mu=1e-4

def do_step_save_frame_spectrum(buf,p):
    frame_num=p[-1]
    np.save('/tmp/frame-%d.npy' % (frame_num,), hgt.dft(buf))
    return gradient_ascent_step(
        buf,p[0],mu,grad=lambda x,v: hgt.d_log_ps_dk(x,v*N)
    )

def do_gradient_method(
params=(vstart,0),
do_step=lambda buf,p: gradient_ascent_step(buf,p[0],mu,
grad=lambda x,v: hgt.d_log_ps_dk(x,v*N)),
extract_k=lambda p: p[0],
x_const=False,
n_steps=10):
    buf=np.zeros(N,dtype=x.dtype)
    if x_const:
        k=np.zeros((len(vstart),n_steps+1))
        k[:,0] = extract_k(params)
        buf[:]=x[h[10]:h[10]+N]
        for i in range(n_steps):
            params=do_step(buf,params)
            if not isinstance(params,tuple):
                params=(params,)
            params+=(0,)
            k[:,i+1]=extract_k(params)
    else:
        k=np.zeros((len(vstart),len(h)+1))
        k[:,0] = extract_k(params)
        for i, n_h in enumerate(h):
            buf[:]=x[n_h:n_h+N]
            for m in range(n_steps):
                params=do_step(buf,params)
                if not isinstance(params,tuple):
                    params=(params,)
            params+=(i,)
            k[:,i+1]=extract_k(params)
    return k

x_const=False
# Why does the gradient tracker "jump" after a while when the signal is stable...?
v_ga=do_gradient_method(
do_step=do_step_save_frame_spectrum,
x_const=x_const,
n_steps=3)

if x_const:
    def _X(k):
        return 20*np.log10(abs(multi_mod_sum_of_cos_dft_k(B,(v0*p)*N,k,A,L,W,N)))
    k=np.arange(0,N,0.1)
    X=_X(k)
    plt.plot(k,X)
    steps=np.arange(v_ga.shape[0])
    for p_ in p:
        k=p_*v_ga[0,:]*N
        plt.plot(k,_X(k))
else:
    plt.specgram(x,NFFT=N,Fs=1,window=sum_of_cos(A,L,W,N,centered_at_time_zero=False),noverlap=N-N_h,sides='onesided')
    print(v_ga*Fs)
    for v_ga_ in v_ga:
        p_x_v_ga=np.multiply.outer(p,v_ga_)
        for v in p_x_v_ga:
            lines=plt.plot(h,v[:-1])
        for i,(h_,v_) in enumerate(zip(h,p_x_v_ga[0,:-1])):
            lines[0].axes.annotate(i,xy=(h_,v_))

plt.figure()
plt.plot(n,x)
plt.title("harmonic chirp signal")

plt.figure()
plt.plot(n[1:]-0.5,np.diff(x))
plt.title("finite-difference of harmonic chirp signal")

plt.show()

# Example of adaptive hill climbing

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import dft_hill_climbing as dhc

def adaptive_dhc(x,v_k,w,mu=1e-6):
    buf=np.zeros_like(w)
    buf[1:]=x[:len(buf)-1]
    v_ks = np.zeros_like(x[len(buf)-1:])
    Xs = np.zeros_like(x[len(buf)-1:],dtype='complex128')
    for n, xn in enumerate(x[len(buf)-1:]):
        buf[:-1] = buf[1:]
        buf[-1] = xn
        v_k = v_k + mu * dhc.dft_bin_pow_dv(buf*w,v_k)
        v_ks[n] = v_k
        Xs[n] = dhc.dft_bin(buf*w,v_k)
    return v_ks, Xs


N=1024
Fs=16e3
f0=1000
f1=1010
f_diff=50
f_k=1008
T=3
NT=T*Fs
w=signal.get_window('hann',N)
w/=np.sum(w)

x0=signal.chirp(np.arange(NT),f0/Fs,NT,f1/Fs)
x1=signal.chirp(np.arange(NT),(f0+f_diff)/Fs,NT,(f1+f_diff)/Fs)
x=x0+x1+np.random.standard_normal(int(NT))

v_ks,Xs=adaptive_dhc(x,f_k/Fs,w,mu=1e-6)

fig,axs=plt.subplots(2,1)
axs[0].plot(np.arange(len(v_ks))/Fs,v_ks*Fs,label='est')
axs[0].plot([0,T],[f0,f1],label='1st chirp')
axs[0].plot([0,T],[f0+f_diff,f1+f_diff],label='2nd chirp')
axs[0].legend()
axs[1].plot(np.arange(NT),x0,label='true x0')
axs[1].plot(np.arange(len(Xs)),np.real(Xs)*2,label='est x0')
err=x0[:len(Xs)] - np.real(Xs)*2
axs[1].plot(np.arange(len(Xs)),err,label='err')
axs[1].legend()


plt.show()

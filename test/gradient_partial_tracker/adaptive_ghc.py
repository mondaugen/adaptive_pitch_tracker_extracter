# Example of adaptive hill climbing

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import dft_hill_climbing as dhc

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

v_ks,Xs,grad=dhc.adaptive_ghc_slow(x,f_k/Fs,w,mu=1e-6)

fig,axs=plt.subplots(3,1)
axs[0].plot(np.arange(len(v_ks))/Fs,v_ks*Fs,label='est')
axs[0].plot([0,T],[f0,f1],label='1st chirp')
axs[0].plot([0,T],[f0+f_diff,f1+f_diff],label='2nd chirp')
axs[0].legend()
axs[1].plot(np.arange(NT),x0,label='true x0')
axs[1].plot(np.arange(len(Xs)),np.real(Xs)*2,label='est x0')
err=x0[:len(Xs)] - np.real(Xs)*2
axs[1].plot(np.arange(len(Xs)),err,label='err')
axs[1].legend()
axs[2].plot(np.arange(len(v_ks))/Fs,grad,label='gradient')
axs[2].legend()


plt.show()

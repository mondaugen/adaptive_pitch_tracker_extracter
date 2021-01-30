# Examine the gradients of different windows

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import signal
from dft_hill_climbing import dft_bin_pow, dft_bin_pow_dv, dft_bin_pow_d2v

WINTYPE=os.environ.get('WINTYPE','hann')
v_min=-0.0002
v_max=0.0002
N_v=1000
N_w=4096
v=np.linspace(v_min,v_max,N_v)
x=signal.get_window(WINTYPE,N_w)
x/=np.sum(x)
X=np.zeros(N_v)
dX=np.zeros(N_v)
d2X=np.zeros(N_v)
for k, v_ in enumerate(v):
    X[k] = dft_bin_pow(x,v_)
    dX[k] = dft_bin_pow_dv(x,v_)
    d2X[k] = dft_bin_pow_d2v(x,v_)

fig,ax=plt.subplots(4,1)
ax[0].plot(v,10*np.log10(X),label='power spectrum 10*log10')
ax[1].plot(v,dX,label='grad power spectrum')
ax[2].plot(v,d2X,label='second derivative power spectrum')
ax[3].plot(v,v/dX,label='step/gradient ratio')
for a in ax:
    a.legend()

plt.show()
    

# Investigate the phase error when time stretching and then compressing so that
# the input and ouput signals match in length and from some time point match in
# content (or they should match).

import numpy as np
import classic_puckette_timestretch as cpt
import matplotlib.pyplot as plt
from scipy import signal

snd=np.fromfile('/tmp/render.f64')
snd+=np.random.standard_normal(len(snd))*1e-6
W=1024
H=256
#N=100000
N=len(snd)
SAMPLE_RATE=16000
D=2
n=np.arange(N)
f0=.1
f1=.2
#x=np.random.standard_normal(N)
#x=signal.chirp(n,f0,N,f1)
#no=np.random.standard_normal(N)
#x=(no>3).astype('float')+no[::-1]*1e-8
x=snd
t=np.arange(0,N,H)
# distrupt the analysis times
#t_d_idcs=np.arange(100,1100)
# distrupt from 0.125 to 0.375
t_d_idcs=np.where((t>(0.25*SAMPLE_RATE))&(t<((D-1.5)*SAMPLE_RATE)))[0]
#t[t_d_idcs]-=np.linspace(H,0,len(t_d_idcs)).astype('int')
# freeze
t[t_d_idcs]=t[t_d_idcs[0]]+np.linspace(0,50,len(t_d_idcs)).astype('int')
print(t[t_d_idcs])
# reset the pvoc after the disruption
#reset_times=[t[t_d_idcs[-1]+1]]
attack_time=D*SAMPLE_RATE
reset_times=[t[np.argmin(np.abs(t-attack_time))-2]]
y=cpt.time_stretch_arb_times(x, t, H, W, reset_times=reset_times)[:N]

fig,axs=plt.subplots(3,1)
axs[0].plot(n,x)
axs[0].set_title('original signal')
axs[1].plot(n,y)
axs[1].set_title('time stretched signal')
axs[1].plot([t[t_d_idcs[0]],t[t_d_idcs[0]]],[0,1],'r')
axs[1].plot([t[t_d_idcs[-1]],t[t_d_idcs[-1]]],[0,1],'r')
axs[1].plot([reset_times[0],reset_times[0]],[0,1],'g')
axs[2].plot(n,np.abs(x-y))
axs[2].set_title('error')

x.tofile('/tmp/x.f64')
y.tofile('/tmp/y.f64')

print(np.sum(np.abs(x-y)))


plt.tight_layout()
plt.show()

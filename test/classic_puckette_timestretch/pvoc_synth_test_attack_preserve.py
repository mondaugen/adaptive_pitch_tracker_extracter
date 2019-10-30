from classic_puckette_timestretch import pvoc_synth
import window_tools
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from time_map_tstretch import attack_avoider

W=1024
H=256
N=500*H
n=np.arange(N)
# chirp frequency
f0=0.01
# stretch factor
S=0.5

attack_times=np.array([1513,3013,4513])

x=signal.chirp(n,f0,N,f0)
x[attack_times]=1
av=attack_avoider(attack_times,-H,H+W,H)

output_times=np.round(np.arange(0,N,S*H)).astype('int')
y=np.zeros(len(output_times)*H,dtype=x.dtype)

wl=window_tools.windowed_lookup(x,W)
pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    lambda t,l: wl.access(t))

h_y=0
for h in output_times:
    atime,reset=av.adjust(h)
    y[h_y:h_y+H]=pv.process(atime,reset)
    h_y+=H

fig,axs=plt.subplots(2,1)

axs[0].plot(n,x)
axs[0].set_title('original')
axs[1].plot(np.arange(len(y)),y)
axs[0].set_title('stretched')

plt.show()


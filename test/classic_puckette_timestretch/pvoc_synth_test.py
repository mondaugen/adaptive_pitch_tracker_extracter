from classic_puckette_timestretch import pvoc_synth
import window_tools
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

W=1024
H=256
N=50*H
n=np.arange(N)
# chirp frequency
f0=0.05
# stretch factor
S=0.5

x=signal.chirp(n,f0,N,f0)

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
    y[h_y:h_y+H]=pv.process(h,False)
    h_y+=H

fig,axs=plt.subplots(2,1)

axs[0].plot(n,x)
axs[0].set_title('original')
axs[1].plot(np.arange(len(y)),y)
axs[0].set_title('stretched')

plt.show()


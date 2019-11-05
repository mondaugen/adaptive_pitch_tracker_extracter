# Simply pitch shift a sine wave

import numpy as np
from scipy import signal
from pitch_shift import pitch_shifter
import window_tools
from classic_puckette_timestretch import pvoc_synth
import matplotlib.pyplot as plt

W=1024
H=256
N=500*H
n=np.arange(N)
# chirp frequency
f0=0.01
# stretch factor
S=1.
# shift factor
P=0.5
x=signal.chirp(n,f0,N,f0)

wl=window_tools.windowed_lookup(x,W)
pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    lambda t,l: wl.access(t))

ps=pitch_shifter(
    lambda t,n: pv.process(int(np.round(t)),False),
    B=H)

y=np.zeros_like(x)
for h in range(0,N,H):
    y[h:h+H]=ps.process(np.ones(H)*P)

fig,axs=plt.subplots(2,1)

axs[0].plot(n,x)
axs[0].set_title('original')
axs[1].plot(np.arange(len(y)),y)
axs[1].set_title('pitch-shifted')

plt.show()

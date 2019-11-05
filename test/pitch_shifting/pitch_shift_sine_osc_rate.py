# Simply pitch shift a sine wave

import numpy as np
from scipy import signal
from pitch_shift import pitch_shifter
import window_tools
from classic_puckette_timestretch import pvoc_synth
import matplotlib.pyplot as plt

from_file=True

W=1024
H=256
if from_file:
    x=np.fromfile('/tmp/me.f64',dtype='float64')
    N=0
    while N < (len(x)-H):
        N+=H
    x=x[:N]
    n=np.arange(N)
    x+=np.random.standard_normal(N)*1e-8
else:
    N=500*H
    n=np.arange(N)
    # chirp frequency
    f0=0.01
    x=signal.chirp(n,f0,N,f0)
# stretch factor
S=1.
# shift factors
min_P=1.5
max_P=.5
P_osc_freq=0.000001
p=signal.chirp(n,P_osc_freq,N,P_osc_freq)
#p_=signal.chirp(n,P_osc_freq,N,P_osc_freq)
#p=np.zeros_like(p_)
#p[p_ > 0] = max_P
#p[p_ <= 0] = min_P

p=0.5*(p+1)
p*=(max_P-min_P)
p+=min_P

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
    y[h:h+H]=ps.process(p[h:h+H])

fig,axs=plt.subplots(2,1)

axs[0].plot(n,x)
axs[0].set_title('original')
axs[1].plot(np.arange(len(y)),y)
axs[1].set_title('pitch-shifted')

x.tofile('/tmp/x.f64')
y.tofile('/tmp/y.f64')

plt.show()

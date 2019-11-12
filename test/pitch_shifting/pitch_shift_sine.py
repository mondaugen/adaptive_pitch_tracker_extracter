# Simply pitch shift a sine wave

import numpy as np
from scipy import signal
from pitch_shift import pitch_shifter
import window_tools
from classic_puckette_timestretch import pvoc_synth
import matplotlib.pyplot as plt
import rel_del_line

REAL_TIME=False
from_file=True

W=1024
H=256

# stretch factor
S=1.
# shift factor
P=1.5

if from_file:
    x=np.fromfile('/tmp/in.f64')
    N=0
    while N < (len(x)-H):
        N += H
    x=x[:N]
    x+=np.random.standard_normal(N)*1e-8
    n=np.arange(N)
else:
    N=500*H
    n=np.arange(N)
    # chirp frequency
    f0=0.01
    f1=0.02
    x=signal.chirp(n,f0,N,f1)

class dummy_access(rel_del_line.access_struct):
    def __call__(self):
        pass

if REAL_TIME:
    L=H+W+1
    wl=rel_del_line.rel_del_line(L,print_warnings=True)
    wl.process(np.random.standard_normal(L)*1e-8)
    wl.time_zero_ago = 0
    def wl_access(t,l):
        af=dummy_access()
        wl.access(t,l,af)
        return af.values
    def ts_access(t,n):
        # we just get samples as close to the beginning of the delay line as possilbe
        t = wl.time_zero_ago-W
        return pv.process(t,False)
else:
    wl=window_tools.windowed_lookup(x,W)
    def wl_access(t,l):
        return wl.access(t)
    def ts_access(t,n):
        return pv.process(int(np.round(t)),False)

pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    wl_access)
#    lambda t,l: wl.access(t))

ps=pitch_shifter(ts_access,B=H)

y=np.zeros_like(x)
for h in range(0,N,H):
    if REAL_TIME:
        wl.process(x[h:h+H])
    y[h:h+H]=ps.process(np.ones(H)*P)

fig,axs=plt.subplots(2,1,sharex=True,sharey=True)

axs[0].plot(n,x)
axs[0].set_title('original')
axs[1].plot(np.arange(len(y)),y)
axs[1].set_title('pitch-shifted')

x.tofile('/tmp/orig.f64')
y.tofile('/tmp/out.f64')

plt.show()

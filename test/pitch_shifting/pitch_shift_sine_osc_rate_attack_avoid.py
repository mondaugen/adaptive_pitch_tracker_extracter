# Simply pitch shift a sine wave

import numpy as np
from scipy import signal
from pitch_shift import pitch_shifter
import window_tools
from classic_puckette_timestretch import pvoc_synth
import matplotlib.pyplot as plt
import rel_del_line
from time_map_tstretch import attack_avoider

REAL_TIME=False
from_file=False
adjust_for_attacks=True

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
    attack_times=np.array([100*H+10,200*H-10,300*H-20])
    n=np.arange(N)
    # chirp frequency
    f0=0.01
    x=signal.chirp(n,f0,N,f0)
    x[attack_times]=1
av=attack_avoider(attack_times,-H,H+W,H)
# stretch factor
S=1.
# shift factors
min_P=1.5
max_P=.5
P_osc_freq=0.00001
p=signal.chirp(n,P_osc_freq,N,P_osc_freq)

p=0.5*(p+1)
p*=(max_P-min_P)
p+=min_P

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
        if adjust_for_attacks:
            atime,reset=av.adjust(int(np.round(t)))
        else:
            atime=int(np.round(t))
            reset=False
        return pv.process(atime,reset)

pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    wl_access)

ps=pitch_shifter(ts_access,B=H)

y=np.zeros_like(x)
for h in range(0,N,H):
    if REAL_TIME:
        wl.process(x[h:h+H])
    y[h:h+H]=ps.process(p[h:h+H])

fig,axs=plt.subplots(2,1)

axs[0].plot(n,x)
axs[0].set_title('original')
axs[1].plot(np.arange(len(y)),y)
axs[1].set_title('pitch-shifted')

x.tofile('/tmp/x.f64')
y.tofile('/tmp/y.f64')

plt.show()

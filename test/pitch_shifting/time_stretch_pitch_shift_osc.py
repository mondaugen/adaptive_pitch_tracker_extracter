# Test that time stretching and pitch shifting can be done in tandem
import numpy as np
from scipy import signal
import pitch_shift
from classic_puckette_timestretch import pvoc_synth 
from window_tools import windowed_lookup
import matplotlib.pyplot as plt

W=2048
H=256
N0=5*W
N1=5*W
N=N0+N1
f0=500/16000
f1=600/16000
ps0=0.5
ps1=2
x=np.concatenate((
    signal.chirp(np.arange(W),f0,W,f0),
    np.zeros(N0-W),
    signal.chirp(np.arange(W),f1,W,f1),
    np.zeros(N1-W),
))
x+=np.random.standard_normal(N)*1e-6

wl=windowed_lookup(x,W)

pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    lambda n: wl.access(n))

ps=pitch_shift.pitch_shifter(lambda t: pv.process(int(np.round(t)),False))
y=np.zeros(N)
for n in np.arange(0,N0-W,H):
    y[n:n+H]=ps.process(np.ones(H)*ps0,np.ones(H)*ps0)
pv.reset_past_output()
ps.set_pos_at_block_start(N0)
for n in np.arange(N0,N-H,H):
    y[n:n+H]=ps.process(np.ones(H)*ps1,np.ones(H)*ps1)

x.tofile('/tmp/in.f64')
y.tofile('/tmp/out.f64')
n=np.arange(N)
plt.plot(n,x,label='in')
plt.plot(n,y,label='out')
plt.legend()
plt.show()

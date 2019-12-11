# Make a constant tone at some point
# time-stretch, then pitch-shift and then do both to see if the pitch-shifter is
# working

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import envelopes
import lfo
import pitch_shift
from classic_puckette_timestretch import pvoc_synth 
from window_tools import windowed_lookup
import matplotlib.pyplot as plt
import time_map_tstretch

show_plot=True
SR=16000
T=16
N=SR*T
n=np.arange(N)
t=n/SR
x=np.zeros(N)
T_start=6
N_start=int(np.round(T_start*SR))
T_end=12
N_end=int(np.round(T_end*SR))
L_tone=N_end-N_start
n_tone=np.arange(L_tone)
f0=1000
x[N_start:N_end]=signal.chirp(n_tone,f0/SR,L_tone,f0/SR)
x+=np.random.standard_normal(N)*1e-6
attack_times=[N_start,N_end]


# first time-stretch factor
ts_0=3/2
# first pitch-shift factor
ps_0=0.75

# prepare control signals
W=1024
H=128
T_out=T
N_out=int(np.ceil(SR*T_out))
y=np.zeros(N_out)

# activation gate
tone_out_T0=1
tone_out_N0=tone_out_T0*SR
tone_out_N1=N_out
gate_0=np.zeros_like(y)
gate_0[tone_out_N0:tone_out_N1]=1

rs=envelopes.region_segmenter(H)
wl=windowed_lookup(x,W)
pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    lambda n: wl.access(n))
av=time_map_tstretch.attack_avoider(
    attack_times,
    -3*H,
    W+2*H,
    H)
aaa=time_map_tstretch.attack_avoid_access(
    lambda t,r: pv.process(int(np.round(t)),r),
    av)
ps=pitch_shift.pitch_shifter(
aaa,
B=H)
adsr=envelopes.gate_to_adsr(200,1,1,200)

# control signals
pos_sig=np.linspace(3*SR,3*SR,N_out)
ps_sig=np.linspace(ps_0,ps_0,N_out)
ts_sig=np.linspace(ts_0,ts_0,N_out)

for n_ in range(0,N_out,H):
    H_=min(N_out-n_,H)
    en,st,ed,ac,sta=adsr.gate_to_adsr_env_start_end_active(gate_0[n_:n_+H_])
    ant,ans,ane,regs=rs.region_segmenter_update(st,ed,ac)
    for s,e in regs:
        if e>s:
            l=e-s
            if st[s] > 0:
                pv.reset_past_output()
                ps.set_pos_at_block_start(pos_sig[n_+s])
            y[n_+s:n_+e]=ps.process(ps_sig[n_+s:n_+e],ts_sig[n_+s:n_+e])*en['adsr'][s:e]

plt.plot(t,x,label='input')
plt.plot(t,y+2,label='output')
plt.legend(loc='upper right')
x.tofile('/tmp/x.f64')
y.tofile('/tmp/y.f64')
if show_plot:
    plt.show()

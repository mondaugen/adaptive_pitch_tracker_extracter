# To test to see how well time-stretching and pitch-shifting works together we:
# Synthesize a chirp
# Then we synthesize a signal which is synthesized from parts of the chirp, but
# is offset, pitch-shifted and time-stretched
# we can't compensate perfectly for the chirp because we have to take more or
# fewer samples in order to give enough points for the interpolation. For
# example, if we are trying to pitch shift down at exactly the inverse rate of
# the chirp (so that the chirp remains at its initial pitch), we compensate the
# pitch-shifting down (over-sampling) by time-compressing. By time-compressing,
# we move further along the chirp, which in the case of an increasing chirp
# would be higher than our pitching-down, so we fail to pitch shift it down
# enough.
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import envelopes
import lfo
import pitch_shift
from classic_puckette_timestretch import pvoc_synth 
from window_tools import windowed_lookup
import matplotlib.pyplot as plt

W=1024
H=128
SR=16000
T=8
N=T*SR
f0=100/SR
f1=200/SR
n=np.arange(N)
x=signal.chirp(n,f0,N,f1)
x+=np.random.standard_normal(N)*1e-6

rs=envelopes.region_segmenter(H)

wl=windowed_lookup(x,W)

pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    lambda n: wl.access(n))

ps=pitch_shift.pitch_shifter(
lambda t: pv.process(int(np.round(t)),False),
B=H)

gate=lfo.chirp(N,1/SR/4,1/SR/4,phase=-0.25).squarewave()
adsr=envelopes.gate_to_adsr(20,1,1,20)
env=np.zeros_like(gate)
starts=np.zeros_like(gate)
ends=np.zeros_like(gate)
active=np.zeros_like(gate)
# a signal from which we get the start position, sampled when an envelope region starts
pos=np.linspace(0,N,N)
# a signal from which we get the pitch shift
ps_sig=np.linspace(1,1,N)#f0/np.linspace(f0,f1,N+1)[:N]
# a signal from which we get the time_stretch
ts_sig=np.linspace(1,0,N)
y=np.zeros_like(ts_sig)
for n_ in range(0,N,H):
    H_=min(N-n_,H)
    en,st,ed,ac,sta=adsr.gate_to_adsr_env_start_end_active(gate[n_:n_+H_])
    env[n_:n_+H_]=en['adsr']
    starts[n_:n_+H_]=st
    ends[n_:n_+H_]=ed
    active[n_:n_+H_]=ac
    ant,ans,ane,regs=rs.region_segmenter_update(st,ed,ac)
    for s,e in regs:
        if e>s:
            l=e-s
            if st[s] > 0:
                pv.reset_past_output()
                ps.set_pos_at_block_start(pos[n_+s])
            y[n_+s:n_+e]=ps.process(ps_sig[n_+s:n_+e],ts_sig[n_+s:n_+e])*en['adsr'][s:e]

plt.figure()
plt.plot(n,gate,label='gate')
plt.plot(n,env+1,label='env')
plt.plot(n,starts,label='starts')
plt.plot(n,ends,label='ends')
plt.plot(n,active+2,label='active')
plt.legend()

plt.figure()
plt.plot(n,y)

plt.figure()
plt.specgram(x,NFFT=W,Fs=SR,noverlap=W-H)

plt.figure()
plt.specgram(y,NFFT=W,Fs=SR,noverlap=W-H)

y.tofile('/tmp/out.f64')

plt.show()

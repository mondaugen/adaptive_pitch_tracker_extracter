# Make a sound consisting of two notes that are sine waves windowed by an ADSR
# Then make a resulting sound by playing these notes at different transpositions
# and speeds (time-stretch).

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import envelopes
import lfo
import pitch_shift
from classic_puckette_timestretch import pvoc_synth 
from window_tools import windowed_lookup
import matplotlib.pyplot as plt

# prepare input signal

SR=16000
note_params=[
    dict(T=2,f=440),
    dict(T=3,f=660)
]

for note in note_params:
    note['N']=note['T']*SR
    note['v']=note['f']/SR
    note['x']=signal.chirp(np.arange(note['N']),note['v'],note['N'],note['v'])
    note['g']=lfo.chirp(note['N'],1/note['N'],1/note['N'],phase=-0.25).squarewave()
    adsr=envelopes.gate_to_adsr(100,1000,0.5,4000)
    en,st,ed,ac,sta=adsr.gate_to_adsr_env_start_end_active(note['g'])
    note['e']=en['adsr']
    note['y']=note['x']*note['e']

N=0
for note in note_params:
    N+=note['N']
y=np.zeros(N)
n=0
for note in note_params:
    y[n:n+note['N']]=note['y']
    n+=note['N']
y+=np.random.standard_normal(N)*1e-6

y.tofile('/tmp/in.f64')

# prepare control signals
W=1024
H=128
T_out=30
N_out=SR*T_out
gate_f0=2/SR
gate_sig=lfo.chirp(N_out,gate_f0,gate_f0,phase=-0.3).squarewave()
pos_sig=lfo.chirp(
    N_out,
    gate_f0*0.5,
    gate_f0*0.5,
    phase=-0.25,
    x_min=0,
    x_max=note_params[0]['N']).squarewave()
ps_sig=np.linspace(1,1,N_out)
ts_sig=np.linspace(4,1,N_out)

rs=envelopes.region_segmenter(H)

wl=windowed_lookup(y,W)

pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    lambda n: wl.access(n))

ps=pitch_shift.pitch_shifter(
lambda t: pv.process(int(np.round(t)),False),
B=H)
adsr=envelopes.gate_to_adsr(200,1,1,200)

z=np.zeros_like(ts_sig)
for n_ in range(0,N_out,H):
    H_=min(N_out-n_,H)
    en,st,ed,ac,sta=adsr.gate_to_adsr_env_start_end_active(gate_sig[n_:n_+H_])
    ant,ans,ane,regs=rs.region_segmenter_update(st,ed,ac)
    for s,e in regs:
        if e>s:
            l=e-s
            if st[s] > 0:
                pv.reset_past_output()
                ps.set_pos_at_block_start(pos_sig[n_+s])
            z[n_+s:n_+e]=ps.process(ps_sig[n_+s:n_+e],ts_sig[n_+s:n_+e])*en['adsr'][s:e]

plt.figure(0)
plt.plot(y)

plt.figure(1)
plt.plot(z)

z.tofile('/tmp/out.f64')
plt.show()

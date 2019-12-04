# Simply gate a synthesized square wave signal

from wavetables import bl_square_synth
import numpy as np
import lfo
import envelopes
import matplotlib.pyplot as plt

# Sample rate
SR=16000
# Length of signal in seconds
T=60
# Block size
B=256
# Gate frequency
gate_f0=5
gate_f1=7.5
# Length in samples
N=int(np.round(T*SR))
# frequency of synthesizer (Hz)
f0=10
f1=70
fm_f0=4.5

blss=bl_square_synth()
f=lfo.chirp(N,fm_f0/SR,fm_f0/SR,x_min=f0/SR,x_max=f1/SR).sinusoid()

gate=lfo.chirp(N,gate_f0/SR,gate_f1/SR).squarewave()
adsr=envelopes.gate_to_adsr(SR*0.01, SR*0.01, 1, SR*0.01)
x=np.zeros(N)
rs=envelopes.region_segmenter(B)
for n in np.arange(0,N-B,B):
    (adsr_dict,
    start,
    end,
    state,
    adsr_states)=adsr.gate_to_adsr_env_start_end_active(gate[n:n+B])
    ant,ans,ane,regs=rs.region_segmenter_update(start,end,state)
    for s,e in regs:
        if e>s:
            l=e-s
            x[n+s:n+e]=blss.synth(f[n+s:n+e],start[s]>0)*adsr_dict['adsr'][s:e]

plt.plot(np.arange(N),x)
plt.show()

x.tofile('/tmp/adsr_region_synth.f64')

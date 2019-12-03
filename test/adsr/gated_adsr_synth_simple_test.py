# Simply gate a synthesized square wave signal

from wavetables import bl_square_synth
import numpy as np
import lfo
import envelopes

# Sample rate
SR=16000
# Length of signal in seconds
T=10
# Block size
B=256
# Gate frequency
gate_f0=0.5
# Length in samples
N=int(np.round(T*SR))

blss=bl_square_synth()

gate=lfo.chirp(N,gate_f0/SR,gate_f0/SR).squarewave()
adsr=envelopes.gate_to_adsr(SR*0.1, SR*0.1, 0.5, SR*1)
x=np.zeros(N)
for n in np.arange(0,N-B,B):
    a=adsr.gate_to_adsr_seq(gate[n:n+B])
    e=adsr.adsr_seq_to_env(a)
    r=(e!>0).astype('float')
    x[n:n+B]=
    

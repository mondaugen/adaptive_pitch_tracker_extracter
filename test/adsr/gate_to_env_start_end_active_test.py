import numpy as np
import matplotlib.pyplot as plt
import common
import envelopes
import lfo

# Sample rate
SR=16000
# Length of signal in seconds
T=3
# Block size
B=256
# Gate frequency
gate_f0=0.5
# Length in samples
N=int(np.round(T*SR))
env=np.zeros(N)
start=np.zeros(N)
end=np.zeros(N)
state=np.zeros(N)

gate=lfo.chirp(N,gate_f0/SR,gate_f0/SR).squarewave()
adsr=envelopes.gate_to_adsr(SR*0.1, SR*0.1, 0.5, SR*.1)

for n in np.arange(0,N-B,B):
    (env[n:n+B],
    start[n:n+B],
    end[n:n+B],
    state[n:n+B])=adsr.gate_to_adsr_env_start_end_active(gate[n:n+B])

n=np.arange(N)

plt.plot(n,env)
common.logic_plot(n,start+1)
common.logic_plot(n,end+2)
common.logic_plot(n,state+3)
common.logic_plot(n,gate+4)

plt.show()


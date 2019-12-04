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
gate_f0=.5
# Length in samples
N=int(np.round(T*SR))
env=np.zeros(N)
start=np.zeros(N)
end=np.zeros(N)
state=np.zeros(N)
adsr_states=np.zeros(N)
attack=np.zeros(N)
decay=np.zeros(N)
sus=np.zeros(N)
rel=np.zeros(N)

gate=lfo.chirp(N,gate_f0/SR,gate_f0/SR).squarewave()
adsr=envelopes.gate_to_adsr(SR*0.1, SR*0.1, 0.5, SR*.1)

for n in np.arange(0,N-B,B):
    (adsr_dict,
    start[n:n+B],
    end[n:n+B],
    state[n:n+B],
    adsr_states[n:n+B])=adsr.gate_to_adsr_env_start_end_active(gate[n:n+B])
    env[n:n+B]=adsr_dict['adsr']
    attack[n:n+B]=adsr_dict['a_ramp']
    decay[n:n+B]=adsr_dict['d_sig']
    sus[n:n+B]=adsr_dict['s_sig']
    rel[n:n+B]=adsr_dict['r_sig']
#(env,
#start,
#end,
#state)=adsr.gate_to_adsr_env_start_end_active(gate)

n=np.arange(N)

plt.plot(n,env,label='envelope')
common.logic_plot(n,start+1,label='starts')
common.logic_plot(n,end+2,label='ends')
common.logic_plot(n,state+3,label='state')
common.logic_plot(n,gate+4,label='gate')
common.logic_plot(n,adsr_states+5,label='adsr')
plt.plot(n,attack+9,label='attack')
plt.plot(n,decay+10,label='decay')
plt.plot(n,sus+11,label='sustain')
plt.plot(n,rel+12,label='release')
plt.legend()

plt.show()


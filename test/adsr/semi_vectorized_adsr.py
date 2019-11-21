import numpy as np
import matplotlib.pyplot as plt
from envelopes import gate_to_adsr
from common import logic_plot

N=1000
thresh=1
gate=np.random.standard_normal(N)
gate=(gate>thresh).astype('float')
gate=np.cumsum(gate)
gate %= 2
#gate=np.zeros(N)
#gate[25:N-25]=1
gate_to_adsr_inst=gate_to_adsr(3,4,0.5,5,decay_min_dB=-24)
adsr_sig=gate_to_adsr_inst.gate_to_adsr_seq(gate)
adsr_env=gate_to_adsr_inst.adsr_seq_to_env(adsr_sig)
logic_plot(np.arange(N),gate,label='gate')
logic_plot(np.arange(N),adsr_sig,color='purple',label='ADSR states')
logic_plot(np.arange(N),adsr_env,color='orange',label='ADSR envelope')
plt.legend()

plt.show()


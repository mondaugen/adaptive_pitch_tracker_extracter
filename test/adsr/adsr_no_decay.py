import numpy as np
import matplotlib.pyplot as plt
import envelopes
import common

N=30
n=np.arange(N)
adsr=envelopes.gate_to_adsr(5,1,1,5)

gate0=np.zeros(N)
gate0[5:20]=1
states0=adsr.gate_to_adsr_seq(gate0)

gate1=np.zeros(N)
gate1[5:8]=1
gate1[9:19]=1
states1=adsr.gate_to_adsr_seq(gate1)

plt.figure()
common.logic_plot(n,states0)
plt.figure()
common.logic_plot(n,states1)

plt.show()

import envelopes
import numpy as np
import matplotlib.pyplot as plt

N=100
n=np.arange(N)
gtta=envelopes.gate_to_time_advance(5,10,0.25,0.1,0.2)
gate=np.zeros(N)
gate[5:15]=1
gate[18:30]=1
gate[41:43]=1
gate[55:75]=1
time=gtta.gate_to_time_sequence(gate)

plt.plot(n,time)
plt.show()

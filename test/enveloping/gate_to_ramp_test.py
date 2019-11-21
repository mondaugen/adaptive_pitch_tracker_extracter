import numpy as np
import matplotlib.pyplot as plt
from envelopes import gate_to_ramp
from common import logic_plot

N=100
thresh=1.
gate=np.random.standard_normal(N)
gate=(gate>1.).astype('float')
gate=np.cumsum(gate)
gate %= 2
ramp_sig=gate_to_ramp(gate)
assert(len(gate)==len(ramp_sig))
logic_plot(np.arange(N),gate)
logic_plot(np.arange(N),ramp_sig,color='purple')

plt.show()


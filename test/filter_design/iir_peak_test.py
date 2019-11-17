import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

f0=0.03
f1=0.04
q0=100
q1=10
w=np.arange(0,0.5,0.001)
b0,a0=signal.iirpeak(f0,q0,fs=1)
b1,a1=signal.iirpeak(f1,q1,fs=1)
w0,h0=signal.freqz(b0,a0,worN=w,fs=1)
w1,h1=signal.freqz(b1,a1,worN=w,fs=1)

plt.plot(w0,20*np.log10(np.abs(h0)),label='f0')
plt.plot(w1,20*np.log10(np.abs(h1)),label='f1')
plt.legend()

plt.show()

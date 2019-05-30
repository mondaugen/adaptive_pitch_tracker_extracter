import ptracker
import matplotlib.pyplot as plt
import numpy as np

# this shows not beating but ringing, which happens when the signal starts
# perhaps why the pitch estimate is delayed and oscillating

sr=16000
f=100
f_=101
T=100000
t=np.arange(T)
comb=ptracker.pitch_check_comb(f/sr,.99)
x=np.cos(f_/sr*t*2*np.pi)
y=comb.proc(x)

plt.plot(t,y)
plt.show()

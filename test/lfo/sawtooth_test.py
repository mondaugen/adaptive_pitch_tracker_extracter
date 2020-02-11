import numpy as np
import lfo
import matplotlib.pyplot as plt

N=10000
n=np.arange(N)
f0=1e-3
f1=1e-3
xmin=-0.2
xmax=0.5
phase = 0.25
c=lfo.chirp(N,f0,f1,x_min=xmin,x_max=xmax,phase=phase)
rup=c.sawtooth()
rdown=c.sawtooth(falling=True)
plt.plot(n,rup,label='up')
plt.plot(n,rdown,label='down')
plt.legend()
plt.show()

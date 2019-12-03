import matplotlib.pyplot as plt
from scipy import signal
import numpy as np

N=1000
n=np.arange(N)
f0=0.01
x0=signal.chirp(n,f0,N,f0,phi=0)
x1=signal.chirp(n,f0,N,f0,phi=180)

plt.plot(n,x0)
plt.plot(n,x1)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import lfo

N=1000
ch=lfo.chirp(N,0.01,0.01)
x_si=ch.sinusoid()
x_sq=ch.squarewave()
n=np.arange(N)

plt.plot(n,x_si)
plt.plot(n,x_sq)

plt.show()

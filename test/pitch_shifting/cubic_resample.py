from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt

M=2

for m in range(M):
    N=30

    n=np.arange(N)
    x=np.random.standard_normal(N)

    t=np.arange(0,N-2,0.2)

    x_cubic_int=interpolate.interp1d(n,x,kind='cubic')

    x_int=x_cubic_int(t)

    plt.plot(n,x,'.')
    plt.plot(t,x_int)

plt.show()

import numpy as np
import matplotlib.pyplot as plt
import vectorized_cubic

N=2**15
P=2500
i=1/16
x=np.arange(N)
y=np.cos(2*np.pi*x/P).astype('float32')
xi=np.arange(1,N-2,i).astype('float32')
xr=np.round(xi).astype('uint32')
yi=vectorized_cubic.lookup4(xi,y)
err=yi-y[xr]
err=20*np.log10(np.abs(err))

fig,ax=plt.subplots(2,1)
ax[0].plot(x,y,label='original')
ax[0].plot(xi,yi,label='interpolated')
ax[1].plot(xr,err,label='err dB')

plt.show()

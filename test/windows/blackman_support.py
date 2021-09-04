import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

lobe_radius=2
oversamp=16
N=2048
w=np.zeros(N*oversamp)
w[(N*oversamp)//2-N//2:(N*oversamp)//2+N//2]=signal.get_window('hann',N)
w/=np.sum(w)
W=np.fft.fft(w)
v=np.arange(0,N,1/oversamp)
fig,ax=plt.subplots(1)
ax.plot(v,20*np.log(np.abs(W)))
ax.set_xlim([0,3])
print(np.abs(W)[lobe_radius*oversamp])
plt.show()

import attack_finder
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np

b,a=attack_finder._shifted_lp(N=16)
w,h=signal.freqz(b,a,whole=True,fs=2)

fig,ax=plt.subplots(2,1)
ax[0].plot(w,20*np.log10(np.abs(h)))
ax[1].plot(w,np.unwrap(np.angle(h)))
plt.show()

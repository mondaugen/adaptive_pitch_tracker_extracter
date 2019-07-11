import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import rab_pitch

w=signal.get_window('hann',16)
x=np.zeros((32,4))
x[8:24,:]=np.random.standard_normal((16,4))*w[:,None]
X=20*np.log10(np.abs(np.fft.rfft(x,axis=0))+1e-6)
for row in X.T:
    plt.plot(np.arange(len(row)),row)
X_peaks=rab_pitch.pick_peaks(X)
print(X_peaks)
for xi,yi in zip(*X_peaks):
    plt.plot(xi,X[xi,yi],'.r')

plt.show()


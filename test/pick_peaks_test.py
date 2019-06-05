import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import dft_freq_est

w=signal.get_window('hann',16)
x=np.zeros(32)
x[8:24]=np.random.standard_normal(16)*w
X=20*np.log10(np.abs(np.fft.rfft(x)))
plt.plot(np.arange(len(X)),X)
X_peaks=dft_freq_est.pick_peaks(X)
plt.plot(X_peaks,X[X_peaks],'.')

plt.show()

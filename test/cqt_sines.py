import cqt
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import common

filename="/tmp/sines.f64"
sr=16000

x=np.fromfile(filename,dtype='float64')
C=cqt.cqt(
lambda N: signal.get_window("hann",N),
1024,
ws=[2*np.pi*440*np.power(2,(n-69)/12)/sr for n in range(60,73)])
x=common.frame(x,256,1024)
X=C(x)

plt.imshow(20*np.log10(np.abs(X)),origin='lower',aspect='auto')

plt.show()

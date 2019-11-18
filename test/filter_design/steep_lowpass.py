import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

filter_N=2
b,a=signal.butter(filter_N,0.1,fs=1)
w,h=signal.freqz(b,a)
l,=plt.plot(w/(2*np.pi),20*np.log10(np.abs(h)))
filter_N=4
b,a=signal.butter(filter_N,0.001,fs=1)
w,h=signal.freqz(b,a)
l,=plt.plot(w/(2*np.pi),20*np.log10(np.abs(h)))
filter_N=4
b,a=signal.cheby1(filter_N,1,0.1,fs=1)
w,h=signal.freqz(b,a)
l,=plt.plot(w/(2*np.pi),20*np.log10(np.abs(h)))
l._axes.set_ylim(-100,5)
plt.show()

import attack_finder
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import librosa

x,sr=librosa.load('sounds/midi.wav',sr=16000)

dsf=1
xds=signal.resample_poly(x,1,dsf)
y=attack_finder.analytic_signal(xds)
t=np.arange(len(xds))/sr

y_unwrap=np.unwrap(np.angle(y))
y_diff=np.diff(y_unwrap)

fig,axs=plt.subplots(3,1)
axs[0].plot(t,np.real(y))
axs[0].set_title("original")
axs[1].plot(t[:-1],y_diff)
axs[1].set_title("instantaneous frequency of analytic signal")
axs[2].plot(t,y_unwrap)
axs[2].set_title("unwraped phase of analytic signal")

plt.show()

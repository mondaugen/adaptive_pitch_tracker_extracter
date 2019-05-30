import matplotlib.pyplot as plt
from scipy import signal
import ptracker
import librosa
import numpy as np

x,sr=librosa.load('sounds/guitar.wav',sr=16000)
x_=signal.decimate(x,10)
t=np.arange(len(x))/sr
t_=np.arange(len(x_))/(sr/10)

p=ptracker.avg_filter(x_*x_,alph=0.99)
dp=ptracker.avg_filter(ptracker.diff_filter(p,1),alph=0.9)

lp=10*np.log10(p)

fig,ax=plt.subplots(3,1)


ax[0].plot(t,x)
ax[1].plot(t_,lp)
ax[2].plot(t_,dp)

plt.show()

import attack_finder
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import librosa

x,sr=librosa.load('sounds/midi.wav',sr=16000)

y=attack_finder.analytic_signal(x)
t=np.arange(len(x))/sr
dsf=1
#xds=signal.decimate(x,dsf,n=8)
xds=signal.resample_poly(x,1,dsf)
sr_ds=sr/dsf
t_ds=np.arange(len(xds))/sr_ds

fig,ax=plt.subplots(4,1)
ax[0].plot(t,np.abs(y))
ax[1].plot(t,np.abs(x))
#xm=attack_finder.local_max(np.abs(xds),N=30)
xm=attack_finder.fast_max(np.abs(xds),alph=0.9999)
xmp=attack_finder.thresh_local_max_samples(np.diff(xm),alph=0.9999,beta=1,N=sr_ds*.125)
#xmp=np.diff(xm)
t_=np.arange(len(xm))/sr_ds
#ax[2].plot(t_,xm)
ax[2].plot(t_[:-1],np.diff(xm))
ax[3].plot(t_[:-1],xmp)

plt.show()

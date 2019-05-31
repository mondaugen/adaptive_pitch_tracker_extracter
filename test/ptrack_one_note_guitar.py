#from ptracking import ptrackers
import ptracker
import numpy as np
import matplotlib.pyplot as plt
import librosa

sr=16000
x,sr=librosa.load("sounds/one_guitar_note.wav",sr)
t=np.arange(len(x))/sr
w0=65.406/sr*2*np.pi # frequency of c2
n_parts=16

parts=ptracker.multi_ptrackers(
x,
w0,
n_parts)

t=np.arange(len(x))/sr

fig,axs=plt.subplots(n_parts,2)
for n,(a,phi,w) in enumerate(parts):
    axs[n,0].plot(t,w)
    axs[n,1].plot(t,a)

plt.show()
    



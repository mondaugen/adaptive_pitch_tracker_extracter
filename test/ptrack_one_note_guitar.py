#from ptracking import ptrackers
import ptracker
import numpy as np
import matplotlib.pyplot as plt
import librosa
from scipy import signal

sr=16000
x=np.fromfile('/tmp/one_guitar_note.f64',dtype='float64')
t=np.arange(len(x))/sr
w0=65.406/sr*2*np.pi # frequency of c2
n_parts=32

parts=ptracker.multi_ptrackers(
x,
w0,
n_parts,
ptrack_args=(
0.9999,
0.9,
1e-3,
1e-4)
)

alph_w=0.999
alph_a=0.99

t=np.arange(len(x))/sr

fig,axs=plt.subplots(n_parts,2)
y=np.zeros(t.shape,dtype='float64')
for n,(a,phi,w) in enumerate(parts):
    axs[n,0].plot(t,w)
    axs[n,1].plot(t,a)
    w_mean=np.mean(w)
    w_=w-w_mean
    w_=signal.lfilter([1-alph_w],[1,-alph_w],w_)
    a_=signal.lfilter([1-alph_a],[1,-alph_a],a)
    w_+=w_mean
    #w_=w
    a_=a
    phi_=np.cumsum(w_)
    y+=a_*np.cos(phi_)

y-=np.mean(y)
y/=np.max(np.abs(y))
y*=0.9
y.tofile('/tmp/synth_guit_nparts=%d.f64' % (n_parts,))

plt.show()
    



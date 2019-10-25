# You should be able to mix frames from the phase vocoder and directly from the input, no?

import numpy as np
from scipy import signal
import common

N=10000
W=32
H=8

# generate a signal
x=np.random.standard_normal(N,)
# extend so the windowing gives a full number of overlaps for the beginning and
# end when analysing at the output rate
x_ext=np.concatenate((np.zeros(W),x,np.zeros(W),np.zeros(H-(2*W+len(x))%H)))
print((len(x_ext)-W)/H)
# generate the analysis and synthesis window
w=signal.get_window('blackmanharris',W)
# compute the output scaling
win_div=1./common.ola_shorten(np.power(w,2),H)
x_ext_fr=common.frame(x_ext,H,W)*w[:,None]
X_ext=np.fft.rfft(x_ext_fr,axis=0)
y_from_X_fr=np.fft.irfft(X_ext,axis=0)

# build a signal that is a mix of the fourier transformed and the raw signal
y_from_X_and_x_fr=np.zeros_like(y_from_X_fr)
idcs=np.arange((len(x_ext)-W)//H)
np.random.shuffle(idcs)
idcs_X=idcs[:len(idcs)//2]
idcs_x=idcs[len(idcs)//2:]
y_from_X_and_x_fr[:,idcs_X]=y_from_X_fr[:,idcs_X]
y_from_X_and_x_fr[:,idcs_x]=x_ext_fr[:,idcs_x]

y_from_X=np.zeros_like(x_ext)
for n,h in enumerate(np.arange(0,len(x_ext)-W,H)):
    y_from_X[h:h+W]+=y_from_X_fr[:,n]*w
    y_from_X[h:h+H]*=win_div
print(np.sum(x-y_from_X[W:W+len(x)]))

y_from_x=np.zeros_like(x_ext)
for n,h in enumerate(np.arange(0,len(x_ext)-W,H)):
    y_from_x[h:h+W]+=x_ext_fr[:,n]*w
    y_from_x[h:h+H]*=win_div
print(np.sum(x-y_from_x[W:W+len(x)]))

y_from_X_and_x=np.zeros_like(x_ext)
for n,h in enumerate(np.arange(0,len(x_ext)-W,H)):
    y_from_X_and_x[h:h+W]+=y_from_X_and_x_fr[:,n]*w
    y_from_X_and_x[h:h+H]*=win_div
print(np.sum(x-y_from_X_and_x[W:W+len(x)]))

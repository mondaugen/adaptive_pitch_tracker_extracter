# This actually doesn't really help

import cqt
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import common

ge=common.get_env

filename=ge("FILENAME","/tmp/guit.f64",None)
min_pitch=36
max_pitch=96
sr=16000

x=np.fromfile(filename,dtype='float64')
C=cqt.cqt(
lambda N: signal.get_window("hann",N),
2048,
ws=[2*np.pi*440*np.power(2,(n-69)/12)/sr
    for n in range(min_pitch,max_pitch+1)])
x=common.frame(x,256,2048)
X=C(x)
X=np.abs(X)

# Weight by adding octave equivalence
W=np.zeros_like(X)
for row in range(W.shape[0]):
    n = row
    cnt=0
    while n < W.shape[0]:
        W[row,:]+=X[n,:]
        n+=12
        cnt+=1
    #if cnt > 0:
    #    W[row,:]/=cnt
X/=np.max(X)
W/=np.max(W)
plt.imshow(20*np.log10(W*X),origin='lower',aspect='auto')


plt.show()


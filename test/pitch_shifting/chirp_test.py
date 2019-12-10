import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

SR=16000
T=3
N=SR*T
f0=100/SR
f1=500/SR
n=np.arange(N)
x=signal.chirp(n,f0,N,f1)
W=1024
H=256
plt.specgram(x,NFFT=W,Fs=SR,noverlap=W-H)
plt.show()

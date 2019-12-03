import bl_square
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

show_plot=False

SR=16000
N=1000000
n=np.arange(N)
B=1024
blss=bl_square.bl_square_synth(max_f=0.25,N=8192,D=0.5)
print(blss.min_f)
fm0=.25/SR
fm1=.25/SR
f0=blss.min_f*4
f1=blss.min_f*4
f=signal.chirp(n,fm0,N,fm1)
f+=1
f*=0.5
f*=f1-f0
f+=f0
y=np.zeros(N).astype('float32')
for b in np.arange(0,N-B,B):
    y[b:b+B]=blss.synth(f[b:b+B])
y.tofile('/tmp/synth.f32')

if show_plot:
    plt.figure(0)
    plt.plot(n,y)
    plt.figure(1)
    plt.specgram(y,NFFT=2**14,Fs=SR)
    plt.show()

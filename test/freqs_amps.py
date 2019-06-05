import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import librosa
import dft_freq_est

sr=16000
x=np.fromfile('/tmp/snd.f64')

n_samps=min(len(x),10000)
x=x[:n_samps]

pa=dft_freq_est.peak_analyzer(
2048,
1024,
1)

fig,ax=plt.subplots(2,1)
w=[]
for n in np.arange(0,n_samps,1):
    x_=x[n:n+pa.N_W+pa.H]
    ph,a,w_,bet=pa.freqs_amps(x_)
    w.append(w_)
    t=[n/sr for _ in w_]
    ax[0].plot(t,w_,'.')

ax[1].plot(np.arange(n_samps)/sr,x)

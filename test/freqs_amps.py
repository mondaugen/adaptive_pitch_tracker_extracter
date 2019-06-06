import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import librosa
import dft_freq_est
import attack_finder

sr=16000
x=np.fromfile('/tmp/snd.f64')

N_W=1024
N_FFT=2048
N_H=256

peak_thresh=-100

n_samps=min(len(x),1e8)
n_samps=((n_samps-N_W)//N_H)*N_H
print(n_samps)

pa=dft_freq_est.peak_analyzer(
N_FFT,
N_W,
N_H)

fig,ax=plt.subplots(2,1)
w=[]
for n in np.arange(0,n_samps,pa.H):
    x_=x[n:n+pa.N_W+pa.H]
    ph,a,w_,bet=pa.freqs_amps(x_,max_n_peaks=64,peak_thresh=peak_thresh)
    w.append(w_)
    t=[(n+(pa.N_W+pa.H)/2)/sr for _ in w_]
    a_=(20*np.log10(a)-peak_thresh)/(-peak_thresh)
    ax[0].scatter(t,w_,s=1,c=a_)

ax[1].plot(np.arange(n_samps)/sr,x[:n_samps])
attack_i=attack_finder.find_attacks(x)
ax[1].scatter(attack_i/sr,x[attack_i],c='k')

plt.show()

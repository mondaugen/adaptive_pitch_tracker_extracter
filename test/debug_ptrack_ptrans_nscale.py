# why can't we transpose with partial tracking?

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import signal
import librosa
import dft_freq_est
import attack_finder
from common import normalize

def cat_to_cols(x,y):
    return np.concatenate((x[:,None],y[:,None]),axis=1)

do_plot=False

sr=16000
x=np.fromfile('/tmp/snd.f64')
x=normalize(x)

N_W=1024
N_FFT=2048
N_H=256
n_max=((len(x)-N_W)//N_H)*N_H

peak_thresh=-140

# intialize the sinusoidal analyser
pa=dft_freq_est.peak_analyzer(
N_FFT,
N_W,
N_H)

# range of peaks we are looking for
pitch_f0=36
f0=440*np.power(2,(pitch_f0-69)/12.)
n_parts=16
part_idx=15
max_detune_cents=200
# amount of transposition when resynthesizing
resynth_ptrans=0.5#np.power(2,-7./12)
print(resynth_ptrans)
f0s=(np.arange(n_parts)+1)*f0
#print(f0s)
w0s=f0s/sr*2*np.pi
w_ranges=np.zeros((n_parts,2))
w_ranges[:,0]=w0s*np.power(2,-max_detune_cents/1200.)
w_ranges[:,1]=w0s*np.power(2,max_detune_cents/1200.)

res=[]
for n in np.arange(0,n_max,N_H):
    x_=x[n:n+N_W+N_H]
    X_=pa.dft(x[n:n+N_W])
    ph,a,dph,da=pa.nearest_peaks(x_,w_ranges,default_a=np.power(10,peak_thresh/20))
    #print('freqs',dph)
    #print('amps',a)
    res.append((X_,ph,a,dph,da))

# synthesize results
(Xk,phk,ak,dphk,dak)=res[10]
(Xk1,phk1,ak1,dphk1,dak1)=res[11]
ph=cat_to_cols(phk,phk1)
a=cat_to_cols(ak,ak1)
dph=cat_to_cols(dphk,dphk1)
da=cat_to_cols(dak,dak1)
_,_,_,_,amp,th=pa.synth_peaks_trans_stretch(ph,a,dph,da,ptrans=1)
_,_,_,_,amp_trans,th_trans=pa.synth_peaks_trans_stretch(ph,a,dph,da,ptrans=resynth_ptrans)
n=np.arange(pa.H+1)
fig,ax=plt.subplots(3,1)
ax[0].plot(n,th[part_idx,:],label='orig')
ax[0].plot(n,th_trans[part_idx,:],label='trans')
ax[0].legend()
ax[1].plot(n[:-1],np.diff(th[part_idx,:]),label='orig')
ax[1].plot(n[:-1],np.diff(th_trans[part_idx,:]),label='trans')
ax[1].legend()
xdata=np.arange(N_FFT/2+1)*sr/N_FFT
ax[2].plot(xdata,20*np.log10(np.abs(Xk)))
ax[2].scatter(dph[part_idx,0]/(2*np.pi)*sr,20*np.log10(a[part_idx,0]))
print(th[part_idx,-1]-th[part_idx,0])
print(th_trans[part_idx,-1]-th_trans[part_idx,0])
plt.show()

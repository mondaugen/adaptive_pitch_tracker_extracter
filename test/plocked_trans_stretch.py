# Find peaks in a signal nearest some target peaks

import numpy as np
import matplotlib
matplotlib.use("Agg")
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
N_H=512
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
n_parts=48
max_detune_cents=200
# amount of transposition when resynthesizing
resynth_ptrans=np.power(2,0./12)
resynth_tstretch=5
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
res_=res#[res[20] for _ in res]

# synthesize results
y=np.zeros((np.round(n_max*resynth_tstretch).astype('int'),))
ph_last=res_[0][1] % 2*np.pi
n=0
_tstretch=resynth_tstretch
for (_,phk,ak,dphk,dak),(__,phk1,ak1,dphk1,dak1) in zip(res_[:-1],res_[1:]):
    phk1 = phk1 % 2*np.pi
    ph=cat_to_cols(ph_last,phk1)
    a=cat_to_cols(ak,ak1)
    dph=cat_to_cols(dphk,dphk1)
    da=cat_to_cols(dak,dak1)
    y_,a_last,ph_last,_,_,_=pa.synth_peaks_trans_stretch(ph,a,dph,da,
    ptrans=resynth_ptrans,
    tstretch=_tstretch)
    _tstretch=resynth_tstretch/(len(y_)/(pa.H*resynth_tstretch))
    _N=len(y_)
    y[n:n+_N]=y_
    n+=_N

out=np.concatenate((x,y))
out.tofile('/tmp/resynth.f64')

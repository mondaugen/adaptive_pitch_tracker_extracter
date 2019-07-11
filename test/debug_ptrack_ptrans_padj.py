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
N_H=2048
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
n_parts=50
part_num=5
max_detune_cents=200
# amount of transposition when resynthesizing
resynth_ptrans=np.power(2,7./12)
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
fig,ax=plt.subplots(6,1)
m=0
th_res=np.zeros(pa.H*(len(res)-1))
th_res_trans=np.zeros(pa.H*(len(res)-1))
y_res=np.zeros(pa.H*len(res))
y_res_trans=np.zeros(pa.H*len(res))
ph_last=res[0][1]
ph_last_trans=res[0][1]
for (Xk,phk,ak,dphk,dak),(Xk1,phk1,ak1,dphk1,dak1) in zip(res[:-1],res[1:]):
    ph=cat_to_cols(ph_last,phk1)
    ph_trans=cat_to_cols(ph_last_trans,phk1)
    a=cat_to_cols(ak,ak1)
    dph=cat_to_cols(dphk,dphk1)
    da=cat_to_cols(dak,dak1)
    n=np.arange(pa.H)
    y_res[n+m],ph_last,amp,th=pa.synth_peaks_trans(ph.copy(),a,dph,da)
    y_res_trans[n+m],ph_last_trans,amp_trans,th_trans=pa.synth_peaks_trans(ph_trans.copy(),a,dph,da,ptrans=resynth_ptrans)
    th_res[n+m]=th[part_num,:]
    th_res_trans[n+m]=th_trans[part_num,:]
    print((ph_last_trans[part_num]-ph[part_num,0])/(ph_last[part_num]-ph[part_num,0]))
    m+=pa.H

th_res=np.unwrap(th_res)
th_res_trans=np.unwrap(th_res_trans)

ax[0].plot(np.arange(len(th_res)),th_res,label='orig')
ax[0].plot(np.arange(len(th_res_trans)),th_res_trans,label='trans')
ax[0].legend()
ax[1].plot(np.arange(len(th_res)-1),np.diff(th_res),label='orig')
ax[1].plot(np.arange(len(th_res_trans)-1),np.diff(th_res_trans),label='trans')
ax[1].legend()

xdata=np.arange(N_FFT/2+1)*sr/N_FFT
ax[2].plot(xdata,20*np.log10(np.abs(Xk)))
ax[2].scatter(dph[part_num,0]/(2*np.pi)*sr,20*np.log10(a[part_num,0]))
err=x[:min(len(x),len(y_res))]-y_res[:min(len(x),len(y_res))]
print(max(x))
print(max(y_res))
ax[3].plot(np.arange(len(err)),x[:len(err)])
ax[4].plot(np.arange(len(err)),y_res[:len(err)])
ax[5].plot(np.arange(len(err)),err)
out=np.concatenate((y_res,y_res_trans,err))
out.tofile('/tmp/y_res.f64')
plt.show()

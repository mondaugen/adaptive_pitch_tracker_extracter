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

sr=16000
x=np.fromfile('/tmp/snd.f64')
x=normalize(x)

N_W=1024
N_FFT=2048
N_H=256
n_max=((len(x)-N_W)//N_H)*N_H

peak_thresh=-100

# intialize the sinusoidal analyser
pa=dft_freq_est.peak_analyzer(
N_FFT,
N_W,
N_H)

# range of peaks we are looking for
pitch_f0=36
f0=440*np.power(2,(pitch_f0-69)/12.)
n_parts=48
max_detune_cents=50
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
    ph,a,dph,da=pa.nearest_peaks(x_,w_ranges)
    #print('freqs',dph)
    #print('amps',a)
    res.append((X_,ph,a,dph,da))

# plot results as an animation
def plot_data_gen(t=0):
    cnt = 0
    while cnt < len(res):
        t = (cnt * N_H) / sr
        cnt += 1
        (X_,ph,a,dph,da)=res[cnt]
        yield (t, 20*np.log10(np.abs(X_)), 20*np.log10(a), dph/(2*np.pi)*sr)

fig, ax = plt.subplots()
ax.grid()
xdata=np.arange(N_FFT/2+1)*sr/N_FFT
ydata=[]

def init():
    pass

def run(data):
    # update the data
    t, X, a, f = data
    ax.clear()
    ax.set_ylim(-160, 3)
    ax.set_xlim(0, sr/2)
    ax.plot(xdata,X,label="%f" % (t,))
    ax.scatter(f,a)
    ax.legend()

fig.dpi=300
ani = animation.FuncAnimation(fig, run, plot_data_gen, blit=False, interval=10,
                              repeat=True, init_func=init)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15,codec='ayuv')
ani.save('/tmp/spec.mkv',writer=writer)
#plt.show()

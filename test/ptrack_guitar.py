import ptracker
import attack_finder
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import librosa
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap

def find_attacks(x,sr):
    xm=attack_finder.fast_max(np.abs(x),alph=0.9999)
    xmp=attack_finder.thresh_local_max_samples(np.diff(xm),alph=0.9999,beta=1,
    N=sr*.125)
    return xmp

x,sr=librosa.load('sounds/midi.wav',sr=16000)

y=attack_finder.analytic_signal(x)

t=np.arange(len(x))/sr

att=np.where(find_attacks(x,sr)>0)[0]

apt_args=(
0.9,
0.9,
0.9,
10,
0.9)
apta=ptracker.apt_array(apt_args,1)

start_n=att[0]-200
y_analysis=y[start_n:start_n+sr]
tracks=apta.proc(y_analysis,
32.703/sr*2*np.pi
#130.81/sr*2*np.pi,
#123.47/sr*2*np.pi,
)

fig,ax=plt.subplots(2,1)

NFFT=1024
h=4
noverlap=NFFT-h
ax[0].specgram(y_analysis, NFFT=1024, Fs=2*np.pi, Fc=0, noverlap=noverlap,
               cmap=None, xextent=None, pad_to=None, sides='default',
               scale_by_freq=None, mode='default', scale='default')

all_a=20*np.log10(np.concatenate(tuple([t[0] for t in tracks])))

out=np.zeros(y_analysis.shape,dtype="float")

for a,phi,corr,err in tracks:
    w=np.diff(np.unwrap(phi))
    t__=np.arange(len(phi))/sr
    t_=t__[:-1]
    #points = np.array([t_, w]).T.reshape(-1, 1, 2)
    #segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #norm = plt.Normalize(all_a.min(), all_a.max())
    #lc = LineCollection(segments, cmap='viridis', norm=norm)
    ## Set the values used for colormapping
    #lc.set_array(20*np.log10(a))
    #lc.set_linewidth(3)
    #line = ax[0].add_collection(lc)
    ax[0].plot(np.arange(len(w))/noverlap,w,'k')
    out+=a*np.cos(phi)
    ax[1].plot(t__,corr)
#fig.colorbar(line, ax=ax)

librosa.output.write_wav("/tmp/isitguitar.wav",out,sr)
librosa.output.write_wav("/tmp/orig.wav",np.real(y_analysis+np.conj(y_analysis)),sr)
plt.show()

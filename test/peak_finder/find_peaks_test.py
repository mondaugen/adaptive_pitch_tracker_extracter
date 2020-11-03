import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import os
from peak_finder import find_peaks

INFILE=os.environ.get('INFILE','/tmp/in.f32')
FS=float(os.environ.get('FS','16000'))
NFFT=int(os.environ.get('NFFT','2048'))
# frame to plot
C=int(os.environ.get('C','50'))
# peak threshold
K=float(os.environ.get('K',3))
x=np.fromfile(INFILE,'float32')

N=len(x)
n=np.arange(N)
t=n/FS

fig,ax=plt.subplots(2)
_,dummyax=plt.subplots(1)
S,f,t,_=dummyax.specgram(x,window=signal.get_window('blackmanharris',NFFT), NFFT=NFFT,noverlap=int(NFFT*.75),Fs=FS)
# normalize S
S/=S.max()
dummyax.plot([t[C],t[C]],[f.min(),f.max()])

ax[0].plot(f,20*np.log10(S[:,C]))
peaks,max_over_min,max_i,min_i,ext_val_fhold,ext_val_bhold=find_peaks(S[:,C],K=10,T=1e-5)
ax[0].plot(f[peaks],20*np.log10(S[peaks,C]),'.')
ax[1].plot(f,max_over_min)
#ax.plot(f[max_i],S[max_i,C],'.')
#ax.plot(f[min_i],S[min_i,C],'.')
#ax.plot(f,ext_val_fhold,label='fhold')
#ax.plot(f,ext_val_bhold,label='bhold')

plt.show()

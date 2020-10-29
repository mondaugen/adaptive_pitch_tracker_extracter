# Plot the waveform as a timeseries on the left and the spectrogram on the right
import matplotlib.pyplot as plt
import numpy as np
import os

INFILE=os.environ.get('INFILE','/tmp/in.f32')
FS=float(os.environ.get('FS','16000'))
x=np.fromfile(INFILE,'float32')

N=len(x)
n=np.arange(N)
t=n/FS

fig,ax=plt.subplots(1,2)

ax[0].plot(t,x)
ax[1].specgram(x,NFFT=2048,noverlap=2048-512,Fs=FS)
fig.suptitle(INFILE)

plt.show()

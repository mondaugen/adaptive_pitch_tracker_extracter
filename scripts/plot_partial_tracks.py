# Loads amplitude and phase tracks and plots:
#   - the tracks as instantaneous frequency vs. time
#   - through a window plotting amplitude versus frequency, the partial
#   amplitudes at a given time point selected in the tracks window.

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
import os
from common import wrap, dB

envget=os.environ.get

PTRACK_TH_A=envget('PTRACK_TH_A','/tmp/ptrack_th_a_out.npz')
AMPSCALE=envget('AMPSCALE','log')

with open(PTRACK_TH_A,'rb') as fd:
    npzfile=np.load(fd)
    files=sorted(npzfile.files)
    A,FS,H,Th=[npzfile[k] for k in files]

ph=np.angle(Th)
F=wrap(ph[:,1:]-ph[:,:-1])/(2.*np.pi)*FS
FMIN=F.min()
FMAX=F.max()
ylim=(FMIN,FMAX)
N=F.shape[-1]
t=np.arange(0,N,H)/FS
t_HFFT=H/FS
fig,axs=plt.subplots(3,1)
pt_ax=axs[0]
slider_ax=axs[1]
psd_ax=axs[2]
for f in F[:,0::H]:
    pt_ax.plot(t,f)

# Vertical line showing current frame shown in partial tracks window
tmp=pt_ax.get_ylim()
ylims=[tmp[0],tmp[0],tmp[1],tmp[1]]
# Store the line so we can move it with the slider
frame_line=pt_ax.axvline(0)

# Slider to choose a frame
sframe=Slider(slider_ax,'Frame (s)',0,t.max(),valinit=0,valstep=t_HFFT)
def update_frame_line(val):
    next_val=val+t_HFFT
    frame_line.set_xdata([val,val])
    frame=int(round(val*FS))
    psd_ax.clear()
    psd_ax.plot(F[:,frame],dB(A[:,frame]),'.')
    psd_ax.set_ylim([-200,dB(A.max())])
    # If specgram has log-scaled y axis, then PSD has log-scaled x-axis
    psd_ax.set_xlim(ylim)
    psd_ax.set_xscale(AMPSCALE)
update_frame_line(0)
sframe.on_changed(update_frame_line)

plt.show()




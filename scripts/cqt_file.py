import cqt
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import common

ge=common.get_env

filename=ge("FILENAME","/tmp/guit.f64",None)
min_pitch=ge("P_MIN",40,int)
max_pitch=ge("P_MAX",76,int)
pitch_res=ge("P_RES",1,float)
H=ge("H",256,int)
W=ge("W",1024,int)
window=ge("WINDOW","blackman",str)
p_range=np.arange(min_pitch,max_pitch+1,pitch_res)
sr=16000

x=np.fromfile(filename,dtype='float64')
C=cqt.cqt(
lambda N: signal.get_window(window,N),
W,
ws=[2*np.pi*440*np.power(2,(n-69)/12)/sr
    for n in p_range])
x=common.frame(x,H,W)
X=C(x)

ax=plt.imshow(20*np.log10(np.abs(X)),origin='lower',aspect='auto')
ax.axes.set_yticks([_ for _ in range(len(p_range))])
ax.axes.set_yticklabels(p_range)

plt.show()


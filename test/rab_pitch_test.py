import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import common
import rab_pitch

ge=common.get_env

filename=ge("FILENAME","/tmp/guit.f64",None)
min_pitch=ge("P_MIN",40,int)
max_pitch=ge("P_MAX",76,int)
pitch_res=ge("P_RES",1,float)
H=ge("H",256,int)
W=ge("W",1024,int)
N=ge("N",2048,int)
# variance of the gaussians
V=ge("V",1e-5,float)
window=ge("WINDOW","blackman",str)
p_range=np.arange(min_pitch,max_pitch+1,pitch_res)
sr=16000

def var_fun(f,n):
    return np.ones_like(n)*V#rab_pitch.default_variance(n)*V
def weight_fun(f,n):
    return np.power(1/(n+1),1)

x=np.fromfile(filename,dtype='float64')
C=rab_pitch.rab_pitch(
[440*np.power(2,(n-69)/12)/sr for n in p_range],
#variance=var_fun,
weight=weight_fun,
N_frame=N,
N_window=W)

C._plot_tables()
plt.show()

x=common.frame(x,H,N)
X=C(x)

ax=plt.imshow(np.log(X),origin='lower',aspect='auto')
ax.axes.set_yticks([_ for _ in range(len(p_range))])
ax.axes.set_yticklabels(p_range)

plt.show()


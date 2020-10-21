import gal_f32
import numpy as np
from scipy import signal
import filters
import matplotlib.pyplot as plt

N=200000
n=np.arange(N)
R=np.array([0.5,-0.5])
P=len(R)
a=filters.r_to_A(R)[-1,:]
x=np.random.standard_normal(N)
y=signal.lfilter([1],a,x)

mu=np.ones(P)*1e-5
gal=gal_f32.gal_f32(P)
galp=gal_f32.gal_f32_proc(y,mu,opt=(1<<0),beta=1e-3,l=1e-1)
gal.proc(galp)

fig,ax=plt.subplots(1,1)
for p in range(P):
    ax.plot(n,galp.R[:,p],label="r_%d" %(p,))
ax.legend()
ax.set_ylim(-1,1)
ax.set_xlim(0,N-1)
plt.show()

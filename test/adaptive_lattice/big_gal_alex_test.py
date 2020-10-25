import filters
import gal_alexander
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

N=10000
P=20
x=np.random.standard_normal(N)
K=np.array(np.random.uniform(-0.999,0.999,P))
a=filters.r_to_A(K)[-1,:]
y=signal.lfilter([1],a,x)
Ef,Eb,D,K=gal_alexander.ngal(y,P,alpha=0.001,beta=0.001,normalize=True)
n=np.arange(N)
fig,ax=plt.subplots(2,1)
for p in range(P):
    ax[1].plot(n,K[1:N+1,p],label="K_%d"%(p,))
ax[1].legend()
ax[0].plot(n,y)
plt.show()

